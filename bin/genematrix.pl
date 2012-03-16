#!/usr/bin/perl

use List::Util 'shuffle';
use Data::Dumper;
use Storable qw(dclone);
use strict;
# use diagnostics;
# use warnings;



die
"\nthis program takes 3 inputs:
1) a taxalist
2) a list of clustered genes, one cluster per line
3) a taxon id from the taxalist which should be the same or a subset of those in the cluster file

e.g. genematrix.pl taxalist-Bacillus BC_100000.1e-20_0.85_2_1.cluster_list BC_100000\n\n"

if ( scalar(@ARGV) != 3 );


#usage e.g. genematrix.pl taxalist-Bacillus BC_100000.1e-10_0.8_2_1.cluster_list BC_100000
# read in a file containing one cluster per line
# read in a list of refseq ids to analyze
# output matrix of shared clusters among ids (clusters x ids)
# output a matrix of phyletic patterns (cluster patterns x ids)

my $maxreps = 0;
$| = 1;

# read in the first ID and make sure it is of the proper format
my %taxahash;
my $taxalist = $ARGV[0];
open(IF, "<$taxalist") || die("$0 couldn't open file $taxalist\n");
while (<IF>) {
  chomp;
  my @taxa = split(/[ \t]/,$_);
  my $tag = shift @taxa;
  if (scalar(@taxa) > 0) {
    @{$taxahash{$tag}} = @taxa;
  } elsif (scalar(@taxa) == 0) {
  	next;
  } else {
    @{$taxahash{$tag}} = $tag;
  }
}
close(IF);


my @tmpids;
my $id;
if ($ARGV[2] =~ /\D{2}\_\d+/) {
    $id = $ARGV[2];
    die("$0 id $id not found\n") unless defined(@{$taxahash{$id}});
    @tmpids = @{$taxahash{$id}};
} else {
    die("failed to understand input argument 3: ", $ARGV[2], "\n");
}

# print Dumper(@tmpids),"\n";
# to de-code the list
my @ids = ();
my %allgenomes;
foreach my $moreid (@tmpids) {
  if (defined(@{$taxahash{$moreid}})) {
    push(@ids,[ @{$taxahash{$moreid}} ]);
    @allgenomes{@{$taxahash{$moreid}}} = (1) x @{$taxahash{$moreid}};
  }
  else {push(@ids,[ $moreid ]); $allgenomes{$moreid}++;}
}

my $preclustered;
if ( -e $ARGV[1] ) {
  $preclustered = $ARGV[1];
} else {die("cluster file not found\n")}

print "getting cogs...\n";
my @cogs = `/home/rec3141/repo/retrieve_sequences.pl cog 0 $preclustered $taxalist $id`;
my @cogmodes = &counter(@cogs);

print "getting functions...\n";
my @funcs = `/home/rec3141/repo/retrieve_sequences.pl func 0 $preclustered $taxalist $id`;
my @funcmodes = &counter(@funcs);

print "getting genes...\n";
my @genes = `/home/rec3141/repo/retrieve_sequences.pl gene 0 $preclustered $taxalist $id`;
my @genemodes = &counter(@genes);

print "getting descriptions...\n";
my @descs = `/home/rec3141/repo/retrieve_sequences.pl desc 0 $preclustered $taxalist $id`;
my @descmodes = &counter(@descs);

# we need a list of which genes are in each cluster (list of lists):
print "making hash...";

my @faadirs = qw|/work/rec3141/gene_clustering/all_genomes/faa /work/rec3141/gene_clustering/all_genomes/scaffold/faa /work/rec3141/gene_clustering/all_genomes/contig/faa /work/rec3141/gene_clustering/all_genomes/metagenomes/faa|;
my %accessionlist;
my %accseen=();
print "reading accession lists...";
foreach my $faadir (@faadirs) {
  open(IF,"<$faadir/accession_list") or die "couldn't open file $faadir/accession_list\n";
  while (<IF>) {
    my ($accgenome,$accgene) = split(/\t/,$_);
    if (exists($allgenomes{$accgenome})) {
      chomp $accgene;
      $accessionlist{$accgenome}{$accgene}++;
    }
  }
  close(IF);
}

#my @cluster_lists;

my %clusterhash;

# open the preclustered file and read into the hashes
my %genehash;
open(CLUST, "<$preclustered") or die "$0 couldn't open file $preclustered\n";
my(@preclusters) = <CLUST>;
#my $precluster=0;
close(CLUST);
#print "read file\n";
my $clusternumber=0;
foreach my $line (@preclusters) {
  $clusternumber++;
  chomp $line;
  my @proteins = split(/ /,$line);

  foreach my $k1 (@proteins) {
    print "DUPLICATE\n" if defined($genehash{$k1});
     $genehash{$k1}{'cluster'} = $clusternumber;
     $genehash{$k1}{'number'} = scalar(@proteins);
#     print "$k1\t$clusternumber\t$genehash{$k1}{'number'}\n";
  }
}
print "done.\n";

#read genes from each faa file into array
my %genelist;
my $gene;
my %singlelist;
my @newids;
my $filecount=0;
for my $i (0 .. $#ids) {
  my $newfileid = join('-',@{$ids[$i]});
  print "$i: $newfileid\n";
  push(@newids,$newfileid);

  #enables to add multiple molecules together e.g. plasmids + chromosomes
  foreach my $fileid (@{$ids[$i]}) {
    foreach my $gene (keys %{$accessionlist{$fileid}}) {
    push(@{$genelist{$newfileid}},$gene);
#    print "$gene\t$genehash{$gene}\n" if defined($genehash{$gene});
    if (defined($genehash{$gene}{'cluster'})) {
	  push(@{$clusterhash{$newfileid}},$genehash{$gene}{'cluster'});
	  push(@{$singlelist{$newfileid}},$gene) if ($genehash{$gene}{'number'} == 1);
	  delete($genehash{$gene});	  
    } else {
    # this should be corrected and not happen, if it does it's a bug
      print "UNCLUSTERED: $gene\t $genehash{$gene}{'cluster'}\n";
      push(@{$clusterhash{$newfileid}},'UNCLUSTERED');
      }
   }
  }

 if ((scalar grep { $_ != 0}  @{$clusterhash{$newfileid}}) == 0) {
  print "$newfileid was not clustered, removed from analysis\n";
  undef $newids[$filecount];
 }

$filecount++;
}

print join(' ',keys %genehash), " were excluded because they were not in $id\n" if (scalar(keys %genehash)>0);
@newids = grep { defined $_ } @newids;

my %orderhash;
foreach my $k1 (@newids) {
my @tmparray = ();
@tmparray = count_unique(@{$clusterhash{$k1}});
#print @tmparray;
@{$orderhash{$k1}} = @tmparray;

# print @tmparray;
# print $k1, "\t", join(' ',sort {$a <=> $b} @{$clusterhash{$k1}}), "\n" if defined(@{$clusterhash{$k1}});
#print $k1, "\t", join("\t",@{$orderhash{$k1}}), "\n";# if defined($orderhash{$k1});
}
# print Dumper(%clusterhash);
# print Dumper(%orderhash);

my @phyletics;
my @numphyletics;
open (PRETTY,">$preclustered.r.csv") or die "$0 couldn't open file $preclustered.r.csv\n";
print PRETTY "cluster\t# genomes\ttotal #\t",join("\t",@newids),"\tCOG\tCOG function\tgene\tdescription\n";
for (my $i=0;$i<$clusternumber;$i++) {
my $j = $i+1;
my @pattern = ();
my @abundance = ();
my $num=0;
my $numgenomes=0;
	foreach my $k1 (@newids) {
		if (defined(@{$orderhash{$k1}}[$i])) {
# 		  print @{$orderhash{$k1}}[$i],"\n";
		  push(@pattern,'1');
		  push(@abundance,@{$orderhash{$k1}}[$i]);
		  $num += @{$orderhash{$k1}}[$i];
		  $numgenomes++;
		} else {
		  push(@pattern,'0');
		  push(@abundance,'0');
		  $num += 0;
		  $numgenomes += 0;
		}
	}
# 	print join("\t",@abundance),"\n";
# 	print join("\t",@pattern),"\n";
	print PRETTY "$j\t$numgenomes\t$num\t",join("\t",@abundance),"\t", "$cogmodes[$i]\t$funcmodes[$i]\t$genemodes[$i]\t$descmodes[$i]\n";
	push(@phyletics,join("\t",@pattern));
	push(@numphyletics,join("\t",@abundance));
}
close(PRETTY);



sub count_unique {
    my @array = @_;
    my @list;
    my %count;
    map { $count{$_}++ } @array;

      #print them out:

#    map {print "$_ = $count{$_}\n"} sort {$a <=> $b} keys(%count);
# for every cluster list how many members it has ${count{$_}}
    map { $list[$_] = $count{$_} } sort {$a <=> $b} keys(%count);

      #or just return the hash:
# print @list;
# print Dumper(@list);
shift(@list); #because 'counts' start at 1 not 0
return @list;

#    return %count;
}

sub count_unique2 {
    my @array = @_;
    my @list;
    my %count;
    map { $count{$_}++ } @array;

      #print them out:

#PRINTME    map {print "$_ = ${count{$_}}\n"} sort {$count{$b} <=> $count{$a}} keys(%count);
#    map {$list[$_]=${count{$_}}} sort {$a <=> $b} keys(%count);

      #or just return the hash:
#return @list;
    return %count;

}


sub counter {

my @array = @_;
my @arraymodes;
for(0 .. $#array) {
 my $choice = 0;
 chomp $array[$_];
 my @send = split(/\t/,$array[$_]);
 my %count = ();
 foreach my $element (@send) {
    $count{$element}++;
  }
  my @mode = sort {$count{$b} <=> $count{$a}} keys(%count);
  $choice = 1 if ( ($mode[0] eq '-') && (defined($mode[1]) ) );
  $arraymodes[$_] = "[$count{$mode[$choice]}] $mode[$choice]";
#  print $arraymodes[$_],"\t";
}
return @arraymodes;
}



