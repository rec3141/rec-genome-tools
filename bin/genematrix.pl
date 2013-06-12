#!/usr/bin/perl
use strict;
use diagnostics;
use warnings;

use Getopt::Long;
use Pod::Usage;
use File::Slurp;

use Data::Dumper;
use RECGenomeTools qw(:All);
sub clusterCounter;

$|++; #flush the buffer

### get (required) options
my ($taxalist, $query_taxid, $cluster_list);

### get (optional) options
my $help = 0;
my $debug = 0;
my $output_dir = '.';
my $annotations = 'cog,cogfunc,gene,desc';
my $justcheckformissing = 0; #just check for missing proteins
my $maxproteins = 0; #passed to retrieve_sequences.pl to tell how many proteins to get details for; 0 for all

my @accdirs = qw|
		/work/rec3141/gene_clustering/all_genomes/faa
		/work/rec3141/gene_clustering/all_genomes/scaffold/faa
		/work/rec3141/gene_clustering/all_genomes/contig/faa
		|; #directories where accession lists live
		# /work/rec3141/gene_clustering/all_genomes/metagenomes/faa
my $hitsdir = '/work/rec3141/gene_clustering/prok'; # dir where hits files live

pod2usage(1) unless scalar(@ARGV);

GetOptions(
		#required
	   't|taxalist=s'  	=> \$taxalist,
	   'q|query=s'		=> \$query_taxid,
	   'l|cluster_list=s'	=> \$cluster_list,
	   #optional
	   'help|h' 		=> \$help,
	   'v|verbose+'  	=> \$debug,
	   'm|maxproteins'	=> \$maxproteins,
	   'a|annotations'	=> \$annotations,
	   'output_dir:s'	=> \$output_dir,
	   'accdirs:s'		=> \@accdirs,
	   'hitsdir:s'		=> \$hitsdir,
	   'justcheckformissing'	=> \$justcheckformissing,
	   );

pod2usage(-verbose => 2) if ($help);


### CHECK OPTIONS ###

### CLUSTER_LIST -- a cluster_list file
die("cluster file not found\n") unless (-e $cluster_list);

### TAXALIST -- a taxalist file
my %taxahash = %{ parse_taxalist( {'taxalist' => $taxalist} ) };
print scalar keys %taxahash, " taxahash keys found\n" if $debug;

### TAXALIST ID
### a taxalist id to limit clusters to
my ($ids_ref,$all_ids_ref) = parse_taxaid({'taxahash' => \%taxahash, 'id' => $query_taxid });
my @ids = @{ $ids_ref }; #list of lists of ids belonging to each genome
my %all_ids = %{ $all_ids_ref }; # $all_ids{'NC_123456'}
print scalar keys %all_ids, " taxaid keys found in $query_taxid\n" if $debug;

### ACCDIRS: locations of ACCESSION LISTs
@accdirs =  split(/,/,join(',',@accdirs));
foreach my $dir (@accdirs) {
	die("couldn't find directory $dir\n") unless ( -d $dir);
}

print "parsing accession_index\n" if $debug;
my ($accgenomes_ref, $genes_ref) = parse_accession_index( {'all_ids' => \%all_ids, 'parse_locus' => 1, 'accdirs' => \@accdirs} );


my %accgenomes = %{$accgenomes_ref}; # %accgenomes{'NC_010814'}{'YP_006429118.1'} = undef
print scalar keys %accgenomes, " genomes in accession_index for $query_taxid\n" if $debug;
my %genes = %{ $genes_ref }; # %accgenomes{'YP_006429118.1'} = 'NC_010814'
print scalar keys %genes, " genes in accession_index for $query_taxid\n" if $debug;

# open the cluster_list file and read into the hashes
my %genehash;
my @preclusters = read_file("$cluster_list");
my %clusterhash;

print "making gene hash\n" if $debug;
my $clusternumber=0;
foreach my $line (@preclusters) {
	chomp $line;
	my @proteins = split(/\s+/,$line);
	
	foreach my $k1 (@proteins) {
		print "DUPLICATE\n" if (defined($genehash{$k1}) and $debug>1);
		$genehash{$k1}{'cluster'} = $clusternumber;
		$genehash{$k1}{'single'} = undef if scalar(@proteins) == 1;
	}
		$clusternumber++;
}
print scalar keys %genehash," genes hashed\n" if $debug;

#read genes from each faa file into array
# my %genelist;
my %singlelist;
my $gene;
my @newids;
my $filecount=0;
for my $i (0 .. $#ids) {
	my $newfileid = join('-',@{$ids[$i]});
	print "$i: $newfileid\n" if $debug;
	push(@newids,$newfileid);

	#enables adding multiple molecules together e.g. plasmids + chromosomes
	foreach my $fileid (@{$ids[$i]}) {
		foreach my $gene (keys %{$accgenomes{$fileid}}) {
# 			push(@{$genelist{$newfileid}},$gene);
			if (defined($genehash{$gene}{'cluster'})) {
				push(@{$clusterhash{$newfileid}},$genehash{$gene}{'cluster'});
				push(@{$singlelist{$newfileid}},$gene) if exists($genehash{$gene}{'single'});
				delete($genehash{$gene});
			} else {
				# this can happen if some genes aren't clustered e.g. when 'rbh' clustered
				print "UNCLUSTERED: $gene\t $genehash{$gene}{'cluster'}\n" if $debug>1;
				push(@{$clusterhash{$newfileid}},-1); #UNCLUSTERED = -1
			}
		}
	}
	
	if ( (scalar grep { $_ != 0}  @{$clusterhash{$newfileid}}) == 0) {
		print "$newfileid was not clustered, removed from analysis\n" if $debug;
		undef $newids[$filecount];
	}
	
	$filecount++;
}

print scalar(keys %genehash), " genes were excluded because they were not in $query_taxid\n" if (scalar(keys %genehash)>0 and $debug);

exit 0 if $justcheckformissing;

undef %genehash;
undef %accgenomes;

@newids = grep { defined $_ } @newids;

# GET ANNOTATIONS
my $syscall = "cluster_info.pl $annotations $maxproteins $cluster_list 0 $taxalist $query_taxid";
print "getting cogs, cog functions, genes, and descriptions:\n$syscall\n" if $debug;
my @infos = qx/$syscall/;
print "$!" if $?;

my (@cogs,@cogfunc,@gene,@desc);
foreach (@infos) {
	chomp;
	my ($info_cog,$info_cogfunc,$info_gene,$info_desc) = split(/\x01/,$_);
	push(@cogs,$info_cog);
	push(@cogfunc,$info_cogfunc);
	push(@gene,$info_gene);
	push(@desc,$info_desc);
}

my @cogmodes = @{ clusterCounter(\@cogs) };
my @funcmodes = @{ clusterCounter(\@cogfunc) };
my @genemodes = @{ clusterCounter(\@gene) };
my @descmodes = @{ clusterCounter(\@desc) };

print scalar @infos," clusters annotated\n" if $debug;

#make a hash of lists of the number of proteins in each cluster in each organism
my %orderhash;
foreach my $k1 (@newids) {
	my $uniqhash_ref = count_unique(\%{$clusterhash{$k1}});
	my %uniqhash = %{$uniqhash_ref};
	foreach my $k2 (grep { $_ >= 0 } sort {$a <=> $b} keys %uniqhash) {
		${$orderhash{$k1}}[$k2] = $uniqhash{$k2};
	}
}

print "writing to $cluster_list.r.csv\n" if $debug;
open (PRETTY,">$cluster_list.r.csv") or die "$0 couldn't open file $cluster_list.r.csv\n";
print PRETTY "cluster\t# genomes\ttotal #\t",join("\t",@newids),"\tCOG\tCOG function\tgene\tdescription\n";
for my $i (0 .. $clusternumber-1) {
 	my $j = $i+1;
	my @pattern = ();
	my @abundance = ();
	my $num=0;
	my $numgenomes=0;
	foreach my $k1 (@newids) {
		if (defined(@{$orderhash{$k1}}[$i])) {
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

	print PRETTY "$j\t$numgenomes\t$num\t",join("\t",@abundance),"\t", "$cogmodes[$i]\t$funcmodes[$i]\t$genemodes[$i]\t$descmodes[$i]\n";
}
close(PRETTY);

print "Run Complete.\n" if $debug;

exit 0;



#####################
#### SUBROUTINES ####
#####################


# ------------------------
#  SUB COUNT_UNIQUE
# ------------------------
# sub count_unique {
#     my @array = @_;
#     my @list;
#     my %count;
#     map { $count{$_}++ } @array;
# 
#     map { $list[$_] = $count{$_} } sort {$a <=> $b} keys(%count);
# 
# 	shift(@list); #because 'counts' start at 1 not 0
# 	return @list;
# }

#------------------------
# SUB CLUSTERCOUNTER
#------------------------
sub clusterCounter {
	my $array_ref = shift;
	my @array = @{ $array_ref };

	my @arraymodes;
  
	foreach my $line (@array) {
		my $str = '[0|0] -'; #default string
		
		if ($line) {
			my @spl = grep defined, split(/\t/,$line);
			my %count = ();
			map { $count{$_}++ } @spl;
	
			my @mode = sort {$count{$b} <=> $count{$a}} keys(%count);
			my $choice=$mode[0];
			foreach (@mode) {
			  ($choice = $_) and last if ($_ !~ m/^(-|hypothetical)/);
			}
			$choice = 'unknown' if $choice eq '';
			$str = "[" . $count{$choice} . "|" . scalar(keys(%count)) . "] $choice";
		}
		push(@arraymodes,$str);
	}
	
  return(\@arraymodes);
}




##############################
####### END           ########
##############################


 __END__
 
=head1 NAME

genematrix - annotates gene cluster

=head1 SYNOPSIS

This program converts a gene cluster list into a spreadsheet containing
gene copy numbers and the following [default] annotations: cog, cogfunc, gene, and desc
obtained by calling cluster_info.pl

=item B<example:>

	genematrix.pl -t taxalist-Bacillus -q BL_110001 -l BL_110001.1e-20_0.7_2_rbh_link.cluster_list


=head2 Required:

-t|taxalist: a file containing the refseq identifier

-q|query: query taxa, a refseq identifier or group id from a taxalist

-l|cluster list: a file containing gene identifiers, one cluster per line

=head1 OPTIONAL ARGUMENTS

-help|h                 print help

-v|verbose              zero for silent, once for normal, more for debugging



-accdirs                directories where accession_index files live

-hitsdir                directory to find blast hit files (.hits)

-justcheckformissing    short circuit and just look for missing blast hits

-output_dir             base directory for output files [./]

=head1 OPTION DETAILS


=head1 OUTPUTS

=item B<*.r.csv>
a gene matrix suitable for importation into R or a spreadsheet program
with the following properties:

 1) matrix containing frequency of protein sequences in clusters (rows) in each genome (columns)

 2) modal annotations for each cluster including cog, cog function, gene label, and gene description
   key to annotations: '[a|b]' where a=number of genes with modal (most frequent) annotation, b=number of different annotations

=head1 RESOURCE UTILIZATION
a cluster list containing 520 genomes and 72218 clusters was converted
to a gene matrix in 10 minutes on 1 CPU using 4Gb of memory

=head1 CHANGES

=item B<2/04/2013:> reformat

=cut
