#!/usr/bin/perl
use lib "$ENV{HOME}/lib/perl5";
use File::Slurp;
use Data::Dumper;

use strict;
use warnings;
use diagnostics;

my $debug=0;
$| = 1; #flush

my @faadirs = qw|/work/rec3141/gene_clustering/all_genomes/faa /work/rec3141/gene_clustering/all_genomes/scaffold/faa /work/rec3141/gene_clustering/all_genomes/contig/faa /work/rec3141/gene_clustering/all_genomes/metagenomes/faa|;
my $prokdir = '/work/rec3141/gene_clustering/prok';
my $justcheckformissing; #define to avoid clustering and just find missing BLAST hits

# Summary:

# This script builds a cluster map of a species genes based on 
# single-link clusters.
#
# It produces a statistical summary that tabulates how many clusters 
# of each size there are in *.cluster_hist
#
# and a file that lists which genes are in  which cluster, *.cluster_list

# this sets the largest bin for the histogram. 
my $max_genes_cluster=1000;

# which output files to produce
my $output_hist = 1;
my $output_list = 1;
my $output_graph = 0;

# check to make sure that the program was given the correct number
# of input args 

if (scalar(@ARGV) != 6) {
  die("\n\nscript takes:
  
  a taxalist file containing the refseq identifier\n
  a refseq identifier or group id from a taxalist
  an expectation cutoff
  a matching fraction cutoff
  a flag (long, short, seed) to determine which gene to consider for the length of match \n
    'seed' tells whether to allow good clusters to recruit shorter fragments on 2nd pass
  a flag (off or on) to require reciprocal hits (but not reciprocal BEST hits)\n
    if 'seed', recruited fragments do not require reciprocal hits, others do

  eg. cluster.pl taxalist-Bacillus NC_000123 1e-20 0.8 seed on

  CHANGED: allow a second pass that allows short matches to join clusters
	    but does not allow clusters to merge\n\n");
}


### TAXALIST
###  a taxalist file containing the refseq identifier
my $taxalist = shift(@ARGV);
die("failed to understand input argument 1: ",$taxalist, "\n") unless ( -e $taxalist);

my %taxahash;
open(IF, "< $taxalist") || die("couldn't open file $taxalist\n");
while (<IF>) {
  chomp;
  my @taxa = split('\s',$_);
  my $tag = shift @taxa;
  next unless $tag;
  if (scalar(@taxa) > 0) {
    @{$taxahash{$tag}} = @taxa;
  } else {
    @{$taxahash{$tag}} = $tag;
  }
}
close(IF);

### TAXA ID
###   a refseq identifier or group id from a taxalist
my @ids;
my $id =  shift(@ARGV);
if ($id =~ /^\D{2}\_\d+$/) {
    if (defined(@{$taxahash{$id}})) {
      # to de-code the list
      foreach my $moreid (@{$taxahash{$id}}) {
	if (defined(@{$taxahash{$moreid}})) {
	  push(@ids,[ @{$taxahash{$moreid}} ]);
	} else {push(@ids,[ $moreid ]);}
      }
    } else {warn("$id not found\n")}
} else {
    die("failed to understand input argument 2: ", $id, "\n");
}

### E-CUTOFF
###   an expectation cutoff
my $e_cutoff = shift(@ARGV);
die("failed to understand input argument 3: ",$e_cutoff, "\n") unless ($e_cutoff =~ /^\d.+\.?e?\d.*$/);

### LENGTH FRACTION
###   a matching fraction cutoff
###   a flag (long or short) to determine which gene to consider for the length of match \n
# read in the minimum length fraction
# 0 no restriction
# 1 by shorter sequence
# 2 by longer sequence
# 3 allow seeding of good fragments

my $mlf = shift(@ARGV);
my $mlf_name = shift(@ARGV);
my $mlf_type;
my $seeded;
die("failed to understand input argument 4: ",$mlf, "\n") unless ($mlf =~ /^\d*\.?\d*$/);
while ($mlf >= 1) {$mlf=$mlf/10};
if ($mlf_name =~ /seed/) {
  $seeded = 1;
  $mlf_type = 2;
}
$mlf_type = 2 if ($mlf_name =~ /long/);
$mlf_type = 1 if ($mlf_name =~ /short/);
$mlf_type = 0 if ($mlf == 0.0);
die("failed to understand input argument 5: ",$mlf_name, "\n") unless defined($mlf_type);

### RECIPROCAL
###   a flag (off or on) to require reciprocal hits (but not reciprocal BEST hits)\n
my $recip_name = shift(@ARGV);
my $recip;
$recip=0 if ($recip_name =~ /off/);
$recip=1 if ($recip_name =~ /on/);
die("failed to understand input argument 6: ",$recip_name, "\n") unless defined($recip);


print "Parameters for this run: $e_cutoff\t$mlf\t$mlf_name\n";

########
#########
########## BEGIN
my $starttime=time;

### genome lists for counting
my @allgenomes;
my %allgenomes;
my $ngenomes = scalar(@ids);
for (my $i=0;$i<$ngenomes;$i++) {
 foreach my $genome (@{$ids[$i]}) {
 push(@allgenomes,$genome);
 $allgenomes{$genome}++;
 }
}
my $nreplicons = scalar(@allgenomes);

my %genelist;
my %seen=();
my %genehash;
print "reading gene list...";
foreach my $faadir (@faadirs) {
open(IF,"<$faadir/accession_list") or die "couldn't open file $faadir/accession_list\n";
while (<IF>) {
 my ($genome,$gene) = split(/\t/,$_);
 if (exists($allgenomes{$genome})) {
  chomp $gene;
  $genelist{$gene} = undef;
  $seen{$genome}++;
  $genehash{$genome}{$gene}=undef;
 }
}
close(IF);
}
my $genesum;
foreach (@allgenomes) {
  die "genome $_ missing from accession_list\n" unless $seen{$_};
  $genesum += $seen{$_};
}
print "found " . scalar(keys %genelist) . " genes in " . (time - $starttime) . "s\n";


# open the significant hits files that were generated by 
# hit_extractor.pl and read each line into a string that is 
# stored in an array 
# 
# at present each line is of the form:
# refseq_query_id refseq_hit_id significance mlf_short mlf_long
#
# eg.
# NP_012345.1 NP_6789.2 1.35e-32 0.23 0.11
#
# the filename is ordered in the opposite direction though,
# Genome_Database.Genome_Queries.hits

my %hoh1; #cluster hash for first pass
my %hoh2; #cluster hash for second pass
$starttime = time;
my %orphaned; #genomes with outdated accession_list
my %redo; #missing hits files
my $i=0; # genome counter 
print "reading files...\n";
my $missing = 0; #missing hits files?
foreach my $id1 (@allgenomes) {
  printf("%s\t %d of %d",$id1, ++$i,$nreplicons);

  foreach my $id2 (@allgenomes) {
    my @lines = ();
    if (-e "./$id1/$id1.$id2.hits") {
      next if ($justcheckformissing); #short circuit
      @lines = read_file("./$id1/$id1.$id2.hits");
    } else {
      print "-";
      $redo{$id1}{$id2}++;
      $missing++;
      next;
    }

    foreach my $hit (@lines) {
      #id2-protein $id1-protein evalue shorter-frag-frac longer-frag-frac
      chomp $hit;
      my ($query_id,$hit_id,$e_val,$mlfs,$mlfl)=split(/\s+/,$hit); #split the hit
      warn("$0: problem in hit line: $hit\n") and next unless (defined($query_id) && defined($hit_id) && defined($e_val) && defined($mlfs) && defined($mlfl));
      
      # find if proteins not in accession_list
      $orphaned{$id1}=undef unless(exists($genelist{$hit_id}));
      $orphaned{$id2}=undef unless(exists($genelist{$query_id}));

      # do NOT save self hits to hash
      # will add these back in later if they have no other hits
      next if ($query_id eq $hit_id);
      next if ( $e_val > $e_cutoff); #skip bad matches (except self hits)

      delete($genehash{$id2}{$query_id});
      delete($genehash{$id1}{$hit_id});

#       my ($ha,$hb) = sort { $a cmp $b} ($query_id,$hit_id);
# 
#       # build a hash of hashes with good matches
#       if (( $mlf_type == 2 ) && ( $mlfl > $mlf )) {
# 	$hoh1{ $ha }{ $hb }++; #if longer is longer than cutoff
#       } elsif (defined($seeded) && ( $mlfs > $mlf )) {
# 	if (exists($hoh2{ $ha }{ $hb })) {
# 	  $hoh2{ $ha }{ $hb } = $e_val if $e_val < $hoh2{ $ha }{ $hb };
# 	} else {
# 	  $hoh2{ $ha }{ $hb } = $e_val; #if longer is not long enough but shorter is
# 	}
#       } elsif (( $mlf_type == 1 ) && ( $mlfs > $mlf )) {
# 	$hoh1{ $ha }{ $hb }++; #if shorter is longer than cutoff
#       } elsif ( $mlf_type == 0 ) {
# 	$hoh1{ $ha }{ $hb }++;  #no length restriction
#       } else {print "$mlfl < $mlfs < $mlf ?\n" if $debug==1;}
    } #end foreach hit
    exists($redo{$id1}{$id2}) ? print "+" : print ".";
  } #end foreach file2
  
  print "done hashing in ". (time-$starttime) . "s\n";
  #could cluster here to save memory on %hoh1 if not reciprocals...?
} #end foreach file1

print "outdated in accession_list: ", join("\t",keys %orphaned), "\n" if (keys %orphaned > 0);
undef %orphaned;

if ($missing>0) {
  open(OF, ">>./$id/$id-blast-redo.run") or die "couldn't open ./$id/$id-blast-redo.run";
    foreach my $k (sort keys %redo) {
      foreach my $kk (sort keys %{$redo{$k}}) {
	  print OF "rm $prokdir/$k/$k.$kk.hits; /home/rec3141/repo/signifier.sh $k $kk\n";
      }
    }
  close(OF);
  print "\n\nRedo the blasts in $id-blast-redo.run\n";

  exit 0 if ($justcheckformissing);
}
undef %redo;


# print Dumper %genehash;
foreach my $key (keys %genehash) {
  print $key,"\t";
  print join("\t",keys %{$genehash{$key}});
  print "\n";
}
exit 0;

# deal with reciprocals
## loop over queries in hash 1 to find reciprocal hits
## remove those hash entries that don't have reciprocals
if ( $recip == 1 ) {
$starttime = time;
print "finding reciprocal hits...";

 foreach my $k1 ( keys %hoh1 ) {
  foreach my $k2 ( keys %{$hoh1{$k1}} ) {
   delete($hoh1{$k1}{$k2}) if ( $hoh1{$k1}{$k2} < 2);
  }
  delete($hoh1{$k1}) unless (keys %{$hoh1{$k1}});
 }
print "done in " . (time - $starttime) . "s\n";
}


# Begin Clustering
print "clustering...";
$starttime = time;
# cluster data
# we need an index for how many clusters have been created:
my $num_clus=0;
# we need a list of which genes are in each cluster (list of lists):
my @cluster_lists;
# we need a hash of which cluster each gene is associated with:
my %gene_clus; #   %gene_clus{'NC_111000'} = 1
# a hash of which genes are in each cluster
my %clus_gene; #   %clus_gene{'1'}{'NC_111000'} = undef, %clus_gene{'1'}{'NC_111111'} = undef

# now we need to loop through the hits and cluster them.  All the hits
# left in the hoh1 hash structure are valid at this point and need to
# be catalogued.

#get each set of blast hits $k1 and $k2
while (my ($k1, undef) = each %hoh1) {
  foreach my $k2 ( keys %{$hoh1{$k1}} ) {
    # if $k1 is already in a cluster
    if ( defined $gene_clus{$k1} ) {
      # and if $k2 is already in a cluster
      if ( defined $gene_clus{$k2} ) {
      # then they need to be merged
      # don't do anything if they're already in the same cluster
        unless ( $gene_clus{$k1} == $gene_clus{$k2} ) {
	  my $todel = $gene_clus{$k2};
          #and add all of the genes in the second cluster to the first
	  foreach my $k3 (keys %{$clus_gene{ $gene_clus{$k2} }}) {
	  $gene_clus{$k3} = $gene_clus{$k1};
	  $clus_gene{$gene_clus{$k1}}{$k3} = undef;
	  }
	  #then delete second cluster
	  delete($clus_gene{ $todel });
        }
        # if gene2 isn't already in a cluster, add it to that of gene1
      } else {
        $gene_clus{ $k2 } = $gene_clus{ $k1 };
	$clus_gene{$gene_clus{$k1} }{ $k2 } = undef;
      }
      # if gene1 isn't in a cluster but gene2 is, add gene1 to it
    } elsif ( defined $gene_clus{$k2} ) {
      $gene_clus{ $k1 } = $gene_clus{ $k2 };
      $clus_gene{ $gene_clus{$k2} }{ $k1 } = undef;
      #if neither gene1 nor gene2 are in a cluster, add them both to a new cluster
    } else {
      $gene_clus{$k1} = $num_clus;
      $gene_clus{$k2} = $num_clus;
      $clus_gene{$num_clus}{$k1} = undef;
      $clus_gene{$num_clus}{$k2} = undef;
      $num_clus++;
    }
    delete($genelist{$k2});
  }
  delete($genelist{$k1});
  delete($hoh1{$k1});
}
print "found " . scalar(keys %clus_gene) . " clusters\n";
undef %hoh1;

print Dumper(%hoh2);
### now for second pass
if (defined $seeded) {
my $seeds;
print "second pass...";
 foreach my $k1 ( keys %hoh2 ) {
  foreach my $k2 ( keys %{$hoh2{$k1}} ) {
   if ( $gene_clus{$k1} ) {
    if ( $gene_clus{$k2} ) {
      # if both are already in clusters, delete them
      delete($hoh2{$k1}{$k2});
    } else {
      # if second has match in a cluster, move it to first position
      # so it becomes $hoh2{fragment}{match}
      $hoh2{$k2}{$k1} = $hoh2{$k1}{$k2};
      $seeds++;
      delete $hoh2{$k1}{$k2};
    }
   } elsif ( $gene_clus{$k2} ) {
      undef; # if first has match in cluster, keep it 
   } else {
    delete $hoh2{$k1}{$k2}; #if neither have match in cluster, delete 
   }
  }
  delete($hoh2{$k1}) unless (keys %{$hoh2{$k1}});
 }

### after filtering, find best hit for each remaining protein
### and add protein to cluster of best hit

 foreach my $k1 ( keys %hoh2 ) {
  foreach my $k2 ( sort { ($hoh2{$k1}{$a} || 1) <=> ($hoh2{$k1}{$b} || 1) } keys %{$hoh2{$k1}} ) {
    $gene_clus{$k1} = $gene_clus{$k2};
    $clus_gene{ $gene_clus{$k2} }{$k1} = undef;
    delete($genelist{$k2});
    last;
  }
  delete($genelist{$k1});
 }
print "$seeds added\n";
}
print Dumper %hoh2;
undef %hoh2;

# add the unclustered genes/singletons to cluster_lists
foreach my $k1 ( sort {$a cmp $b} keys %genelist ) {
 push(@cluster_lists, ( [$k1] ));
}
undef %genelist;

print "done in " . (time - $starttime) . "s\n";
$starttime = time;
print "making cluster list...";

#make @cluster_lists
map {
  push(
    @cluster_lists, 
    [ sort {$a cmp $b} keys %{ $clus_gene{$_} } ]
  )
} keys %clus_gene;

#first sort the cluster list by order of genes
@cluster_lists = sort {"@$a" cmp "@$b"} @cluster_lists;
#then sort the cluster list by number of genes in the cluster
@cluster_lists = sort {@{$b} <=> @{$a}} @cluster_lists;

print "done in " . (time - $starttime). "s\n";

if ( $debug == 2 ) {
  print "NUM_CLUSTERS:",scalar(@cluster_lists),"\n";
  # We've built the clusters!  Time to analyze:
  foreach my $k1 ( @cluster_lists ) {
    if ( defined @{$k1} ) {
      print "CLUSTER: @{$k1}\n";
    } else {
      print "CLUSTER: UNDEFINED\n";
    }
  }
}

$starttime = time;
# find genes in master gene list but not in clustered list
# these genes didn't have any good hits or
# these genes were excluded by the blastall or formatdb programs
# they tend to be short and may be gene fragments
# unhashed => singletons

mkdir($id) unless (-d $id);

# OUTPUT CLUSTER LIST
if ($output_list == 1) {
open(OF, ">./$id/$id.$e_cutoff\_$mlf\_$mlf_type\_$recip.cluster_list") or die "couldn't open ./$id/$id.$e_cutoff\_$mlf\_$mlf_type\_$recip.cluster_list";
for my $cluster ( @cluster_lists ) {
 print OF "@$cluster\n"; #join(' ', sort {$a cmp $b} @{$cluster}),"\n";
}
close(OF);
print "output ./$id/$id.$e_cutoff\_$mlf\_$mlf_type\_$recip.cluster_list\n";
}



# MAKE HISTOGRAM
my @hist;
$num_clus=scalar(@cluster_lists);
my $num_genes=0;
my $max_genes_found=0;
foreach my $k1 ( @cluster_lists ) {
    $num_genes=scalar(@{$k1});
    $hist[$num_genes]++;
}

for (my $i=0;$i<scalar(@hist);$i++) {
  $hist[$i] = 0 unless defined($hist[$i]);
}

print "Number of clusters = ". ($num_clus - $hist[1]) . "\n";
my $clustsum; 
for(my $i=2;$i<=$max_genes_cluster;$i++) {$hist[$i] ? $clustsum+=$hist[$i]*$i : $clustsum+=0};
print "Number of genes in clusters = $clustsum\n"; 
print "Number of singletons = $hist[1]\n";
print "Total genes = " . ($hist[1]+$clustsum) . "\n";

# OUTPUT HISTOGRAM
if ($output_hist == 1) {
open(OF, ">>./$id/$id.cluster_hist") or die "couldn't open ./$id/$id.cluster_hist";
#print OF "genes\tcluster\n";
my $total_clus=0;
my $over_clus=0;
print OF "$id\t$e_cutoff\t$mlf\t$mlf_name\t$recip_name\t" , $seeded ? "seeded\t" : "noseed\t";
for (my $ii = 1; $ii < scalar(@hist); $ii++) {
  if ($ii > $max_genes_cluster) {
    $over_clus += $hist[$ii];
    next unless ($ii == scalar(@hist)-1);
    print OF "\t$over_clus";
    $total_clus=$total_clus+$hist[$ii];
  } elsif (defined($hist[$ii])) {
    print OF "\t$hist[$ii]";
    $total_clus=$total_clus+$hist[$ii];
  } else {print OF "\t0";}
}
print OF "\n";
close(OF);
warn("Number of clusters incongruent $total_clus $num_clus\n") if ($total_clus != $num_clus);
print "output ./$id/$id.cluster_hist\n";
}

exit 0