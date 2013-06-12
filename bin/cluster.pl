#!/usr/bin/perl
#use lib "$ENV{HOME}/lib/perl5";
use File::Slurp;
use Data::Dumper;
use List::Util 'shuffle';

sub print_redos;
sub DoClust;
sub procMCL;
sub finishMCL;
sub procHITS;
sub procSeeds;
sub procRBH;

use strict;
use warnings;
use diagnostics;

my $debug=0;
$| = 1; #flush

# my $scratchdir = '/scratch/rec3141/';
my @faadirs = qw|/work/rec3141/gene_clustering/all_genomes/faa /work/rec3141/gene_clustering/all_genomes/scaffold/faa /work/rec3141/gene_clustering/all_genomes/contig/faa|; 
# /work/rec3141/gene_clustering/all_genomes/metagenomes/faa|;
my $prokdir = '/work/rec3141/gene_clustering/prok';
my $justcheckformissing = 0; #define to avoid clustering and just find missing BLAST .hits files
my $justcheckfororphaned = 0; #define to check for BLAST hits not in accession_index, indicating an obsolete protein
my $authorize_redo = 0;# allow in-line redo of missing or outdated BLASTs

Summary:

# This script builds a cluster map of a species genes based on 
# single-link clusters.
#
# It produces a statistical summary that tabulates how many clusters 
# of each size there are in *.cluster_hist
#
# and a file that lists which genes are in  which cluster, *.cluster_list

# this sets the largest bin for the histogram. 
my $max_genes_cluster=10000;

# which output files to produce
my $output_hist = 1;
my $output_list = 1;
my $output_graph = 0;

# check to make sure that the program was given the correct number
# of input args 

if (scalar(@ARGV) < 7) {
  die("\n\nscript takes:
  
  1) a taxalist file containing the refseq identifier\n
  2) a refseq identifier or group id from a taxalist
  3) an expectation cutoff
  4) a matching length fraction cutoff
  5) a flag (long, short, seed) to determine which gene to consider for the length of match \n
    'none' turns off match length cutoff (overrides #4 and sets cutoff to 0)
    'long' requires #4 to match longer sequence
    'short' requires #4 to match shorter sequence 
    'seed' tells whether to allow good clusters to recruit shorter fragments on 2nd pass
       (implies 'long')
  6) a flag to require reciprocoal hits
      'on' to require reciprocal hits
        (if #5='seed', recruited fragments do not require reciprocal hits, others do)
      'off' to NOT require reciprocal hits     
      'rbh' to use only Reciprocal Best/Top Hits (RBHs)
  7) a flag to determine the type of clustering
      'link' to perform SINGLE-LINKAGE clustering
      'mcl' to perform MCL clustering using the Markov Clustering algorithm
  8) REQUIRED iff #7 is mcl: an inflation parameter for the MCL algorithm.
       Default is 2.0. MCL author suggests trying 1.4, 2, 4, and 6

  eg. cluster.pl taxalist-Bacillus NC_000123 1e-20 0.8 seed on link

  CHANGED:  1/20/2013: enable RBH filtering
	    3/30/2012: enable MCL clustering
	    2/01/2012: allow a second pass that allows short matches to join clustered seeds\n\n");
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
if ($mlf_name =~ /seed/i) {
  $seeded = 1;
  $mlf_type = 2;
}
$mlf_type = 2 if ($mlf_name =~ /long/i);
$mlf_type = 1 if ($mlf_name =~ /short/i);
$mlf_type = 0 if ($mlf == 0.0 or $mlf_name =~ /none/i);
die("failed to understand input argument 5: ",$mlf_name, "\n") unless defined($mlf_type);

### RECIPROCAL
###   a flag (off or on) to require reciprocal hits (but not reciprocal BEST hits)\n
my $recip_name = shift(@ARGV);
my $recip;
my $do_rbh = 0;
$recip=0 if ($recip_name =~ /off/i);
$recip=1 if ($recip_name =~ /on/i);
$recip=2 if ($recip_name =~ /rbh/i);
die("failed to understand input argument 6: ",$recip_name, "\n") unless defined($recip);

### CLUSTERING
###    a flag to determine the type of clustering
my $clusttype_name = shift(@ARGV);
my $mcl;
$mcl=0 if ($clusttype_name =~ /link/i);
$mcl=1 if ($clusttype_name =~ /mcl/i);
die("failed to understand input argument 7: ",$clusttype_name,"\n") unless defined($mcl);


my $inflation_name;
my $inflation;
if ($mcl>0) {
  $inflation_name = shift(@ARGV);
  $inflation = undef;
  $inflation = 0+$inflation_name if ($inflation_name =~ /[0-9.]+/);
  die("failed to understand input argument 8: ",$inflation_name,"\n") unless defined($inflation);
} else{$inflation=0}


if ($mcl>0) {
	print "Parameters for this run: $e_cutoff\t$mlf\t$mlf_name\t$recip_name\t$clusttype_name\t$inflation\n";
} else {
	print "Parameters for this run: $e_cutoff\t$mlf\t$mlf_name\t$recip_name\t$clusttype_name\n";
}

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
print "reading gene list...";
foreach my $faadir (@faadirs) {
  open(IF,"<$faadir/accession_index") or die "couldn't open file $faadir/accession_index\n";
  while (<IF>) {
    chomp;
    my @line = split(/\s+/,$_);
    my $genome = shift(@line);
    next unless exists($allgenomes{$genome});
    foreach my $gene (@line) {
      $genelist{$gene} = 1;
      $seen{$genome}++;
    }
  }
  close(IF);
}

my $genesum;
foreach (@allgenomes) {
  print "genome $_ missing from accession_index\n" unless $seen{$_};
  $genesum += $seen{$_};
}
print "found " . scalar(keys %genelist) . " genes in " . (time - $starttime) . "s\n";


# open the significant hits files that were generated by 
# blast_extract.pl and read each line into a string that is 
# stored in an array 
# 
# at present each line is of the form:
# refseq_query_id refseq_hit_id significance mlf_short mlf_long
#
# eg.
# NP_012345.1 NP_678901.2 1.35e-32 0.23 0.11
#
# the filename is ordered in the opposite direction though,
# Genome_Database.Genome_Queries.hits

# to do in-progress clustering
# cluster data
# we need an index for how many clusters have been created:
my $num_clus=1; #don't start at 0 because we check for true/false
# we need a hash of which cluster each gene is associated with:
my %gene_clus; #   %gene_clus{'NC_111000'} = 1
# a hash of which genes are in each cluster
my %clus_gene; #   %clus_gene{'1'}{'NC_111000'} = undef, %clus_gene{'1'}{'NC_111111'} = undef


if($mcl>0) {
  open MCL, "| mcl - --abc --abc-neg-log -I $inflation -o ./$id/$id.$e_cutoff\_$mlf\_$mlf_type\_$recip\_$clusttype_name\_$inflation.cluster_list" or die "can't open output: $!";
}

my %hoh1; #cluster hash for first pass %hoh1{'NC_000000'}{'NC_111111'} = number of times seen (1 or 2)
my %hoh2; #cluster hash for second pass %hoh2{'NC_000000'}{'NC_111111'} = 1e-30
my %redo; #missing hits files and genomes with outdated accession_list
my $i=0; # genome counter 
my $missing = 0; #missing hits files?

print "reading files and clustering...\n";
$starttime = time;

@allgenomes = shuffle(@allgenomes);

foreach my $id1 (@allgenomes) {
	printf("%s\t %d of %d",$id1, ++$i,$nreplicons);
	FILE2: foreach my $id2 (@allgenomes) {
		if (-e "$prokdir/$id1/$id1.$id2.hits") {
		  next if ($justcheckformissing); #short circuit
		  my @list = read_file("$prokdir/$id1/$id1.$id2.hits");
		  my $newlist_ref;
		  $recip > 1 ? $newlist_ref = procRBH( \@list ) : $newlist_ref = \@list;
		  procHITS( ($id1,$id2,$newlist_ref) ); #process hits file
		} else {
		  print "-";
		  $redo{$id1}{$id2}++;
		  $missing++;
		  next;
		}
		
		if (exists($redo{$id1}{$id2})) {
			print $redo{$id1}{$id2};
			&print_redos;
			%redo=();
		} else {
				print ".";
		}
	} #end foreach file2
  
	unless(($justcheckformissing || $justcheckfororphaned)) {
	  print "done hashing at ". (time-$starttime) . "s, ";
	  &DoCluster;
	  print "done clustering at ". (time-$starttime) . "s";
	}
	print "\n";
} #end foreach file1

print "\n\nRedo the blasts in $id-blast-redo.run\n" if (($missing>0) && !($authorize_redo));
exit 0 if ($justcheckformissing || $justcheckfororphaned);

# mcl will take over from here if MCL
&finishMCL if ($mcl>0);


# Finish Clustering
$starttime = time;
print "found " . scalar(keys %clus_gene) . " clusters\n";

#what's left in %hoh1 at this point?
#non reciprocal hits if $recip==1
#are these represented in %hoh2? yes, the reciprocals are in there
# print Dumper(%hoh1);
undef %hoh1;

### now for second pass
&procSeeds if (defined $seeded);

### get ready for printing
$starttime = time;
print "making cluster list...";

# we need a list of which genes are in each cluster (list of lists):
my @cluster_lists;

#make @cluster_lists
map {
  push(
    @cluster_lists, 
    [ sort {$a cmp $b} keys %{ $clus_gene{$_} } ]
  )
} keys %clus_gene;
undef %clus_gene;

# add the unclustered genes/singletons to @cluster_lists
unless ($recip>1) {
  foreach my $k1 ( sort {$a cmp $b} keys %genelist ) {
  push(@cluster_lists, ( [$k1] )) if defined($genelist{$k1});
  }
}
undef %genelist;

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


mkdir($id) unless (-d $id);

# OUTPUT CLUSTER LIST
if ($output_list == 1) {
open(OF, ">./$id/$id.$e_cutoff\_$mlf\_$mlf_type\_$recip\_$clusttype_name.cluster_list") or die "couldn't open ./$id/$id.$e_cutoff\_$mlf\_$mlf_type\_$recip\_$clusttype_name.cluster_list";
for my $cluster ( @cluster_lists ) {
 print OF "@$cluster\n"; #join(' ', sort {$a cmp $b} @{$cluster}),"\n";
}
close(OF);
  print "output ./$id/$id.$e_cutoff\_$mlf\_$mlf_type\_$recip\_$clusttype_name.cluster_list\n";
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
print "(Approximate) Number of genes in clusters = $clustsum\n"; 
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

exit 0;

sub procMCL {
  my $query = shift;
  my $hit = shift;
  my $eval = shift;
  $genelist{$query} = undef;
  $genelist{$hit} = undef;
  print MCL "$query $hit $eval\n";
return 1;
}

sub finishMCL {
 foreach my $k1 ( sort {$a cmp $b} keys %genelist ) {
   procMCL($k1,$k1,0) if ($genelist{$k1});
 }
 undef %genelist;

 close(MCL); #waits for child process to exit
 exit 0;

}

sub print_redos {
	open(OF, ">>./$id/$id.$e_cutoff\_$mlf\_$mlf_type\_$recip\_$clusttype_name\_$inflation.blast-redo.run") or die "couldn't open ./$id/$id.$e_cutoff\_$mlf\_$mlf_type\_$recip\_$clusttype_name\_$inflation.blast-redo.run";
	foreach my $k (sort keys %redo) {
      foreach my $kk (sort keys %{$redo{$k}}) {
		  print OF "rm $prokdir/$k/$k.$kk.hits; /home/rec3141/repo/blast_extract.pl $k $kk\n";
		  if ($authorize_redo) {
		  	unlink("$prokdir/$k/$k.$kk.hits");
			my @output = `/home/rec3141/repo/blast_extract.pl $k $kk`;
			print join("\n",@output) if ($debug>0);
			print '!';
		  }
		}
	}
	close(OF);
return 1;
}


sub procHITS {
		my $id1 = shift;
		my $id2 = shift;
		my @list = @{(shift)};
		 HITS: foreach my $hit (@list) {
		  #id2-protein $id1-protein evalue shorter-frag-frac longer-frag-frac
		  chomp $hit;
		  my ($query_id,$hit_id,$e_val,$mlfs,$mlfl)=split(/\s+/,$hit); #split the hit
		  
		  # verify input data
		  unless (defined($query_id) && 
				  defined($hit_id) && 
				  defined($mlfl) &&
				  ($mlfl >= 0) &&
				  defined($mlfs) &&
				  ($mlfs >= 0) && 
				  defined($e_val) &&
				  ($e_val >= 0)
				  ) {
	 					warn("\n$0: problem in hit line: $hit\n") if ($debug > 0);
						$redo{$id1}{$id2} = '+';
						$missing++;
						(last HITS) if ($justcheckfororphaned);
					}

		  # find if proteins not in accession_list						
			unless (  exists($genelist{$hit_id}) &&
					  exists($genelist{$query_id}) ) {
	 					warn("\n$0: orphaned: $id1 $id2 $hit\n") if ($debug > 0);
						$redo{$id1}{$id2} = '?';
						$missing++;
						(last HITS) if ($justcheckfororphaned);
					}

			  next if ($justcheckfororphaned); #short circuit
			  next if ($query_id eq $hit_id); #skip self hits (will add back in later if they have no other hits)
			  next if ( $e_val > $e_cutoff); #skip bad matches
			  
			  my ($ha,$hb) = sort { $a cmp $b} ($query_id,$hit_id);
		
			  # build a hash of hashes with good matches
			  # or pass straight to MCL
			if (( $mlf_type == 2 ) && ( $mlfl > $mlf )) {
				$mcl>0 ? procMCL($query_id,$hit_id,$e_val) : $hoh1{ $ha }{ $hb }++; #if longer is longer than cutoff
			} elsif (defined($seeded) && ( $mlfs > $mlf )) {
				if ($mcl>0) {
					  procMCL($query_id,$hit_id,$e_val);
				} else {
				  if (exists($hoh2{ $ha }{ $hb })) {
					$hoh2{ $ha }{ $hb } = $e_val if ($e_val < $hoh2{ $ha }{ $hb });
				  } else {
					$hoh2{ $ha }{ $hb } = $e_val; #if longer is not long enough but shorter is
				  }
				}
			} elsif (( $mlf_type == 1 ) && ( $mlfs > $mlf )) {
				$mcl>0 ? procMCL($query_id,$hit_id,$e_val) : $hoh1{ $ha }{ $hb }++; #if shorter is longer than cutoff
			} elsif ( $mlf_type == 0 ) {
				$mcl>0 ? procMCL($query_id,$hit_id,$e_val) : $hoh1{ $ha }{ $hb }++;  #no length restriction
			} elsif ($debug>1) {print "$mlfl < $mlfs < $mlf ?\n"}
		} #end foreach hit; end HITS
return 1;
}

sub DoCluster {
  #get each set of blast hits $k1 and $k2
  while (my ($k1, undef) = each %hoh1) {
    foreach my $k2 ( keys %{$hoh1{$k1}} ) {
      # wait until reciprocal hits found in order to cluster
      next if ($recip && $hoh1{$k1}{$k2} < 2);
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
				$clus_gene{ $gene_clus{$k1} }{$k3} = undef;
			}
			#then delete second cluster
			delete($clus_gene{ $todel });
		  }
		  # if $k2 isn't already in a cluster, add it to that of $k1
		} else {
		  $gene_clus{ $k2 } = $gene_clus{ $k1 };
		  $clus_gene{$gene_clus{$k1} }{ $k2 } = undef;
		}
		# if $k1 isn't in a cluster but $k2 is, add $k1 to it
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
      ($genelist{$k2} = undef) if defined($gene_clus{$k2});
      delete($hoh1{$k1}{$k2});
    }
    ($genelist{$k1} = undef) if defined($gene_clus{$k1});
    delete($hoh1{$k1}) unless (keys %{$hoh1{$k1}}); #?? double check
  }
return 1;

}


sub procSeeds {
print "second pass...";
# potential fragments are in %hoh2{'NC_000000'}{'NC_111111'} = 1e-30
# sorted by refseq id
# want to arrange as %hoh2{fragment}{seed}
# so later we can sort by values and choose the best 'seed' protein

 foreach my $k1 ( keys %hoh2 ) {
  foreach my $k2 ( keys %{$hoh2{$k1}} ) {
   if ( $gene_clus{$k1} ) {
    if ( $gene_clus{$k2} ) {
      # if both are already in clusters, delete the match
      delete($hoh2{$k1}{$k2});
    } else {
      # if only first is in a cluster, move it to second position
      # so it becomes $hoh2{fragment}{seed}
      $hoh2{$k2}{$k1} = $hoh2{$k1}{$k2};
      delete $hoh2{$k1}{$k2};
    }
   } elsif ( $gene_clus{$k2} ) {
      undef; # if only second is in a cluster, keep it 
   } else {
    delete $hoh2{$k1}{$k2}; #if neither are in clusters, delete the match
   }
  }
  delete($hoh2{$k1}) unless (keys %{$hoh2{$k1}});
 }

	### after filtering, find best hit for each remaining protein
	### and add protein to cluster of best hit
	my $seeds=0;

 foreach my $k1 ( keys %hoh2 ) {
  foreach my $k2 ( sort { ($hoh2{$k1}{$a} || 1) <=> ($hoh2{$k1}{$b} || 1) } keys %{$hoh2{$k1}} ) {
	
    $seeds++;
    $gene_clus{$k1} = $gene_clus{$k2};
    $clus_gene{ $gene_clus{$k2} }{$k1} = undef;
    print "k1\t$k1\tk2\t$k2\tgc1\t$gene_clus{$k1}\tgc2\t$gene_clus{$k2}\n" unless (defined($k1) && defined($k2) && defined($gene_clus{$k1}) && defined($gene_clus{$k2}));
    last;
	#($genelist{$k2} = undef); #should already be undef'd
  }
  ($genelist{$k1} = undef);
 }
print "$seeds genes added to seeds\n";
undef %hoh2;
return 1;
}

sub procRBH {
  my @list = @{(shift)};
  my @newlist;
  my %rbh;
  my %listhash;
  foreach my $hit (@list) {
    #id2-protein $id1-protein evalue shorter-frag-frac longer-frag-frac
    chomp $hit;
    my ($query_id,$hit_id,$e_val,$mlfs,$mlfl)=split(/\s+/,$hit); #split the hit

    next if ($query_id eq $hit_id); #skip self hits (will add back in later if they have no other hits)
    
    my ($ha,$hb) = sort { $a cmp $b} ($query_id,$hit_id);
    # build a hash of hashes with good matches
    # or pass straight to MCL
    if (( $mlf_type == 2 ) && ( $mlfl > $mlf )) {
	    $rbh{ $ha }{ $hb } = $e_val; #if longer is longer than cutoff
	    $listhash{$ha}{$hb} = $hit;
    } elsif (( $mlf_type == 1 ) && ( $mlfs > $mlf )) {
	    $rbh{ $ha }{ $hb } = $e_val; #if shorter is longer than cutoff
	    $listhash{$ha}{$hb} = $hit;
    } elsif ( $mlf_type == 0 ) {
	    $rbh{ $ha }{ $hb } = $e_val;  #no length restriction
	    $listhash{$ha}{$hb} = $hit;
    } elsif ($debug>1) {print "$mlfl < $mlfs < $mlf ?\n"}
  } #end foreach hit; end HITS

 foreach my $k1 ( keys %rbh ) {
  foreach my $k2 ( sort { ($rbh{$k1}{$a} || 1) <=> ($rbh{$k1}{$b} || 1) } keys %{$rbh{$k1}} ) {
    push(@newlist,$listhash{$k1}{$k2});
#     print $listhash{$k1}{$k2},"\n";
    last;
  }
 }
return \@newlist;
}