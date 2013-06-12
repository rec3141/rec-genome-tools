#!/usr/bin/perl
use strict;
use warnings;
use diagnostics;

use Getopt::Long;
use Pod::Usage;
use File::Slurp;
use File::Copy;
use Data::Dumper;
use GenomeTools qw(:All);
use Storable qw(dclone);

sub print_redos;
sub doCluster;
sub procMCL;
sub finishMCL;
sub procHITS;
sub procSEED;
sub procRBH;
sub outputClusters;
sub makeHistogram;
sub checkOptions;
sub setupCluster;
sub reconcile;

$| = 1; #flush buffer

### get (required) options
my ($taxalist, $query_taxid, $ecutoff, $mlf, $recip_flag);

### get (optional) options
my $help = 0;
my $debug = 0;
my $clusttype_flag = 'link';
my $mcl = 0;
my $inflation_flag = 2;
my $mlf_flag = 'long';
my $output_dir = ".";
my $output_hist = 1; #output histogram
my $output_list = 1; #output cluster list
my $output_graph = 0; #output cluster graph

my @accdirs = qw|
		/work/rec3141/gene_clustering/all_genomes/faa
		/work/rec3141/gene_clustering/all_genomes/scaffold/faa
		/work/rec3141/gene_clustering/all_genomes/contig/faa
		|; #directories where accession lists live
		# /work/rec3141/gene_clustering/all_genomes/metagenomes/faa
my $hitsdir = '/work/rec3141/gene_clustering/prok'; # dir where hits files live
my $justcheckformissing = 0; #define to avoid clustering and just find missing BLAST .hits files
my $justcheckfororphaned = 0; #define to check for BLAST hits not in accession_index, indicating an obsolete protein
my $authorize_redo = 0;# allow in-line redo of missing or outdated BLASTs

GetOptions(
		#required
	   't|taxalist=s'  	=> \$taxalist,
	   'q|query=s'		=> \$query_taxid,
	   'e|ecutoff=f'	=> \$ecutoff,
	   'l|fraglength=f'	=> \$mlf,
	   'r|reciprocal=s'	=> \$recip_flag,
	   #optional
	   'help|h' 		=> \$help,
	   'v|verbose+'  	=> \$debug,
	   'w|whichfrag:s'	=> \$mlf_flag,
	   'a|algorithm:s'	=> \$clusttype_flag,
	   'i|inflation:f'	=> \$inflation_flag,
	   'output_hist!'	=> \$output_hist,
	   'output_list!'	=> \$output_list,
	   'output_graph'	=> \$output_graph,
	   'output_dir:s'	=> \$output_dir,
	   'accdirs:s'		=> \@accdirs,
	   'hitsdir:s'		=> \$hitsdir,
	   'justcheckformissing'	=> \$justcheckformissing,
	   'justcheckfororphaned'	=> \$justcheckfororphaned,
	   'authorize-redo'			=> \$authorize_redo
	   );

 pod2usage(-verbose => 2) if $help;


######### CHECK OPTIONS ##############

my (%taxahash,%allgenelist);
my (@subject_taxids, @query_taxids, @query_groupids) = ();
my ($mlf_type, $seeded, $recip, $inflation,	$cluster_outfile);

checkOptions();

if ($debug) {
	if ($mcl) {
		print "Parameters for this run: $0 -e $ecutoff -l $mlf -w $mlf_flag -r $recip_flag -a $clusttype_flag -i $inflation\n";
	} else {
		print "Parameters for this run: $0 -e $ecutoff -l $mlf -w $mlf_flag -r $recip_flag -a $clusttype_flag\n";
	}
}

########## OK TO BEGIN ############

my $starttime=time;
# we need a hash of which cluster each gene is associated with:
my %gene_cluster; #   %gene_cluster{'NC_111000'} = 'NC_111000'
# a hash of which genes are in each cluster
my %cluster_gene; #   %cluster_gene{'NC_111123'}{'NC_111000'} = undef, %cluster_gene{'NC_111123'}{'NC_111111'} = undef

my %hoh1; #cluster hash for first pass %hoh1{'NC_000000'}{'NC_111111'} = number of times seen (1 or 2)
my %hoh2; #cluster hash for second pass %hoh2{'NC_000000'}{'NC_111111'} = 1e-30
my %redo; #missing hits files and genomes with outdated accession_list
my $missing; #missing hits files?
my @cluster_lists; # list of which genes are in each cluster (list of lists):
my $subject_taxid;
my %genelist;

foreach my $subject (@query_groupids) {

	$subject_taxid = $subject;
	print "$subject_taxid\n" if $debug>1;

	# open pipe to MCL if that's what we're using
	$cluster_outfile = "$output_dir/$query_taxid/$query_taxid.$subject_taxid.$ecutoff\_$mlf\_$mlf_type\_$recip_flag\_$clusttype_flag";
	if ($mcl) {
		$cluster_outfile .= "\_$inflation";
		next if ( -e "$cluster_outfile.cluster_list");
		open MCL, "| mcl - --abc --abc-neg-log -I $inflation -o $cluster_outfile.cluster_list" or die "can't open output: $!";
	} else {
		next if ( -e "$cluster_outfile.cluster_list");
		open(OF, ">$cluster_outfile.cluster_list") or die "couldn't open $cluster_outfile.cluster_list"; close(OF);
	}

	print "reading files and clustering $subject_taxid and $query_taxid, output to $cluster_outfile.cluster_list\n" if $debug;
	$starttime = time;
	%genelist = %{ dclone(\%allgenelist) };

	%gene_cluster = %cluster_gene = %hoh1 = %hoh2 = %redo = ();
	$missing = 0;
	
	if ( exists($taxahash{$subject_taxid}) ) {
		my ($subject_taxids_ref,$subject_allids_ref) = parse_taxaid( {'taxahash' => \%taxahash, 'id' => $subject_taxid } );
		@subject_taxids = map { @$_ } @{ $subject_taxids_ref };
		die("could not find subject taxid: ", $subject_taxid, "\n") unless scalar( @subject_taxids );
	} else {
		@subject_taxids = @query_groupids;
	}

	setupCluster( (\@subject_taxids, \@query_taxids ) ); #cluster these ids

	unless( $justcheckformissing || $justcheckfororphaned ) {
		print("done hashing at ". (time-$starttime) . "s, ") if $debug;
		doCluster();
		print("done clustering at ". (time-$starttime) . "s") if $debug;
	} #end unless
	
	print "\n" if $debug;

	print "\n\nRedo the blasts in blast-redo.run\n" if ($missing && !($authorize_redo) && $debug);
	exit 0 if ($justcheckformissing || $justcheckfororphaned);
	
	# mcl will take over from here if MCL
	if ($mcl) {
		finishMCL();
		next;
	}
	
	# Finish Clustering
	$starttime = time;
	print "found " . scalar(keys %cluster_gene) . " clusters\n" if $debug;
	
	#what's left in %hoh1 at this point? non reciprocal hits if $recip==1
	%hoh1 = (); #to free up memory
	
	### now for second pass
	procSEED() if $seeded;
	
	# OUTPUT CLUSTER LIST
	@cluster_lists = ();
	outputClusters() if $output_list;
	
	# MAKE HISTOGRAM
	makeHistogram() if ($debug || $output_hist);
}


# wrap it up: check if all needed cluster files are there
# if not, exit
# if so, do reconciliation (clustering of clusters)

foreach my $subject (@query_groupids) {
	
	my $outfile = "$output_dir/$query_taxid/$query_taxid.$subject.$ecutoff\_$mlf\_$mlf_type\_$recip_flag\_$clusttype_flag";

	if ($mcl) {
		$outfile .= "\_$inflation";
	}
	print "$outfile.cluster_list incomplete\n" and exit 0 unless ( -s "$outfile.cluster_list");

}

reconcile();
print "Run Complete.\n" if $debug;
exit 0;


################ SUBROUTINES ####################

##############################
####### CHECK OPTIONS ########
##############################

sub checkOptions {
	### TAXALIST: a taxalist file containing the refseq identifier
	%taxahash = %{ parse_taxalist( {'taxalist' => $taxalist} ) };

	unless (-d $query_taxid) { mkdir($query_taxid) or die "couldn't make directory $query_taxid"}

	my %query_all_ids;
	if ( exists( $taxahash{$query_taxid} ) ) {
		my ($query_taxids_ref,$query_allids_ref) = parse_taxaid( {'taxahash' => \%taxahash, 'id' => $query_taxid } );
		@query_taxids = map { @$_ } @{ $query_taxids_ref };
		@query_groupids = @{$taxahash{$query_taxid}};
		%query_all_ids = %{ $query_allids_ref };
		die("could not find query id: $query_taxid\n") unless scalar(@query_taxids);
	} else {
		@query_groupids = ($query_taxid);
		@query_taxids = ($query_taxid);
		$query_all_ids{$query_taxid} = undef;
	}
	
	print "reading accession index...\n" if $debug;
	my ($accgenomes_ref, $genes_ref) = parse_accession_index( {'all_ids' => \%query_all_ids, 'parse_locus' => 1} );
	my %accgenomes = %{$accgenomes_ref};

	foreach (keys %{ $genes_ref }) { $allgenelist{$_} = $_ };
	
	foreach (@query_taxids) {
	  warn("WARNING: genome $_ missing from accession_index\n") if ( !(exists($accgenomes{$_})) && $debug);
	}
	print "found " . scalar(keys %allgenelist) . " genes in accession_index files\n" if $debug;

	
	### ACCDIRS: locations of ACCESSION LISTs
	@accdirs =  split(/,/,join(',',@accdirs));
	foreach my $dir (@accdirs) {
		die("couldn't find directory $dir\n") unless ( -d $dir);
	}
	

	### E-CUTOFF: an expectation cutoff
	die("invalid -ecutoff: $ecutoff\n") unless ($ecutoff =~ /^\d.+\.?e?\d.*$/);
	
	### LENGTH FRACTION: a matching fraction cutoff
	# read in the minimum length fraction
	# 0 no restriction
	# 1 by shorter sequence
	# 2 by longer sequence
	# 3 allow seeding of good fragments
	
	die("invalid -fraglength: ",$mlf, "\n") unless ($mlf =~ /^\d*\.?\d*$/);
	while ($mlf >= 1) {$mlf=$mlf/10};
	$seeded = 0;
	if (lc($mlf_flag) eq 'seed') {
	  $seeded = 1;
	  $mlf_type = 2;
	}
	$mlf_type = 2 if (lc($mlf_flag) eq 'long');
	$mlf_type = 1 if (lc($mlf_flag) eq 'short');
	$mlf_type = 0 if ($mlf == 0.0 or lc($mlf_flag) eq 'none');
	die("invalid -whichmlf: ",$mlf_flag, "\n") unless defined($mlf_type);
	
	### RECIPROCAL: a flag to require reciprocal hits or reciprocal BEST hits
	$recip=0 if (lc($recip_flag) eq 'off');
	$recip=1 if (lc($recip_flag) eq 'on');
	$recip=0 if (lc($recip_flag) eq 'rbh');
	die("invalid -reciprocal: ",$recip_flag, "\n") unless defined($recip);
	
	### CLUSTERING: a flag to determine the clustering algorithm
	$mcl=0 if (lc($clusttype_flag) eq 'link');
	$mcl=1 if (lc($clusttype_flag) eq 'mcl');
	die("invalid clustering -algorithm: ",$clusttype_flag,"\n") unless defined($mcl);
	
	if ($mcl) {
	  $inflation = undef;
	  $inflation = 0 + $inflation_flag if ($inflation_flag =~ /[0-9.]+/);
	  die("invalid -inflation parameter: ",$inflation_flag,"\n") unless defined($inflation);
	} else { $inflation = 0 }

	return 1;
}

##############################
####### SETUPCLUSTER  ########
##############################

sub setupCluster {
	my @subjects = @{ shift(@_) };
	my @queries = @{ shift(@_) };

	foreach my $id1 (@subjects) {
	print "$subject_taxid [$id1] " if $debug;
		FILE2: foreach my $id2 (@queries) {
			if (-e "$hitsdir/$id1/$id1.$id2.hits") {
			  next if ($justcheckformissing); #short circuit
			  my @list = read_file("$hitsdir/$id1/$id1.$id2.hits");
			  my $newlist_ref;
			  if ($recip_flag eq 'rbh' ) {
				$newlist_ref = procRBH( \@list );
			  } else {
				$newlist_ref = \@list;
			  }
			  procHITS( ($id1,$id2,$newlist_ref) ); #process hits file
			} else {
			  print "-" if $debug;
			  $redo{$id1}{$id2}++;
			  $missing++;
			  next;
			}
			
			if (exists($redo{$id1}{$id2})) {
				print $redo{$id1}{$id2} if $debug;
				print_redos();
				%redo=();
			} else {
				print "." if $debug;
			}
		} #end foreach file 2
		print "\n" if $debug;
	} #end foreach file1
	return 1;
}

##############################
####### PROCRBH       ########
##############################

sub procRBH {
  my @list = @{ (shift @_) };
  my @newlist;
  my %rbh;
  my %listhash;
  foreach my $hit (@list) {
    #id2-protein $id1-protein evalue shorter-frag-frac longer-frag-frac
    chomp $hit;
    my ($query_id,$hit_id,$e_val,$short,$long)=split(/\s+/,$hit); #split the hit

	next unless ($e_val < $ecutoff);
    next if ($query_id eq $hit_id); #skip self hits
	$rbh{ $hit_id }{ $query_id } = $e_val; #if longer is longer than cutoff
	$listhash{$hit_id}{$query_id} = $hit;

   } #end foreach hit; end HITS

 #sort them and keep only the BEST hits
 foreach my $k1 ( keys %rbh ) {
  foreach my $k2 ( sort { ($rbh{$k1}{$a} || 1) <=> ($rbh{$k1}{$b} || 1) } keys %{$rbh{$k1}} ) {
    push(@newlist,$listhash{$k1}{$k2});
	print $listhash{$k1}{$k2},"\n" if ($debug > 2);
    last;
  }
 }
return \@newlist;
}

##############################
####### PRINT_REDOS   ########
##############################

sub print_redos {
	open(OF, ">>$cluster_outfile.blast-redo.run") or die "couldn't open output file $cluster_outfile.blast-redo.run\n";
	foreach my $k (sort keys %redo) {
      foreach my $kk (sort keys %{$redo{$k}}) {
		  print OF "rm $hitsdir/$k/$k.$kk.hits; blast_extract.pl $k $kk\n";
		  if ($authorize_redo) {
		  	unlink("$hitsdir/$k/$k.$kk.hits");
			my @output = `blast_extract.pl $k $kk`;
			print join("\n",@output) if ($debug>1);
			print '!' if $debug;
		  }
		}
	}
	close(OF);
return 1;
}

##############################
####### PROCHITS      ########
##############################

sub procHITS {
		my $id1 = shift;
		my $id2 = shift;
		my @list = @{ shift(@_) };
		 HITS: foreach my $hit (@list) {
		  #id2-protein $id1-protein evalue shorter-frag-frac longer-frag-frac
		  chomp $hit;
		  my ($query_id,$hit_id,$e_val,$short,$long)=split(/\s+/,$hit); #split the hit
		  
		  # verify input data
		  unless (defined($query_id) && 
				  defined($hit_id) && 
				  defined($long) &&
				  ($long >= 0) &&
				  defined($short) &&
				  ($short >= 0) && 
				  defined($e_val) &&
				  ($e_val >= 0)
				  ) {
	 					warn("\n$0: problem in hit line: $hit\n") if ($debug > 1);
						$redo{$id1}{$id2} = '+';
						$missing++;
						(last HITS) if ($justcheckfororphaned);
					}

		  # find if proteins not in accession_list						
			unless (  exists($genelist{$hit_id}) &&
					  exists($genelist{$query_id}) ) {
	 					warn("\n$0: orphaned: $id1 $id2 $hit\n") if ($debug > 1);
						$redo{$id1}{$id2} = '?';
						$missing++;
						(last HITS) if ($justcheckfororphaned);
					}

			  next if ($justcheckfororphaned); #short circuit
			  next if ($query_id eq $hit_id); #skip self hits (will add back in later if they have no other hits)
			  next if ( $e_val > $ecutoff); #skip bad matches
			  
			  my ($ha,$hb);
			  if ($recip_flag eq 'rbh') {
				($ha,$hb) = ($hit_id, $query_id);
			  } else {
			  	($ha,$hb) = sort ($query_id, $hit_id);
			  }
			  # build a hash of hashes with good matches
			  # or pass straight to MCL
			if (( $mlf_type == 2 ) && ( $long > $mlf )) {
				$mcl ? procMCL($query_id,$hit_id,$e_val) : $hoh1{ $ha }{ $hb }++; #if longer is longer than cutoff
			} elsif ($seeded && ( $short > $mlf )) {
				if ($mcl) {
					  procMCL($query_id,$hit_id,$e_val);
				} else {
				  if (exists($hoh2{ $ha }{ $hb })) {
					$hoh2{ $ha }{ $hb } = $e_val if ($e_val < $hoh2{ $ha }{ $hb });
				  } else {
					$hoh2{ $ha }{ $hb } = $e_val; #if longer is not long enough but shorter is
				  }
				}
			} elsif (( $mlf_type == 1 ) && ( $short > $mlf )) {
				$mcl ? procMCL($query_id,$hit_id,$e_val) : $hoh1{ $ha }{ $hb }++; #if shorter is longer than cutoff
			} elsif ( $mlf_type == 0 ) {
				$mcl ? procMCL($query_id,$hit_id,$e_val) : $hoh1{ $ha }{ $hb }++;  #no length restriction
			} elsif ($debug > 2) {print "$long < $short < $mlf ?\n"}
		} #end foreach hit; end HITS
return 1;
}

##############################
####### PROCMCL       ########
##############################

sub procMCL {
  my $query = shift;
  my $hit = shift;
  my $eval = shift;
  $genelist{$query} = undef;
  $genelist{$hit} = undef;
  print MCL "$query $hit $eval\n";
return 1;
}

##############################
####### FINISHMCL     ########
##############################

sub finishMCL {
 foreach my $k1 ( sort {$a cmp $b} keys %genelist ) {
   procMCL($k1,$k1,0) if defined($genelist{$k1});
 }
 undef %genelist;

 close(MCL); #waits for child process to exit
 return 1;

}


##############################
####### RECONCILE     ########
##############################

sub reconcile {
	my %reco = ();
	%hoh1 = %gene_cluster = %cluster_gene = %genelist = ();
	%genelist = %{ dclone(\%allgenelist) };
	$cluster_outfile = "$output_dir/$query_taxid/$query_taxid.$ecutoff\_$mlf\_$mlf_type\_$recip_flag\_$clusttype_flag";
	$cluster_outfile .= "\_$inflation" if $mcl;

	print "$cluster_outfile.cluster_list.tmp exists\n" and exit 0 if ( -e "$cluster_outfile.cluster_list.tmp" );

	open(OF,">","$cluster_outfile.cluster_list.tmp");close(OF);

	print "reconciling cluster_lists...\n" if $debug;
	
	my $print_recip_flag = $recip_flag;
	if ($recip_flag eq 'rbh') {
		$recip = 1;
		$recip_flag = 'RBH';
	}

	$starttime = time;
	print "genome\t\%reco\t\%hoh1pre\t\%hoh1post\t\%cluster_gene\t\%gene_cluster\ttime\n" if $debug;
	foreach my $subject (@query_groupids) {

		print "$subject\t" if $debug;
		my $infile = "$output_dir/$query_taxid/$query_taxid.$subject.$ecutoff\_$mlf\_$mlf_type\_$print_recip_flag\_$clusttype_flag";
		$infile .= "\_$inflation" if $mcl;

		my @clusters = read_file("$infile.cluster_list",chomp=>1);

		#ok, instead of saving every $hit and $key,
		#save references to anonymous hash of %genelist key
		#access using keys or each
		foreach (@clusters) {
			next unless $_;
			my @prots = split;
# 			my $key = shift(@prots); #old safe
			my $keysafe = shift(@prots);
			my $key = \$allgenelist{$keysafe};
			my @nohit = ();
			foreach my $hitsafe (@prots) {
				my $hit = \$allgenelist{$hitsafe};
				if ($recip) {
				  if (exists($reco{$hit}{$key})) {
# 					$hoh1{$hit}{$key} = 2;
					$hoh1{$hitsafe}{$keysafe} = 2;
					delete($reco{$hit}{$key});
				  } else {
					$reco{$key}{$hit} = undef;
				  }
				  
				} else {$hoh1{$key}{$hit}++;}
			} #end foreach @prots

			foreach my $hit (keys %{$reco{$key}}) {
			  delete($reco{$hit}{$key}) if exists($reco{$hit}{$key});
			  delete($reco{$hit}) unless keys %{$reco{$hit}};
			}
		} #end foreach @clusters

		print scalar(keys %reco),"\t",scalar(keys %hoh1),"\t" if $debug;
		doCluster();
		print scalar(keys %hoh1),"\t",scalar(keys %cluster_gene),"\t",scalar(keys %gene_cluster),"\t",(time - $starttime),"\n" if $debug;
	}
		print Dumper %cluster_gene;
		print Dumper %gene_cluster;

	$recip_flag = $print_recip_flag;
	@cluster_lists = ();

	outputClusters() if $output_list;

	unlink("$cluster_outfile.cluster_list.tmp");
	return 1;
}

##############################
####### DOCLUSTER     ########
##############################

sub doCluster {
  # $gene_cluster{'NP_111123'} = 'NP_000123';
  # $cluster_gene{'NP_000123'}{'NP_111123'} = undef
  
  foreach my $k1 ( keys %hoh1 ) {
    foreach my $k2 ( keys %{$hoh1{$k1}} ) {
	  
      # if $recip then wait until reciprocal hits found in order to cluster
      next if ($recip && $hoh1{$k1}{$k2} < 2);
      
      my $gc1 = \$gene_cluster{$k1};
      my $gc2 = \$gene_cluster{$k2};
      
      # if $k1 is already in a cluster
	  if ( $gene_cluster{$k1} ) {
		# and if $k2 is already in a cluster
		if ( $gene_cluster{$k2} && ($recip_flag ne 'rbh') ) {
		# then they need to be merged
		# unless we want unmerged first-pass RBH clusters
		# don't do anything if they're already in the same cluster
		  if ( $gene_cluster{$k1} ne $gene_cluster{$k2} ) {
			#and add all of the genes in the second cluster to the first
			foreach my $k3 (keys %{ $cluster_gene{ $gene_cluster{$k2} }}) {
				$gene_cluster{$k3} = $gene_cluster{$k1};
				$cluster_gene{ $gene_cluster{$k1} }{$k3} = undef;
			}
			$gene_cluster{$k2} = $gene_cluster{$k1};
			$cluster_gene{ $gene_cluster{$k1} }{$k2} = undef;
			#then delete second cluster
			delete($cluster_gene{ $k2 });
		  } else {undef;}
		  # if $k2 isn't already in a cluster, add it to that of $k1
		} else {
		  $gene_cluster{ $k2 } = $gene_cluster{ $k1 };
		  $cluster_gene{ $gene_cluster{$k1} }{ $k2 } = undef;
		}
		# if $k1 isn't in a cluster but $k2 is, add $k1 to it
		# unless first-pass RBH because we want to keep them separate
	  } elsif ( $gene_cluster{$k2} && ($recip_flag ne 'rbh') ) {
			$gene_cluster{ $k1 } = $gene_cluster{ $k2 };
			$cluster_gene{ $gene_cluster{$k2} }{ $k1 } = undef;
			
	  #if neither gene1 nor gene2 are in a cluster, add them both to a new cluster
      } else {
		$gene_cluster{$k1} = $gene_cluster{$k2} = $k1;
		$cluster_gene{$k1}{$k2} = undef;
      }
      $genelist{$k2} = undef if $gene_cluster{$k2};
      delete($hoh1{$k1}{$k2});
    }
    ($genelist{$k1} = undef) if $gene_cluster{$k1};
    delete($hoh1{$k1}) unless (keys %{$hoh1{$k1}}); #?? double check
  }
return 1;

}

##############################
####### PROCSEEDS     ########
##############################

sub procSEED {
print "second pass..." if $debug;
# potential fragments are in %hoh2{'NC_000000'}{'NC_111111'} = 1e-30
# sorted by refseq id
# want to arrange as %hoh2{fragment}{seed}
# so later we can sort by values and choose the best 'seed' protein

 foreach my $k1 ( keys %hoh2 ) {
  foreach my $k2 ( keys %{$hoh2{$k1}} ) {
   if ( $gene_cluster{$k1} ) {
    if ( $gene_cluster{$k2} ) {
      # if both are already in clusters, delete the match
      delete($hoh2{$k1}{$k2});
    } else {
      # if only first is in a cluster, move it to second position
      # so it becomes $hoh2{fragment}{seed}
      $hoh2{$k2}{$k1} = $hoh2{$k1}{$k2};
      delete $hoh2{$k1}{$k2};
    }
   } elsif ( $gene_cluster{$k2} ) {
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
    $gene_cluster{$k1} = $gene_cluster{$k2};
    $cluster_gene{ $gene_cluster{$k2} }{$k1} = undef;

    if ($debug > 1) {
    	print "k1\t$k1\tk2\t$k2\tgc1\t$gene_cluster{$k1}\tgc2\t$gene_cluster{$k2}\n" unless (defined($k1) && defined($k2) && defined($gene_cluster{$k1}) && defined($gene_cluster{$k2}));
	}
    last;
  }
  ($genelist{$k1} = undef);
 }
print "$seeds genes added to seeds\n" if $debug;
%hoh2 = ();

return 1;
}

##############################
####### OUTPUTCLUSTER   ######
##############################

sub outputClusters {
	### get ready for printing
	$starttime = time;
	print "making cluster list..." if $debug;
	
	#make @cluster_lists
	map {
	  push(
		@cluster_lists, 
		[$_, sort {$a cmp $b} keys %{ $cluster_gene{$_} } ] 
	  )
	} keys %cluster_gene;
	%cluster_gene = ();


	#first sort the cluster list by order of genes
# 	@cluster_lists = sort {"@$a" cmp "@$b"} @cluster_lists;

	#then sort the cluster list by number of genes in the cluster
	@cluster_lists = sort {@{$b} <=> @{$a}} @cluster_lists;

	# add the unclustered genes/singletons to @cluster_lists
	unless (lc($recip_flag) eq 'rbh') {
	  foreach my $k1 ( sort {$a cmp $b} keys %genelist ) {
	  push(@cluster_lists, ( [$k1] ) ) if defined($genelist{$k1});
	  }
	}
 	%genelist = ();
	
	print "done in " . (time - $starttime). "s\n" if $debug;

	open(OF, ">","$cluster_outfile.cluster_list.tmp") or warn "couldn't open $cluster_outfile.cluster_list.tmp";
	for my $cluster ( @cluster_lists ) {
		print OF "@$cluster\n";
	}
	print OF "\n"; #this marks file in case of empty cluster_list
	close(OF);
	move("$cluster_outfile.cluster_list.tmp","$cluster_outfile.cluster_list");
	print "output $cluster_outfile.cluster_list\n" if $debug;
	
	return 1;
}

##############################
####### MAKEHISTOGRAM   ######
##############################

sub makeHistogram {

	my @hist;
 	my $num_clus=scalar(@cluster_lists);
	my $num_genes=0;
	my $max_genes_cluster = 10000; # this sets the largest bin for the histogram.
	foreach my $k1 ( @cluster_lists ) {
		$num_genes=scalar(@{$k1});
		$hist[$num_genes]++;
	}
	
	for (my $i=0;$i<$max_genes_cluster;$i++) {
	  $hist[$i] = 0 unless defined($hist[$i]);
	}
	
	if ($debug) {
		print "Number of clusters = ". ($num_clus - $hist[1]) . "\n";
		my $clustsum; 
		for(my $i=2;$i<=$max_genes_cluster;$i++) {$hist[$i] ? $clustsum+=$hist[$i]*$i : $clustsum+=0};
		print "(Approximate) Number of genes in clusters = $clustsum\n"; 
		print "Number of singletons = $hist[1]\n";
		print "Total genes = " . ($hist[1]+$clustsum) . "\n";
	}
	
	if ($output_hist) {
		my $outfile = "$output_dir/$query_taxid/$query_taxid.cluster_hist";
		open(OF, ">>$outfile") or warn "couldn't open histogram file: $outfile";
		my $total_clus = 0;
		my $over_clus = 0;
		print OF "$query_taxid\t$subject_taxid\t$ecutoff\t$mlf\t$mlf_flag\t$recip_flag\t" , $seeded ? "seeded\t" : "noseed\t";
		for (my $ii = 1; $ii < $max_genes_cluster; $ii++) {
		  if ($ii > $max_genes_cluster) {
			$over_clus += $hist[$ii];
			next unless ($ii == $max_genes_cluster-1);
			print OF "\t$over_clus";
			$total_clus=$total_clus+$hist[$ii];
		  } elsif (defined($hist[$ii])) {
			print OF "\t$hist[$ii]";
			$total_clus=$total_clus+$hist[$ii];
		  } else {print OF "\t0";}
		}
		print OF "\n";
		close(OF);
		warn("Number of clusters incongruent $total_clus $num_clus\n") if ($total_clus != $num_clus and $debug);
		print "output $outfile\n\n" if $debug;
	}
	return 1;
}



##############################
####### END           ########
##############################


 __END__
 
=head1 NAME

distcluster - distributed genome clustering tool

=head1 SYNOPSIS

This program clusters the genes in a genome using either single-link clustering or mcl clustering. 
It can be used in parallel

=item B<example:>

	distcluster.pl -t taxalist-Bacillus -s NC_000123 -q SR_110001 -e 1e-20 -f 0.8 -r seed


=head2 Required:

-t|taxalist: a file containing the refseq identifier

-q|query: query taxa, a refseq identifier or group id from a taxalist

-e|ecutoff: a BLAST expectation value cutoff (of the form '0.001' or '1e-3')

-l|fraglength: minimum fragment length cutoff (between 0 and 1)

-r|reciprocal: require reciprocal BLAST hits (one of 'off','on','rbh','seed')

see help page (--help) for more options


=head1 OPTIONAL ARGUMENTS

-help|h                 print help

-v|verbose              zero for silent, once for normal, more for debugging


-a|algorithm            clustering algorithm [link]

-w|whichfrag            which fragment to consider in length cutoff [long] 

-i|inflation            inflation parameter for mcl clustering [2]

-accdirs                directories where accession_index files live

-hitsdir                directory to find blast hit files (.hits)

-justcheckformissing    short circuit and just look for missing blast hits

-justcheckfororphaned   short circuit and just look for orphaned blast hit files

-authorize-redo         allow in-place re-blasting of missing/orphaned .hits files

-nooutput_list          don't output cluster list

-nooutput_hist          don't output histogram

-output_graph           output cluster graph [no]

-output_dir             base directory for output files [./]

=head1 OPTION DETAILS

=item B<-whichfrag>

	a flag to determine which gene to consider for the length of match

	`none' turns off match length cutoff (overrides -fraglength and sets cutoff to 0)

	`long' requires -fraglength to match longer sequence

	`short' requires -fraglength to match shorter sequence 

	`seed' tells whether to allow good clusters to recruit shorter fragments on 2nd pass (implies `long')

=item B<-reciprocal>

	a flag to require reciprocoal hits

	`on' to require reciprocal hits (if -fraglength=seed, recruited fragments do not require reciprocal hits, others do)

	`off' to NOT require reciprocal hits     

	`rbh' to use only Reciprocal Best/Top Hits (RBHs)

=item B<-algorithm>

	a flag to determine the type of clustering

	`link' to perform SINGLE-LINKAGE clustering

	`mcl' to perform MCL clustering using the Markov Clustering algorithm (external)

=item B<-inflation>

	an inflation parameter for the MCL algorithm. MCL author suggests trying 1.4, 2, 4, and 6

=head1 OUTPUTS

=item B<*.cluster_list>
a list of clusters

=item B<*.cluster_hist>
a histogram of cluster sizes

=head1 RESOURCE UTILIZATION
520 genomes were semi-clustered in 5 hours on 36 CPUs using 1200Mb of memory each
final clustering

=head1 CHANGES

=item B<1/25/2013:> enable distributed clustering

=item B<1/20/2013:> enable RBH filtering

=item B<3/30/2012:> enable MCL clustering

=item B<2/01/2012:> enable a second pass that allows short matches to join clustered seeds

=cut
