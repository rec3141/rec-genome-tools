package GenomeTools;
use strict;
use Exporter qw(import);

our $VERSION     = 1.00;
our @ISA         = qw(Exporter);
our @EXPORT_OK   = qw(parse_taxaid  parse_taxalist  parse_accession_index  parse_proks  parse_faa  parse_ffn  parse_ptt  count_unique  parse_cluster_index  check_for_clusters  get_proks  get_hits  get_from_genes  get_from_ncbi);
our %EXPORT_TAGS = ( DEFAULT => [qw( parse_taxaid  parse_taxalist  parse_accession_index  parse_proks  get_proks )],
                 All    => [qw( parse_taxaid  parse_taxalist  parse_accession_index  parse_proks  parse_faa  parse_ffn  parse_ptt  count_unique  parse_cluster_index  check_for_clusters  get_proks  get_hits  get_from_genes  get_from_ncbi )]);

###################
### SUBROUTINES ###
###################


#------------------------
# SUB PARSE_TAXALIST
#------------------------
sub parse_taxalist {
	my %params = %{ shift(@_) };
	my $taxalist = $params{'taxalist'};
	my %_taxahash = ();
	
	die("$0 could not find -taxalist: ", $taxalist, "\n") unless ( -e $taxalist);
	
	# read in the taxalist to get list of genomes
	open(IF, "<$taxalist") || die("couldn't open file $taxalist\n");
	while (<IF>) {
	  chomp;
	  my @taxa = split(/\s/,$_);
	  my $tag = shift @taxa;
	  next unless $tag;
	  push(@{$_taxahash{$tag}},@taxa);

	}
	close(IF);
	
	return \%_taxahash;
}

#------------------------
# SUB PARSE_CLUSTER_INDEX
#------------------------
sub parse_cluster_index {
	my %params = %{ shift(@_) };
	my $opt_cluster_index = $params{'opt_cluster_index'};
	my %_wanted_clusters = ();
	my $cluster_file = $params{'cluster_file'};
	
	if ($opt_cluster_index =~ m/\d+/) {
		if ($opt_cluster_index =~ m/,/) {
		map { $_wanted_clusters{$_}++ } split(/,/,$opt_cluster_index);
	  } else {
		$_wanted_clusters{$opt_cluster_index}++;
	  }
	} else {
	  die("$0 failed to understand input argument: ", $opt_cluster_index, "\n");
	}

	if ( defined($_wanted_clusters{0}) ) {
		open(IF1, "< $cluster_file")  || die("couldn't open file $cluster_file: $!\n");
		while (<IF1>) {}; #get last line number into $.
		map {$_wanted_clusters{$_}++} (1..$.);
		delete $_wanted_clusters{0};
		close(IF1);
	}

	return \%_wanted_clusters;
}

#------------------------
# SUB CHECK_FOR_CLUSTERS
#------------------------
sub check_for_clusters {
	my %params = %{ shift(@_) };
	my %wanted_clusters = %{ $params{'wanted_clusters'} };
	my $opt_maxvals = $params{'opt_maxvals'};
	my $cluster_file = $params{'cluster_file'};
	
	my $clus_num = scalar(keys %wanted_clusters);

	my %_prot_seen = ();
	my $clus_found = 0;
	open(IF1, "< $cluster_file")  || die("couldn't open file $cluster_file: $!\n");
	while (<IF1>) {
	  next unless defined($wanted_clusters{$.});
	  $clus_found++;
	  chomp;
	  my @protline = split(/\s/,$_);
	  # to limit the number of responses (for speeds sake)
	  if ( ($opt_maxvals > 0) && (scalar(@protline)>$opt_maxvals) ) {
		my $i=0;
		while ($i < $opt_maxvals) {
		  my $tmp = $protline[rand @protline];
		  unless ($_prot_seen{$tmp}) { 
		  	$_prot_seen{$tmp}++; 
		  	$i++; 
		  }
		}
	  } else { map { $_prot_seen{$_}++ } @protline}
	
	}
	close(IF1);
	
	die("couldn't find all requested clusters\n") if ($clus_found != $clus_num);

	return \%_prot_seen;
}

#------------------------
# SUB PARSE_TAXAID
#------------------------
sub parse_taxaid {
	my %params = %{ shift(@_) };
	my %taxahash = %{ $params{'taxahash'} };
	my $id = $params{'id'};
	my %_all_ids = (); # %_all_ids{'NC_000123'}=undef
	my @_ids=(); # [ [NC_000123 NC_000124] [NC_123123 NC_123321] ]

	#die if not properly formatted id
	die("$0 invalid id: ", $id, "\n") unless ($id =~ /^\D{2}\_\d+$/);
	#die if id not found
	die("$0: id $id not found\n") unless scalar( @{$taxahash{$id}} );
	
	map { $_all_ids{$_} = undef } @{ $taxahash{$id} };
	
	foreach my $tag ( @{ $taxahash{$id} } ) {
		if ( exists($taxahash{$tag}) ) {
			map { $_all_ids{$_} = undef } @{$taxahash{$tag}};
			push(@_ids,[ @{ $taxahash{$tag} } ]);
		} else {
			push(@_ids,[ $tag ]);
		}
	}
	
	return( (\@_ids, \%_all_ids) );
}

#------------------------
# SUB PARSE_ACCESSION_INDEX
#------------------------
# returns a hash ref of all seen proteins by genome, and an array ref of seen proteins
sub parse_accession_index {
	my %params = %{ shift(@_) };
	my %all_ids = %{ $params{'all_ids'} };
	my $parse_locus;
	if ( exists( $params{'parse_locus'} ) ) {
		$parse_locus = $params{'parse_locus'};
	} else { $parse_locus = undef }
	
	my @accdirs;
	if ( exists $params{'accdirs'} ) {
		@accdirs = @{ $params{'accdirs'} }
	} else {
		@accdirs = qw|/work/rec3141/gene_clustering/all_genomes/faa /work/rec3141/gene_clustering/all_genomes/scaffold/faa /work/rec3141/gene_clustering/all_genomes/contig/faa|;
	}
	
	my %prot_seen = ();
	if ( exists( $params{'prot_seen'} ) ) {
		%prot_seen = %{ $params{'prot_seen'} }
	} else { %prot_seen = () }

	my %_accgenomes=();
	my %_genes=();

	# print "reading accession index...";
	foreach my $accdir (@accdirs) {
	  open(IF,"<$accdir/accession_index") or die "couldn't open file $accdir/accession_index\n";
	  while (<IF>) {
		chomp;
		my @acc = split(/\s+/,$_);
		my $genome = shift(@acc);
		next unless exists($all_ids{$genome});
		map { $prot_seen{$_} || ($parse_locus) ? $_accgenomes{$genome}{$_} = undef : undef } @acc;
		map { $prot_seen{$_} || ($parse_locus) ? $_genes{$_}{'taxonid'} = $genome : undef } @acc;
	  }
	  close(IF);
	}
	return( (\%_accgenomes, \%_genes) );

}

#------------------------
# SUB PARSE_PROKS
#------------------------
sub parse_proks {

	my %_ncbi = ();
	
	my $ncbifile = "/work/rec3141/gene_clustering/all_genomes/prokaryotes.txt";
	open( IF, "<$ncbifile" ) or die("couldn't open file $ncbifile\n");
	while (<IF>) {
		next if $.==1; #skip header line
		my ($taxaname, undef , $projectid, $group ,$subgroup,undef,undef,$chromosomes,undef,$plasmids,undef,$wgs) = split( "\t", $_ );
	
		unless ($chromosomes =~ m/-/) {
		foreach my $replicon ( split(/,/,$chromosomes) ) {
			my ($accession,undef) = split(/\./,$replicon);
			$_ncbi{$accession}{'taxaname'} = $taxaname;
			$_ncbi{$accession}{'refseq'} = $accession;
			$_ncbi{$accession}{'replicon'} = 'chromosome';
			$_ncbi{$accession}{'group'} = $group;
			$_ncbi{$accession}{'subgroup'} = $subgroup;
			$_ncbi{$accession}{'projectid'} = $projectid;
		}}
	
		unless ($plasmids =~ m/-/) {
		foreach my $replicon ( split(/,/,$plasmids) ) {
			my ($accession,undef) = split(/\./,$replicon);
			$_ncbi{$accession}{'refseq'} = $accession;
			$_ncbi{$accession}{'taxaname'} = $taxaname;
			$_ncbi{$accession}{'replicon'} = 'plasmid';
			$_ncbi{$accession}{'group'} = $group;
			$_ncbi{$accession}{'subgroup'} = $subgroup;
			$_ncbi{$accession}{'projectid'} = $projectid;
		}}
		
		unless ($wgs =~ m/-/) {
			$wgs =~ m/(\w{4})(\d{2})/;
			my $accession = 'NZ_' . $1 . '00000000';
			$_ncbi{$accession}{'refseq'} = $accession;
			$_ncbi{$accession}{'taxaname'} = $taxaname;
			$_ncbi{$accession}{'replicon'} = 'plasmid';
			$_ncbi{$accession}{'group'} = $group;
			$_ncbi{$accession}{'subgroup'} = $subgroup;
			$_ncbi{$accession}{'projectid'} = $projectid;
		}
			
	}
	close IF;

	return(\%_ncbi);
}

#------------------------
# SUB GET_HITS
#------------------------
sub get_hits {
# honestly i'm not sure what this is supposed to be useful for....
# maybe for subsetting a cluster into MCL??
# ps not sure it works properly
	my $taxonids_ref = shift;
	my @ids = keys %{$taxonids_ref};
	my $cluster_list_ref = shift;
	my @_proteins = @{$cluster_list_ref};
	my $genes_ref = shift;
	my %genes = %{ $genes_ref };
	my $prokdir = '/work/rec3141/gene_clustering/prok';
	foreach my $subject (@ids) {
		foreach my $query (@ids) {
			my %evalues = ();
			open(IF,"<$prokdir/$subject/$subject.$query.hits") or die("couldn't open hits file $prokdir/$subject/$subject.$query.hits\n");
			while(<IF>) {
				chomp;
				my @line = split(/\s+/,$_);
				next unless (exists( $genes{$line[0]} ) && exists( $genes{$line[1]} ));
				$evalues{$line[0]}{$line[1]} = $line[2];
			}
			close(IF);
		
			my @w = keys(%evalues);
			foreach my $e (keys %evalues) {
				push (@w,keys %{$evalues{$e}});
			}
			
			my %union = ();
			my %isect = ();
			foreach my $f (@w) { $union{$f} = 1 }
			foreach my $g (@_proteins) {
				if ( $union{$g} ) { $isect{$g} = 1 }
				$union{$g} = 1;
			}		

			foreach my $i (keys %isect) {
				foreach my $j (keys %isect) {
					if (exists($evalues{$i}{$j})) {print join("\t",$i,$j,$evalues{$i}{$j}),"\n";}
				}
			}
		} #end foreach my query
	} #end foreach my subject
	
	return;
} #end sub

#------------------------
# SUB GET_DIRS
#------------------------
sub get_dirs {

	my $type = shift;
	my @ids = @{ shift(@_) };

	my %_dirs;

	my @moldirs = qw|/work/rec3141/gene_clustering/all_genomes /work/rec3141/gene_clustering/all_genomes/scaffold /work/rec3141/gene_clustering/all_genomes/contig|;
	foreach my $id (@ids) {
		my $dir = undef;
		foreach my $moldir (@moldirs) {
			next unless ( -e "$moldir/$type/$id.$type");
			$dir = $moldir;
			last;
		}
		
		if (defined($dir)) {
			$_dirs{$id} = $dir;
		} else {
			die ("missing $type file for $id\n");
		}		
	}
	return(\%_dirs);
}

#------------------------
# SUB PARSE_FAA
#------------------------
sub parse_faa {
	my %params = %{ shift(@_) };
	my %genes = %{ $params{'genes'} };
	my @ids = keys %{ $params{'accgenomes'} };
	my @outputs = @{ $params{'outputs'} };
	
	my $dirs_ref = get_dirs( ('faa',\@ids) );
	my %dirs = %{ $dirs_ref };
	my %_pids = ();
	
	foreach my $id (@ids) {
		my $faaio = Bio::SeqIO->new(-file => "$dirs{$id}/faa/$id.faa", -format=>"fasta");
		# >gi|168697871|ref|ZP_02730148.1| hypothetical protein GobsU_00005 [Gemmata obscuriglobus UQM 2246]

		while (my $faa_obj = $faaio->next_seq) {
			my @display_ids = split(/\|/,$faa_obj->display_id);
			chomp @display_ids;
			my $refseq = $display_ids[3];
			my $pid = $display_ids[1];
			next unless exists($genes{$refseq});
			
			$genes{$refseq}{'pid'} = $pid;
			$_pids{$pid} = $refseq;
			# $pids{'330444846'} = 'ZP_08308501.1'

			if (grep {$_ =~ m/fastaaa/} @outputs ) {
				# $genes{'ZP_08308501.1'}{'seq_obj'} = {seq_obj}
				$faa_obj->display_id( $faa_obj->primary_id . "mol\|$id" );
				push(@{ $genes{$refseq}{'fastaaa'} },$faa_obj);
			}

		} #end while sequence

	} #end each taxaid
	
	return(\%_pids);
	
} #end sub parse_faa
      
#------------------------
# SUB PARSE_PTT
#------------------------
sub parse_ptt {
	my %params = %{ shift(@_) };
	my @ids = keys %{ $params{'accgenomes'} };
	my %pids = %{ $params{'pids'} };
	my @outputs = @{ $params{'outputs'} };
	my %genes = %{ $params{'genes'} };
	my %ncbi = %{ $params{'ncbi'} };
	my $dirs_ref = get_dirs( ('ptt',\@ids) );
	my %dirs = %{ $dirs_ref };
	my $parse_locus = $params{'parse_locus'};
	
	my %_locs = ();
	my %_loci = ();
	
	foreach my $id (@ids) {	
		# process ptt file
		open(PTT,"<","$dirs{$id}/ptt/$id.ptt") || die("couldn't open file $dirs{$id}/ptt/$id.ptt\n");

		PTTLINE: while (<PTT>) {
			if ($. == 1) {
				my @line = split(/\t/,$_);
				$ncbi{$id}{'length'} = $line[-1];
			}
 			next unless ($_ =~ m/^\d+\.\.\d+/);
			chomp;
			my $pttline = $_;
			
			# Erwinia pyrifoliae Ep1/96 chromosome, complete genome - 1..4026322
			# Location        Strand  Length  PID     Gene    Synonym Code    COG     Product

			my ($pttloc, $strand, $length, $pttpid, $pttgene, $locus, $code, $pttcog, $pttdesc) = split(/\t/,$pttline);

			next unless exists($pids{$pttpid});
			my $refseq = $pids{$pttpid};
			next unless defined($refseq);
			
			$genes{$refseq}{'strand'} = $strand if (grep { $_ =~ m/strand/ } @outputs);
			$genes{$refseq}{'gene'}   = $pttgene if (grep { $_ =~ m/gene/ } @outputs);
			$genes{$refseq}{'desc'}   = $pttdesc if (grep { $_ =~ m/desc/ } @outputs);
			$genes{$refseq}{'length'} = $length if (grep { $_ =~ m/length/ } @outputs);
			$genes{$refseq}{'locus'}= $locus if (grep { $_ =~ m/locus/ } @outputs);
			$_loci{lc($locus)} = $refseq if ($parse_locus);

			my($locs,$loce) = split(/\.\./,$pttloc);
			$genes{$refseq}{'locs'} = $locs;
			$genes{$refseq}{'loce'} = $loce;
			push(@{ $_locs{$id}{$locs}{$loce} },$refseq);
			
			if ( grep { $_ =~ m/cog/ } @outputs ) {
				if ($pttcog =~ m/(COG\d+)([A-Z]+)/) {
					$genes{$refseq}{'cog'} = $1;
					$genes{$refseq}{'cogfunc'} = $2;
				} else {
					$genes{$refseq}{'cog'} = $pttcog;
					$genes{$refseq}{'cogfunc'} = $pttcog;
				}
			}
		}
		close(PTT);

	} #end foreach my $id
	
	return( (\%_locs,\%_loci) );
} # end parse_ptt

#------------------------
# SUB PARSE_FFN
#------------------------
sub parse_ffn {
	my %params = %{ shift(@_) };
	my @ids = keys %{ $params{'accgenomes'} };
	my %genes = %{ $params{'genes'} };
	my %locs = %{ $params{'locs'} };
	my $dirs_ref = get_dirs( ('ptt',\@ids) );
	my %dirs = %{ $dirs_ref };


	foreach my $id (@ids) {
		#get nucleotide sequences if so requested
		my %ffnhash=(); #temporary for each id
		my $ffn;
		my $ffnio = Bio::SeqIO->new(-file => "$dirs{$id}/ffn/$id.ffn",-format=>"fasta");
		while (my $ffn_obj = $ffnio->next_seq) {
			my ($genome,$location);
			if ($ffn_obj->display_id =~ /^gi/) {
			  (undef,undef,undef,$genome,$location) = split(/\|/,$ffn_obj->display_id);
			} elsif ($ffn_obj->display_id =~ /^ref/) {
			  (undef,$genome,$location) = split(/\|/,$ffn_obj->display_id);
			} else {warn "weird display id: ",$ffn_obj->display_id,"\n";}

			chomp $location;
			if ($location !~ /^:(.+)$/) {
				warn "no location given: ",$ffn_obj->display_id,"\n";
				next;
			}

			my @pieces = ( $1 =~ /(c?\d+\-\d+)/g );
			
			# some genes are in pieces
			my @ffn = ();
			my $complement;
			foreach my $piece (@pieces) {
				my @abc = ($piece =~ m/(c?)(\d+)\-(\d+)/);
				$complement = shift(@abc); # all pieces have 'c'
				push(@ffn,@abc);
			}
			
			@ffn = reverse(@ffn) if ($complement eq 'c');

# 			print"$id\t@ffn\n";
			next unless exists( $locs{ $id }{ $ffn[0] }{ $ffn[-1] } );
			my @refseqs = @{ $locs{ $id }{ $ffn[0] }{ $ffn[-1] } };
			foreach my $refseq (@refseqs) {
				my $pid = $genes{$refseq}{'pid'};
				$ffn_obj->accession_number($id);
				$ffn_obj->primary_id($pid);
				$ffn_obj->display_id("gi\|$pid\|ref\|$refseq\|mol\|$id");
	
				push(@{ $genes{$refseq}{'fastant'} },$ffn_obj);
			}
		} #end each sequence
	} # end taxaid

	return;
} #end parse_ffn

#------------------------
# SUB COUNT_UNIQUE
#------------------------
sub count_unique {
  my @array = @{ shift @_ };
  my %count;
  foreach (@array) { $count{$_}++ };
  return \%count;
}


1; #loaded OK

