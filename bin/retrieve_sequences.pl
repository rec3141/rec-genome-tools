#!/usr/bin/perl -w
use lib "/opt/sharcnet/bioperl/1.6.1/lib/perl5";
use Bio::Perl;

use Bio::Seq;
use Bio::SeqIO;
use Bio::AlignIO;
use Data::Dumper;
use strict;
use diagnostics;

sub get_evalues;
sub make_hash;
sub count_unique;

# Summary:

# This script can be called to retrieve gene sequences and other characteristics belonging to any 
# specific cluster(s)

# clusters are currently specified by their index (line number) in the 
# cluster list file, separated by commas

# a cluster number 0 indicates all clusters in file

# this needs to be generalized to N clusters.
$|++;

if (scalar(@ARGV) != 6) {
  die("\nscript takes 6 inputs:\n
    a parameter to output for each pid, one of 'fastant', 'fastaaa', 'cog', 'gene', 'func', 'taxid', 'projectid', 'taxaname', 'id', 'replicon', 'loc', 'desc', or 'evalue'
    maximum number of values to output per genome (for sequences) or cluster (for others); 0 for all (can be slow)
    a cluster index in the file, 0 for all, comma-separated for multiple
    a filename containing clusters of protein ids
    a taxalist file
    a taxalist id to limit output clusters to
    eg. ./retrieve_sequences.pl cog 42 NC_000123.cluster_list taxalist-Archaea AR_110001

    WARNING: 'fastant' sequences from DRAFT genomes with many contigs/scaffolds may not be available!\n");
}

### OUTPUT
### a parameter to output for each pid, one of 'fastant', 'fastaaa', 'cog', 'gene', 'func', 'taxid', 'projectid', 'taxaname', 'id', 'replicon', 'loc', 'desc', or 'evalue'

my @paramchoice = qw|fastant fastaaa cog gene func taxid projectid taxaname id replicon loc desc evalue|;
my $output_parameter = shift(@ARGV);
my $output_format;
if (grep $_ eq $output_parameter, @paramchoice) {
  if ($output_parameter =~ s/nexus//) {
    $output_format = 'nexus';
  } elsif ($output_parameter =~ s/fasta//) {
    $output_format = 'fasta';
  } else {$output_format = 'NA'}
} else {
    die("$0 failed to understand input argument 1: ", $output_parameter, "\n");
}

### MAXVALS
###
my $maxvals = shift(@ARGV);
die("$0 failed to understand input argument 2: ", $maxvals, "\n") unless $maxvals =~ m/^\d+/;

### CLUSTERS
### a cluster index in the file, 0 for all, comma-separated for multiple
### read in cluster index (just an integer or CSV list of integers)
my %clus_hash;
my $clus_request = shift(@ARGV);
if ($clus_request =~ /\d+/) {
  if ($clus_request =~ /,/) {
    map {$clus_hash{$_}++} split(/,/,$clus_request);
  } else {
    $clus_hash{$clus_request}++;
  }
} else {
  die("$0 failed to understand input argument 3: ", $clus_request, "\n");
}

### CLUSTER FILE
### a filename containing clusters of protein ids
my $clus_file = shift(@ARGV);
die("$0 failed to understand input argument 4: ", $clus_file, "\n") unless ( -e $clus_file);
open(IF1, "< $clus_file")  || die("couldn't open file $clus_file: $!\n");
if ($clus_hash{0}) {while (<IF1>) {}; delete $clus_hash{0}; map {$clus_hash{$_}++} (1..$.);}
close(IF1);

my $clus_num = scalar(keys %clus_hash);
my $clus_found = 0;
my %prot_seen; #proteins seen
open(IF1, "< $clus_file")  || die("couldn't open file $clus_file: $!\n");
while (<IF1>) {
#   print "\n$.\t";
  next unless defined($clus_hash{$.});
  $clus_found++;
  chomp;
  my @protline = split(/\ /,$_);
  if (($maxvals>0) && (scalar(@protline)>$maxvals) && ($output_format eq 'NA')) {
    my $i=0;
    while ($i<$maxvals) {
      my $tmp = $protline[rand @protline];
      unless ($prot_seen{$tmp}) {$prot_seen{$tmp}++; $i++;}
    }
  } else { map {$prot_seen{$_}++} @protline}

  last if $clus_num==$clus_found;
}
die("couldn't find all requested clusters\n") if ($clus_found < $clus_num);
close(IF1);

### TAXALIST
### a taxalist file

my $hit;
my $taxalist = shift(@ARGV);
die("$0 failed to understand input argument 5: ", $taxalist, "\n") unless ( -e $taxalist);

# read in the taxalist to get list of genomes
my %taxahash;
open(IF, "<$taxalist") || die("couldn't open file $taxalist\n");
while (<IF>) {
  chomp;
  my @taxa = split('\s',$_);
  my $tag = shift @taxa;
  next unless $tag;
  if (scalar(@taxa) > 0) {
    @{$taxahash{$tag}} = @taxa;
  } else {
    @{$taxahash{$tag}} = $tag;
#     $taxahash{$tag} = $tag;
  }
}
close(IF);

### TAXALIST ID
### a taxalist id to limit clusters to
my @ids=();
my %all_ids=();
my $id = shift(@ARGV);
if ($id =~ /^\D{2}\_\d+$/) {
    if (defined(@{$taxahash{$id}})) {
      foreach my $moreid (@{$taxahash{$id}}) {
	if (defined(@{$taxahash{$moreid}})) {
	  push(@ids,[ @{$taxahash{$moreid}} ]);
	  map {$all_ids{$_}++} @{$taxahash{$moreid}};
	} else {push(@ids,[ $moreid ]); $all_ids{$moreid}++;}
      }
    } else {
      die("id $id not found in $0\n");
    }
} else {
    die("$0 failed to understand input argument 6: ", $id, "\n");
}

my @faadirs = qw|/work/rec3141/gene_clustering/all_genomes/faa /work/rec3141/gene_clustering/all_genomes/scaffold/faa /work/rec3141/gene_clustering/all_genomes/contig/faa /work/rec3141/gene_clustering/all_genomes/contig/faa|;
my %accgenomes=();
# print "reading accession lists...";
foreach my $faadir (@faadirs) {
  open(IF,"<$faadir/accession_list") or die "couldn't open file $faadir/accession_list\n";
  while (<IF>) {
    chomp;
    my ($accgenome,$accgene) = split(/\t/,$_);
    next unless exists($all_ids{$accgenome});
    $accgenomes{$accgenome}{$accgene}++ if $prot_seen{$accgene};
  }
  close(IF);
}

my @intersection = ();
my %count = ();
foreach my $k1 (keys %accgenomes,keys %all_ids) { 
  $count{$k1}++;
}

foreach my $k2 (keys %count) {
    $count{$k2} > 1 ? push(@intersection,$k2) : undef;
}



###############
# get info on genomes
# accessions are unique per molecule
# taxids are shared by molecules of the same strain in different projects
# projectids are shared by all molecules of the same project

my %ncbihash;
my $ncbifile = "/work/rec3141/gene_clustering/all_genomes/summary.txt";
open( IF, "<$ncbifile" ) or die("couldn't open file $ncbifile\n");
while (<IF>) {
  my ( $access, $genbankacc, $length, $taxid, $projectid, $taxaname, $replicon ) = split( "\t", $_ );

  next if $access =~ m/Accession/;
  my ($accession) = split(/\./,$access);
  $ncbihash{$accession}{'taxaname'} = $taxaname;
  $ncbihash{$accession}{'replicon'} = $replicon;
  $ncbihash{$accession}{'taxid'}    = $taxid;
  $ncbihash{$accession}{'projectid'} = $projectid;
}
close IF;

### MAKE THE PTT HASH
my %ptthash;
%ptthash = make_hash(@intersection);

### OUTPUT THE HASH
#loop over requested clusters

open(IF1, "< $clus_file")  || die("couldn't open file $clus_file: $!\n");
  while (<IF1>) {
  next unless $clus_hash{$.};
  chomp;

  # foreach my $cluster (@th1) {
  my @th2 = split(/\ /,$_);
  my @output = ();

  #output appropriate to request
  if ($output_parameter eq 'nt' || $output_parameter eq 'aa') {
    foreach (@th2) {
      my $sequence = Bio::SeqIO->new(-fh => \*STDOUT, '-format' => 'fasta' );
      $sequence->write_seq($ptthash{$_}{$output_parameter}) if $ptthash{$_}{$output_parameter};
    }

  } elsif ($output_parameter eq 'evalue') {
    my @evalue_ids;
    foreach (@th2) {
      push(@evalue_ids,$ptthash{$_}{'id'});
    }
    my @evalues = get_evalues(\@evalue_ids,\@th2);
    print "\n";
  } else {
    foreach (@th2) {
      push(@output,$ptthash{$_}{$output_parameter}) if ($ptthash{$_}{$output_parameter});
    }
    print join("",@output) . "\n";
  }  
}

exit 0;

###----------------------
### SUBROUTINES 
###----------------------
sub get_evalues {
  my @evals=();
  my @all_ids = @{$_[0]}; #evalue_ids
  my %seen = ();
  @ids = grep { ! $seen{$_} ++ } @all_ids;
  my @proteins = @{$_[1]}; #protein ids

  foreach my $subject (@ids) {
    $subject =~ s/[\s\n\r\f]//g;
    chomp $subject;
    foreach my $query (@ids) {
      $query =~ s/[\s\n\r\f]//g;
      chomp $query;

      my %hash = ();
      open(IF,"<$subject/$subject.$query.hits") or die("couldn't open hits file $subject/$subject.$query.hits\n");
      while(<IF>) {
	chomp;
	my @protline = split(/\s/,$_);
	$hash{$protline[0]}{$protline[1]}=$protline[2];
      }
      close(IF);

      my @w = keys(%hash);
      foreach my $e (keys %hash) {
	push (@w,keys %{$hash{$e}});
      }
      my @x = @proteins;

      my %union = ();
      my %isect = ();
      foreach my $f (@w) { $union{$f} = 1 }
      foreach my $g (@x) {
	if ( $union{$g} ) { $isect{$g} = 1 }
	$union{$g} = 1;
      }

      foreach my $i (keys %isect) {
	foreach my $j (keys %isect) {
	  if (exists($hash{$i}{$j})) {print join("\t",$i,$j,$hash{$i}{$j}),"\n";}
	}
      }
    } #end foreach my query
 } #end foreach my subject

 return(@evals);
} #end sub

#----------------------

sub make_hash {
  my @ids = @_;
  my @moldirs = qw|/work/rec3141/gene_clustering/all_genomes /work/rec3141/gene_clustering/all_genomes/scaffold /work/rec3141/gene_clustering/all_genomes/contig|;
    foreach my $id (@ids) {
      my $dir = undef;
      foreach my $moldir (@moldirs) {
	next unless (( -e "$moldir/faa/$id.faa"));
	next unless (( -e "$moldir/ptt/$id.ptt"));
	next unless (( -e "$moldir/ffn/$id.ffn"));
	$dir = $moldir;
	last;
      }
      unless (defined($dir)) {
       warn ("missing files for $id\n");
       next;
      }
#       warn "id: $id\tdir: $dir\n";

      my %faahash=(); #temporary for each id
      my $faaio = Bio::SeqIO->new(-file => "$dir/faa/$id.faa",-format=>"fasta");

      # >gi|168697871|ref|ZP_02730148.1| hypothetical protein GobsU_00005 [Gemmata obscuriglobus UQM 2246]
      while (my $faa_obj = $faaio->next_seq) {
	my (@display_ids) = split('\|',$faa_obj->display_id);
	chomp @display_ids;
	next unless defined($accgenomes{$id}{$display_ids[3]});
	delete($accgenomes{$id}{$display_ids[3]});
	# $ptthash{'ZP_08308501.1'}{'seq_obj'} = {seq_obj}
	$ptthash{$display_ids[3]}{'aa'} = $faa_obj;
	$ptthash{$display_ids[3]}{'aa'}->display_id($ptthash{$display_ids[3]}{'aa'}->primary_id . "mol\|$id");

	# $faahash{'330444846'}{'ref'} = 'ZP_08308501.1'
	$faahash{$display_ids[1]}{'ref'} = $display_ids[3];

	last if scalar(values %{$accgenomes{$id}}) == 0;
      }

      #get nucleotide sequences if so requested
      my %ffnhash=(); #temporary for each id
      if ($output_parameter eq 'nt') {
	my $ffn;
	my $ffnio = Bio::SeqIO->new(-file => "$dir/ffn/$id.ffn",-format=>"fasta");
	while (my $ffn_obj = $ffnio->next_seq) {
	  my ($refseq,$location);
	  if ($ffn_obj->display_id =~ /^gi/) {
	      (undef,undef,undef,$refseq,$location) = split('\|',$ffn_obj->display_id);
	  } elsif ($ffn_obj->display_id =~ /^ref/) {
	      (undef,$refseq,$location) = split('\|',$ffn_obj->display_id);
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
	  push(@{$ffnhash{ $ffn[0] }{ $ffn[-1] }},$ffn_obj);
	}
      } # end .ffn

      # process ptt file
      open(PTT,"<$dir/ptt/$id.ptt") || die("couldn't open file $dir/ptt/$id.ptt\n");
      PTTLINE: while (<PTT>) {
	chomp;
	my $pttline = $_;
# 	print $pttline,"\n";
	if ($. == 1) {
	  # Enterococcus faecium 1,231,501 genomic scaffold supercont1.16 - 1..1230
	  my @line = split(/\s+/,$pttline);
	  my ($b1) = grep { $line[$_] =~ m/(?:whole|genomic|plasmid|chromosome)/ } 0..$#line;
	  if ($b1) {
	    $ncbihash{$id}{'taxaname'} = join(' ',@line[0..($b1-1)]) unless ($ncbihash{$id}{'taxaname'});
	    $ncbihash{$id}{'replicon'} = join(' ',@line[$b1..-1]) unless ($ncbihash{$id}{'replicon'});
	  } else {
	    $ncbihash{$id}{'taxaname'} = 'NA' unless ($ncbihash{$id}{'taxaname'});
	    $ncbihash{$id}{'replicon'} = 'NA' unless ($ncbihash{$id}{'replicon'});
	  }
	  $ncbihash{$id}{'taxid'}    = 'NA' unless ($ncbihash{$id}{'taxid'});
	  $ncbihash{$id}{'projectid'} = 'NA' unless ($ncbihash{$id}{'projectid'});
	  next PTTLINE;
	}
	next unless ($pttline =~ m/^\d+\.\.\d+/);

	my ($pttloc, undef, undef, $pttpid, $pttgene, undef, undef, $pttcog, $pttdesc) = split(/\t/,$_);
	next unless (defined($faahash{$pttpid}{'ref'}));
	my $ref = $faahash{$pttpid}{'ref'};

	$ptthash{$ref}{'pid'} = $pttpid . "\t";
	$ptthash{$ref}{'gene'} = $pttgene . "\t";
	$ptthash{$ref}{'desc'} = $pttdesc . "\t";
	$ptthash{$ref}{'id'} = "\n" . $id;

	$ncbihash{$id}{'taxid'}     ? $ptthash{$ref}{'taxid'} = "\n" . $ncbihash{$id}{'taxid'}		: $ptthash{$ref}{'taxid'} = "\nNA";
	$ncbihash{$id}{'taxaname'}  ? $ptthash{$ref}{'taxaname'} = "\n" . $ncbihash{$id}{'taxaname'}	: $ptthash{$ref}{'taxaname'} = "\nNA";
	$ncbihash{$id}{'replicon'}  ? $ptthash{$ref}{'replicon'} = $ncbihash{$id}{'replicon'} . "\t"	: $ptthash{$ref}{'replicon'} = "NA\t";
	$ncbihash{$id}{'projectid'} ? $ptthash{$ref}{'projectid'} = "\n" . $ncbihash{$id}{'projectid'}	: $ptthash{$ref}{'projectid'} = "\nNA";

	my($locs,$loce) = split(/\.\./,$pttloc);
	$ptthash{$ref}{'locs'} = $locs;
	$ptthash{$ref}{'loce'} = $loce;
# 	($ptthash{$ref}{'locs'},$ptthash{$ref}{'loce'}) = split(/\.\./,$pttloc);
	$ptthash{$ref}{'loc'} = "\n" . $id . "\t" . $locs . "-" . $loce . "\t" . $ptthash{$ref}{'desc'};
	if ($output_parameter eq 'nt') {
	    if (defined($ffnhash{ $locs }{ $loce })) {
	      if (scalar(@{$ffnhash{ $locs }{ $loce }})>1) {
		$ptthash{$ref}{'nt'} = Bio::Seq->new(-seq => 'N'x128, -format => "fasta", -accession_number => $id, -primary_id => $pttpid, -display_id => "gi\|$pttpid\|ref\|$ref\|mol\|$id");
	      } else {
		$ptthash{$ref}{'nt'} = shift(@{$ffnhash{$locs}{$loce}});
     		$ptthash{$ref}{'nt'}->accession_number($id);
		$ptthash{$ref}{'nt'}->primary_id($pttpid);
		$ptthash{$ref}{'nt'}->display_id("gi\|$pttpid\|ref\|$ref\|mol\|$id");
	      }

	  } else {
	    warn "can't find ptt location in ffn: $id\t$pttpid\t$ref\t$locs\t$loce\n";
	  }
	} 

	if ($pttcog =~ m/(COG\d*)([A-Z]*)/) {
	  $ptthash{$ref}{'cog'} = $1 . "\t";
	  $ptthash{$ref}{'func'} = $2 . "\t";
# 	} elsif ($pttcog =~ m/[\s-]*/) {
# 	  $ptthash{$ref}{'cog'} = "COG0000\t";
# 	  $ptthash{$ref}{'func'} = "X\t";
	} else {
	  $ptthash{$ref}{'cog'} = "-\t";
	  $ptthash{$ref}{'func'} = "-\t";
	}
      }
      close(PTT);
#     print Dumper(%ptthash);
    } #end foreach my $id
   
  return %ptthash;
} # end sub

#----------------------

sub count_unique {
  my @array = @_;
  my %count;
  map { $count{$_}++ } @array;
  return %count;
}

#####################################
# Here is a process to get the nucleotide sequence given the RefSeq ID:

#RefSeq->PID(.faa)->Loc(.ptt)->Seq(.ffn)

#[merz@wha780 program]$ head -5 NC_000964.ptt 
#Bacillus subtilis subsp. subtilis str. 168, complete genome - 1..4215606
#4176 proteins
#Location        Strand  Length  PID     Gene    Synonym Code    COG     Product
#410..1750       +       446     16077069        dnaA    BSU00010        -       COG0593L        chromosomal replication initiation protein
#1939..3075      +       378     16077070        dnaN    BSU00020        -       COG0592L        DNA polymerase III subunit beta

#[merz@wha780 program]$ head -5 NC_000964.faa
#>gi|16077069|ref|NP_387882.1| chromosomal replication initiation protein [Bacillus subtilis subsp. subtilis str. 168]
#MENILDLWNQALAQIEKKLSKPSFETWMKSTKAHSLQGDTLTITAPNEFARDWLESRYLHLIADTIYELT
#GEELSIKFVIPQNQDVEDFMPKPQVKKAVKEDTSDFPQNMLNPKYTFDTFVIGSGNRFAHAASLAVAEAP
#AKAYNPLFIYGGVGLGKTHLMHAIGHYVIDHNPSAKVVYLSSEKFTNEFINSIRDNKAVDFRNRYRNVDV
#LLIDDIQFLAGKEQTQEEFFHTFNTLHEESKQIVISSDRPPKEIPTLEDRLRSRFEWGLITDITPPDLET

#[merz@wha780 program]$ head -5 NC_000964.ffn
#>ref|NC_000964.3|:410-1750
#ATGGAAAATATATTAGACCTGTGGAACCAAGCCCTTGCTCAAATCGAAAAAAAGTTGAGCAAACCGAGTT
#TTGAGACTTGGATGAAGTCAACCAAAGCCCACTCACTGCAAGGCGATACATTAACAATCACGGCTCCCAA
