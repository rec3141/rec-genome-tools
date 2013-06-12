#!/usr/bin/perl
use strict;
use warnings;
use diagnostics;

use lib '/home/rec3141/repo';
use RECGenomeTools qw(:All);
use Data::Dumper;
use Bio::DB::Taxonomy;

die("
   This script is used to replace a genome name or id with another available in prokaryotes.txt
   Options:
    'refseq', 'taxaname', 'replicon', 'group' (~phylum), 'subgroup' (~order), 'projectid', 'taxonid' (NCBI)
    'refseqmulti' in 'from' will split multiple refseq accessions on '-' (e.g. for headers where NC_123456-NC_234567)

    a taxalist file can be provided as 'from' to replace taxon ids (e.g. SF_000001)
			
    to convert RefSeq Genome Accession to Organism Name:     	id2id.pl refseq taxaname <file>|<stdin>
    to convert RefSeq Project ID to RefSeq Genome Accession: 	id2id.pl projectid refseq <file>|<stdin>
    to convert Organism Name to Taxonomy ID			id2id.pl taxaname taxonid <file>|<stdin>
    to convert taxalist IDs to Organism Name:			id2id.pl taxalist-file taxaname <file>|<stdin>

") if scalar(@ARGV) < 2;


my $ncbi_ref = parse_proks(); #parse prokaryotes.txt into %ncbi
my %ncbi = %{$ncbi_ref};

my @ids = qw/refseq refseqmulti taxaname replicon group subgroup projectid taxalist taxonid/;

my $from = lc(shift(@ARGV));
die("try a different 'from' identifier\n") unless grep { $from =~ m/^$_/i } @ids;

my $to = lc(shift(@ARGV));
die("try a different 'to' identifier\n") unless grep { $to =~ m/^$_/i } @ids;

my %taxahash;
if ($from eq 'taxalist') {
	die("taxalist file $from not found\n") unless (-e $from);
	my %inputs;
	$inputs{'taxalist'}=$from;
	$from='taxid';
	my $taxahash_ref = parse_taxalist(\%inputs);
	%taxahash = %{$taxahash_ref};
	
	foreach my $k1 (sort keys %taxahash) {
		foreach my $k2 (@{$taxahash{$k1}}) {
			next if defined($ncbi{$k2}{'taxid'});
			$ncbi{$k2}{'taxid'}=$k1;
			last;
		}
	}
}

my $db;
if ( ($to eq 'taxonid') || ($from eq 'taxonid') ) {
  my $idx_dir = '/work/rec3141/gene_clustering/all_genomes/idx/';
  my ( $nodefile, $namesfile ) = ( 'nodes.dmp', 'names.dmp' );
  $db = new Bio::DB::Taxonomy( -source    => 'flatfile',
				  -nodesfile => $nodefile,
				  -namesfile => $namesfile,
				  -directory => $idx_dir
  );
  die "Couldn't load taxonomic database\n\n" unless $db;

  foreach my $k1 (sort keys %ncbi) {
    $ncbi{$k1}{'taxonid'} = join(",",$db->get_taxonids( $ncbi{$k1}{'taxaname'} ) );
  }

}

my $dosplit;
if ($from eq 'refseqmulti') {
	$dosplit++;
	$from='refseq';
}

#print Dumper(%ncbi);
#exit 0;

while (<>) {
	chomp;
	my $line = $_;
	next unless $line;
	my $print = $line;
	if ($dosplit) {
		$line =~ s/\s+//g;
		my @splits = split(/-/,$line);
		foreach my $split (@splits) {
			if (exists($ncbi{$split})) {
				$line = $split;
				last;
			}
		}
	}

	foreach my $key ( grep { defined($ncbi{$_}{$from}) && defined($ncbi{$_}{$to}) } keys %ncbi) {

#		if ($line =~ s/^(\Q$ncbi{$key}{$from}\E)$/$ncbi{$key}{$to}/) {
		if ($line =~ m/^(\Q$ncbi{$key}{$from}\E)$/) {
			$print .= "\t" . $ncbi{$key}{$to};
			#last;
		}
	}

	if ( ($print eq $line) && ($to eq 'taxonid' || $from eq 'taxonid' ) ) {
	  if ($db->get_taxonids($line)) {
	    my $taxonids = join(",",$db->get_taxonids($line));
	    $print = $taxonids;
	  }
	}
	print "$print\n";
}

exit 0;
