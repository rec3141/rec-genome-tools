#!/usr/bin/perl
use lib '/work/rec3141/gene_clustering/parser/BioPerl-1.6.0';

use warnings;
use Bio::Seq;
use Bio::SeqIO;
use Bio::AlignIO;

use strict;


my $inputfilename = $ARGV[0];

my $in  = Bio::AlignIO->newFh(-file => $inputfilename , '-format' => 'fasta');
my $out = Bio::AlignIO->newFh('-format' => 'nexus', -show_symbols => '0', -show_endblock => '0');

    # World's shortest Fasta<->pfam format converter:
    print $out $_ while <$in>;
