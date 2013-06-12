#!/usr/bin/perl -w
use strict;

use Bio::DB::Taxonomy;
#use Data::Dumper;

# use warnings;
 use diagnostics;

#print("\nthis program updates the taxonomy database\n\n");

my $idx_dir = '/work/rec3141/gene_clustering/all_genomes/idx';

my ( $nodefile, $namesfile ) = ( "nodes.dmp", "names.dmp" );

my $db = new Bio::DB::Taxonomy( -source    => 'flatfile',
                                -nodesfile => $nodefile,
                                -namesfile => $namesfile,
                                -directory => $idx_dir
);
#print("done");
