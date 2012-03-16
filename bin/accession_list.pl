#!/usr/bin/perl

use strict;

# read in all .faa files
# concatenate a list of genome\tgene

my @dirs = qw|/work/rec3141/gene_clustering/all_genomes/faa /work/rec3141/gene_clustering/all_genomes/contig/faa /work/rec3141/gene_clustering/all_genomes/scaffold/faa /workk/rec3141/gene_clustering/all_genomes/metagenomes/faa|;

foreach my $dir (@dirs) {
unlink("$dir/accession_list");

open(OUT,">>$dir/accession_list") or die "couldn't open output file";
my @files = glob("$dir/*.faa");

foreach my $file (@files) {
print "$file\n";
 my @path = split(/\//,$file);
 my ($refseq,undef) = split(/\./,$path[-1]);
 open(IF,"<$file") or die "couldn't open input file";
 while (<IF>) {
   chomp;
   next unless ($_ =~ s/^\>//);
    if ($_ =~ m/^gi\|\d+\|(?:ref|gb)\|(.+)\|$/) {
     print OUT "$refseq\t$1\n";
    } else {
     print OUT "$refseq\t$_\n";
    }
  }
}

close(OUT)
}
