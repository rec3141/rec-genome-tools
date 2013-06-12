#!/usr/bin/perl

use strict;

# read in all .faa files
# concatenate a list of genome\tgene

my @dirs = qw|/work/rec3141/gene_clustering/all_genomes/faa /work/rec3141/gene_clustering/all_genomes/contig/faa /work/rec3141/gene_clustering/all_genomes/scaffold/faa|;

foreach my $dir (@dirs) {
 unlink("$dir/accession_index");

 open(OUT,">>$dir/accession_index") or die "couldn't open output file";
 my @files = glob("$dir/*.faa");

foreach my $file (@files) {
# print "$file\n";
 open(IF,"<",$file) or die "couldn't open input file";
 my @path = split(/\//,$file);
 my ($refseq,undef) = split(/\./,$path[-1]);
 print OUT $refseq;
 while (<IF>) {
   next unless ($_ =~ s/^\>//);
   chomp;
    my (undef,undef,undef,$acc) = split(/\|/,$_);
    my ($acc2,undef) = split(/\s/,$_);
#    if (defined($acc)) {
       print OUT " $acc";
#    } elsif (defined($acc2)) {
#       print OUT " $acc2";
#    } else {
#       print OUT " $refseq.UNKNOWN"
#    }
#     if ($_ =~ m/^gi\|\d+\|(?:ref|gb)\|(.+)\|/) {
# #      print OUT "$refseq\t$1\n";
#        print OUT " $1";
#     } else {
# #      print OUT "$refseq\t$_\n";
#        print OUT " $_";
#     }
  }
  print OUT "\n";
}

close(OUT)
}
