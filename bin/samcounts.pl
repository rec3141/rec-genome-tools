#!/usr/bin/perl -w
use strict;

my %hash;
my %genomes;
my %metagenomes;

my $dir = shift(@ARGV);
$dir ? opendir(DIR, "$dir") : opendir(DIR, ".");
my @files = readdir(DIR);
closedir(DIR);

foreach my $file (@files) {
  next unless $file =~ m/\.sam$/;
  print $file,"\n";
  my ($genome, $metagenome) = split(/\./,$file);
  open(FILE,"<$file");
  while(<FILE>) {
    next unless $_ =~ m/gi/;
    $hash{$genome}{$metagenome}++;
    $genomes{$genome}=1;
    $metagenomes{$metagenome}=1;
  }
  close(FILE);
}

open(OUT,">samcounts.csv");
print OUT "genomes\t",join("\t",keys %metagenomes),"\n";
foreach my $k1 (keys %genomes) {
  print OUT "$k1";
  foreach my $k2 (keys %metagenomes) {
    exists($hash{$k1}{$k2}) ? print OUT "\t",$hash{$k1}{$k2}-1 : print OUT "\t0";
  }
  print OUT "\n";
}
close(OUT);
