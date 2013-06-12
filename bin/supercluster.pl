#!/usr/bin/perl
use strict;
use warnings;
# use diagnostics;

## this program finds clusters associated with superclusters
## i.e. subsets of loose clusterings (e.g. e-value 1e-3, cutoff 0.6)
## made up of tighter clusterings (e.g. e-value 1e-30, cutoff 0.8)

unless (scalar(@ARGV)==3) {die("\nthis program takes 3 nested cluster_file files and outputs
a single file to map them all. the first input should be the less stringent
'supercluster' file, the second input should be a moderately stringent 'cluster'
file, and the last should be the very stringent 'subcluster' file. Each subsequent file
should be as or more stringent in both clustering parameters to ensure subsetability. 
Output is the mapping of each higher-level cluster to the most stringent, which is not printed (given as line number)

A 1-to-1 mapping is not guaranteed if 'seeds' were used in clustering, or if mcl clustering was used.
In this case clusters are mapped to the most frequent higher cluster

   ./supercluster.pl super_cluster_list med_cluster_list sub_cluster_list > cluster_map\n");
}

my $superfile = shift(@ARGV);
my %superc;
my $clusfile = shift(@ARGV);
my %clusters;
my $suberfile = shift(@ARGV);
my %suberc;


open(IF1,"<$clusfile") or die("couldn't open $clusfile\n");
while(<IF1>) {
  chomp;
  my @line = split;
  map { $clusters{$_} = $. } @line;
}
close(IF1);

open(IF2,"<$superfile") or die("couldn't open $superfile\n");
while(<IF2>) {
  chomp;
  my @line = split;
  map { $superc{$_} = $.} @line;
}
close(IF2);


open(IF4,"<$suberfile.r.csv") or die("couldn't open $suberfile.r.csv\n");
my @genem = <IF4>;
close(IF4);

open(IF3,"<$suberfile") or die("couldn't open $suberfile\n");
print "_supercluster\t_subcluster\t",$genem[0];
while(<IF3>) {
  chomp;
  my @line = split;
  my %seenc;
  my %seensuperc;
  foreach my $sub (@line) {
    $seenc{$clusters{$sub}}++;
    $seensuperc{$superc{$sub}}++;
  }
  my @cluster1 = sort {$seenc{$b} <=> $seenc{$a}} keys %seenc;
  my @cluster2 = sort {$seensuperc{$b} <=> $seensuperc{$a}} keys %seensuperc;
  print $cluster2[0],"\t",$cluster1[0],"\t",$genem[$.];
}
