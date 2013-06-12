#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;

# this program finds the best matches between one or more binary search patterns ('10000111\n')
# and a list of observed patterns

my $patternfile = shift(@ARGV);
my $observedfile = shift(@ARGV);

open(PAT,"<$patternfile") or die("couldn't open $patternfile\n");
while(<PAT>) {
}

my %clustmap;
my $clusthash;
open(OBS,"<$obervedfile") or die("couldn't open $observedfile\n");
while(<OBS>) {
  chomp;
  next unless $_ =~ /^\d/;
  my @line = split;
  @{$clustmap{$.}} = splice(@line,0,3);
  foreach @line {push(@{$clusthash{$.}},$_) unless m/\D/;}
}

print Dumper(%clusthash);
