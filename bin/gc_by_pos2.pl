#!/usr/bin/perl

use strict;
use warnings;
use diagnostics;
use Data::Dumper;
use List::Util qw(sum);

sub countgc;
sub printpos;

my $filename = $ARGV[0];
my $bygene = 0; #1 if for each gene, 0 if for whole genome
my $sequence;
my @position;
my %pos;
my $header;
while (<>) {
  chomp;
  if ($_ =~ m/^>/) {
    if ($sequence) {
      if (length($sequence)%3==0) {
	countgc();
	printpos($header) if ($bygene);
      }
    }
    $sequence = '';
    $header=$_;
    chomp $header;
  } else {
    $sequence .= $_; 
  }
}

countgc();
printpos($filename);

sub countgc {
  $sequence =~ s/[^ATGC]/0/g;
  $sequence =~ tr/ATGC/0011/;
  my $codons = length($sequence)/3;
  my $mask1 = '100' x $codons;
  my $mask2 = '010' x $codons;
  my $mask3 = '001' x $codons;
  $pos{0}{'GC'} += ($mask1 & $sequence =~ tr/1//);
  $pos{1}{'GC'} += ($mask2 & $sequence =~ tr/1//);
  $pos{2}{'GC'} += ($mask3 & $sequence =~ tr/1//);
  $pos{0}{'AT'} += ($mask1 & $sequence =~ tr/0//);
  $pos{1}{'AT'} += ($mask2 & $sequence =~ tr/0//);
  $pos{2}{'AT'} += ($mask3 & $sequence =~ tr/0//);

  return 1;
}

sub printpos {
  my $head = shift;
  my @gc;
  @gc = map { $gc[$_] = $pos{$_}{'GC'} } 0..2;
  my @at;
  @at = map { $at[$_] = $pos{$_}{'AT'} } 0..2;

  print $gc[0]/($gc[0]+$at[0]),
  "\t",$gc[1]/($gc[1]+$at[1]),
  "\t",$gc[2]/($gc[2]+$at[2]),
  "\t",$head,
  "\n";

  %pos = ();

  return 1;
}
