#!/usr/bin/perl
use strict;
use warnings;
use diagnostics;

use Data::Dumper;

die
"\nthis program takes 3 inputs:
1) a path to a directory of .fasta alignments to concatenate 
2) a taxalist id
3) a taxalist

e.g. fastaconcat.pl ./trees BL_110001 taxalist-Bacilli-with-plasmids

and outputs:
1) a concatenated fasta file named with RefSeq IDs" 
if ( scalar(@ARGV) != 3 );

# read in the taxalist to get list of genomes, to check if id is a plasmid
my $taxalist = $ARGV[2];
my %genomehash;
my %taxahash;
my %plasmids;
open(IF, "<$taxalist") || die("couldn't open file $taxalist\n");
while (<IF>) {
  chomp;
  my @taxa = split(/\s+/,$_);
  my $tag = shift @taxa;
  next unless $tag;
    if (scalar(@taxa) > 0) {
    @{$taxahash{$tag}} = @taxa;
  } else {
    @{$taxahash{$tag}} = $tag;
    $genomehash{$tag} = undef;
  }

  my $tag2 = shift @taxa;
  if ($tag =~ m/^\D{2}\_0+/) {
	  foreach my $plasmid (@taxa) {
		  $plasmids{$plasmid} = $tag2;
		  $genomehash{$tag2}=undef;
	  }
  }
}
close(IF);

my @tmpids;
my @ids;
my $id;
if ($ARGV[1] =~ /(\D{2}\_\d+)/) {
    $id = $1;
    die("$0 id $id not found\n") unless defined(@{$taxahash{$id}});
    @tmpids = @{$taxahash{$id}};
} else {
    die("$0 failed to understand input argument 1: ", $ARGV[0], "\n");
}


# to de-code the list
 foreach my $moreid (@tmpids) {
  if (defined(@{$taxahash{$moreid}})) {
	push(@ids, @{$taxahash{$moreid}});
	map {$genomehash{$_} = undef} @ids;
  }
  else {
    push(@ids, $moreid );
    $genomehash{$moreid}=undef;

    }
 }


my @faadirs = qw|/work/rec3141/gene_clustering/all_genomes/faa /work/rec3141/gene_clustering/all_genomes/scaffold/faa /work/rec3141/gene_clustering/all_genomes/contig/faa /work/rec3141/gene_clustering/all_genomes/metagenomes/faa|;
my %genelist;
print "reading gene list...";
foreach my $faadir (@faadirs) {
  open(IF,"<$faadir/accession_list") or die "couldn't open file $faadir/accession_list\n";
  print "$faadir...";
  while (<IF>) {
  chomp;
  my ($genome,$gene) = split(/\t/,$_);
  $genelist{$gene}=$genome if exists($genomehash{$genome});
  }
  close(IF);
}
print "done\n";

my %hash;
my $dir = $ARGV[0];
print "DIR: $dir\n";
#for each fasta file
my @files = <$dir/*.fasta>;
#print join("\n",@files);
# >gi|152976215|ref|YP_001375732.1| ribonuclease III [Bacillus cereus subsp. cytotoxis NVH 391-98]
foreach my $file (@files) {
open(FILE,"<$file") or die;
print "$file\n";
  my $name;
  while(<FILE>) {
    chomp;
    if (m/^>/) { 
		m/^>gi\|\d+\|(ref|gb)\|(.*)\|mol/;
		if (defined($genelist{$2})) {
			if (defined($plasmids{$genelist{$2}})) {
				$name = $plasmids{$genelist{$2}};
			} else {$name= $genelist{$2}}
		} else {warn "$0 couldn't find accession $2 in $file\n\tline=$_\n";}
    } elsif (m/^$/) {
    	$name="";
    } else {
    	$hash{$name} .= $_;
    }
  }
close(FILE);
}
# print Dumper(%hash);

my %union;
my %isect;
foreach my $e (keys %hash, @ids) { 
	$union{$e}++ && $isect{$e}++ 
};

open(OUT,">$dir/$id.concat.faa");
foreach my $k (keys %isect) {
print OUT ">$k\n$hash{$k}\n";
}
close(OUT);
