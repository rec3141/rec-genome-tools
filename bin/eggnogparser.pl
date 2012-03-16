#!/usr/bin/perl

use strict;
use warnings;

#1 read in taxalist
#2 parse orthgroups.mapping.txt by taxaid->cog->protein id for taxalist
#3 print table of abundances as taxaid x cog

my %eggnog;
my %groups;
my %taxa;
my @taxalist;

my $taxafile = $ARGV[0];
open(TAXA,"<$taxafile") or die "couldn't open file\n";
while(<TAXA>) {
  chomp;
  my @line = split(/\t/,$_);
  $eggnog{$line[0]}{'name'} = $line[3];
  print "$line[0]\t$line[3]\n";
  push(@taxalist,$line[0]);
}
close(TAXA);

# print "processing ", join(' ',keys(%eggnog)),"\n";
##	Protein name start_position  end_position    orthologous_group       orthologous_group description
#	9685.ENSFCAP00000011039 1       363     meNOG04000      Leucine-Rich repeat protein SHOC-2


open(EGGNOG,"</work/rec3141/gene_clustering/eggNOG/orthgroups.mapping.txt") or die "couldn't open file\n";
while(<EGGNOG>) {
  next unless ($_ =~ m/^/);
  my($id,$start,$end,$group,$desc) = split(/\t/,$_);
  my($taxaid,$proteinid) = split(/\./,$id);
  next unless defined($eggnog{$taxaid});
  $eggnog{$taxaid}{$group}++;
  ($groups{$group} = $desc) unless (defined($groups{$group}));
}
close(EGGNOG);

print "genomes\tproteins\t";
foreach my $taxon (@taxalist) {
#   print "$eggnog{$taxon}{'name'}\t";
print "$taxon\t";
}
print "eggNOG\tdescription\n";


my @printarray;

foreach my $group (sort(keys %groups)) {
  my @printline;
  my $sum;
  my $genomes;
  foreach my $taxon (@taxalist) {
    if (defined($eggnog{$taxon}{$group})) {
      $sum += $eggnog{$taxon}{$group};
      push(@printline,$eggnog{$taxon}{$group});
      $genomes += 1;
      }
    else {
      $sum += 0;
      push(@printline,'0');
      $genomes += 0;
      }
 }
 unshift(@printline,$sum);
 unshift(@printline,$genomes);
 push(@printline,$group);
 push(@printline,$groups{$group});
 push(@printarray,[ @printline ]);
#  print join("\t",@printarray),"\n";
}

#sort by $sum, then $genomes, then $group
my $eggnogcolumn = scalar(@taxalist)+1;
for my $list_ref ( sort { ($b->[1] <=> $a->[1]) || ($b->[0] <=> $a->[0]) || ($a->[$eggnogcolumn] cmp $b->[$eggnogcolumn])} @printarray ) {
    print join("\t",@$list_ref);
}