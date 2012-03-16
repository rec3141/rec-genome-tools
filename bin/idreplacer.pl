#!/usr/bin/perl
use strict;
use warnings;
use diagnostics;

use Data::Dumper;

die("
   This script is used to replace a genome name or id with another available in lproks_1 or lproks_2
   Options: RefSeqID, ProjectID, TaxonomyID, Name, GenbankAcc, RefSeqAcc, GenBankAccList, RefSeqAccList

          convert RefSeq Ids to Organism Names:                    idreplace.pl RefSeqID Name <file>|<stdin>
          convert RefSeq Project ID to RefSeq Genome Accessions:   idreplace.pl RefSeqID RefSeqAcc <file>|<stdin>


") if scalar(@ARGV) < 2;

# ## Columns:     "RefSeq project ID"     "Project ID"    "Taxonomy ID"   "Organism Name" "Super Kingdom" "Group" "Genome Size"   "GC Content"    "Number of Chromosomes" "Number of Plasmids"    "Released date" "Modified date" "List of GenBank accessions (comma separated)"  "List of RefSeq accessions (comma separated)"   "List of publications (comma separated)"        "List of Center/Consortium (pipe separated)"
# 49725   30807   551115  'Nostoc azollae' 0708   Bacteria        Cyanobacteria   5.532   38.3    1       2       03/06/2009      11/21/2011 19:14:22     CP002059,CP002060,CP002061      NC_014248,NC_014249,NC_014250   20628610        US DOE Joint Genome Institute (JGI-PGF)|DOE Joint Genome Institute|JGI-PGF
# ## Columns:     "RefSeq project ID"     "Project ID"    "Taxonomy ID"   "Organism Name" "Super Kingdom" "Group" "Sequence availability" "RefSeq Accession"      "Genbank Accession"     "Number of Contigs"     "Number of proteins on contigs" "Genome Size"   "GC Content"    "Released date" "Center Name"   "Ceneter URL"
# 55729   33011   592010  Abiotrophia defectiva ATCC 49176        Bacteria        Firmicutes      wgs assembly    NZ_ACIN00000000 ACIN00000000    26      3291    3.4774  37.1    03/17/09        Genome Sequencing Center (GSC) at Washington University (WashU) School of Medicine      http://www.genome.wustl.edu/sub_genome_group.cgi?GROUP=3&SUB_GROUP=4

my %hash;
my %ids = (
  refseqid   => [0, 0], 
  projectid  => [1, 1],
  taxonomyid => [2, 2], 
  name       => [3, 3],
  genbankacc => [12, 8],
  refseqacc  => [13, 7],
  genbankacclist => [12, 8],
  refseqacclist  => [13, 7],
);


my $dir = '.';
system("wget -N ftp://ftp.ncbi.nih.gov/genomes/genomeprj/lproks_1.txt");
system("wget -N ftp://ftp.ncbi.nih.gov/genomes/genomeprj/lproks_1.txt");

my @dbs = ("$dir/lproks_1.txt","$dir/lproks_2.txt");

my $from = lc(shift(@ARGV));
die("try a different 'from' identifier\n") unless grep { /$from/i } keys %ids;

my $to = lc(shift(@ARGV));
die("try a different 'to' identifier\n") unless grep { /$to/i } keys %ids;

# parse databases of completed microbial genomes

for (my $i=0;$i<2;$i++) {
 my $db = $dbs[$i];
 open(DB,"<$db") or die;
 while(<DB>) {
   next if m/^#/;
   chomp;
   my @row = split(/\t/,$_);

   my @convfrom;
   # if it's an accession, check for multiple
   if ($from =~ m/acc/) {
    @convfrom = split(',',$row[$ids{$from}[$i]]);
   } else {
    push(@convfrom,$row[$ids{$from}[$i]]);
   }
   next unless ($convfrom[0] =~ m/[\w\d]+/);
#   next if $convfrom[0] eq '-';
#   next if $convfrom[0] eq ' ';

   my @convto;
   if ($to =~ m/acc/) {
     @convto = split(',',$row[$ids{$to}[$i]]);
   } else {
     @convto = $row[$ids{$to}[$i]];
   }
   next unless ($convto[0] =~ m/[\w\d]+/);
#   next if $convto[0] eq '-';
#   next if $convto[0] eq ' ';

   foreach my $cfrom (@convfrom) {
     if ($to =~ m/list/) {
       $hash{$cfrom} = join(' ',@convto);
     } else {
       $hash{$cfrom} = $convto[0];
     }
   }
 }
 close(DB);
}

while (<>) {
  foreach my $key (keys %hash) {
   s/\Q$key\E/$hash{$key}/;
  } 
   print;
}

#print Dumper(%hash);

