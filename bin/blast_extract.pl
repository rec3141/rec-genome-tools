#!/usr/bin/perl

use File::Copy;
use strict;

# Summary:
# this script will extract all hits exceeding a particular
# threshold from a blast xml output and store them 
# in a new file.  It includes both the refseq ID for the
# query and hit gene, as well as the E value and the length
# of the matching fraction for both the shorter and longer
# gene

# use the SearchIO module 
use Bio::SearchIO;

# use this listener extension to avoid allocating memory for 
# high scoring segment pairs

use Bio::SearchIO::FastHitEventBuilder;



# check to make sure we've been supplied two input arguments  

if (scalar(@ARGV) != 2) {
  die("script takes 2 input arguments - 2 RefSeq IDs\n");
}

my $subject = shift(@ARGV);
die "wrong input argument 1, should be a RefSeq ID for subject genome\n" unless $subject =~ m/^[A-Z]{2}\_[\w\d]+/;
my $query   = shift(@ARGV);
die "wrong input argument 2, should be a RefSeq ID for query genome\n" unless $query =~ m/^[A-Z]{2}\_[\w\d]+/;

my $basedir;
if ( -e "/work/rec3141/gene_clustering" ) { $basedir = "/work" }
elsif ( -e "/gwork/rec3141/gene_clustering" ) {$basedir = "/gwork" }
else {die "could not find working directory\n"}
my $workdir="$basedir/rec3141/gene_clustering/all_genomes";
my $prokdir="$basedir/rec3141/gene_clustering/prok";
my @dirlist = ("faa","scaffold/faa","contig/faa","metagenomes/faa");

my $sub_dir;
foreach my $dir (@dirlist) {
#	print "$workdir/$dir/$subject.faa\n";
	if ( -f "$workdir/$dir/$subject.faa") {
		$sub_dir = $dir;
		last;
	}
}
if ($sub_dir) {
	print "found $sub_dir/$subject.faa\n";
} else {die "couldn't find $subject.faa: $!\n"};

my $query_dir;
foreach my $dir (@dirlist) {
# 	print "$workdir/$dir/$query.faa\n";
	if ( -f "$workdir/$dir/$query.faa") {
		$query_dir = $dir;
		last;
	}
}

if ($query_dir) {
	print "found $query_dir/$query.faa\n";
} else {die "couldn't find $query.faa: $!\n"};


mkdir("$prokdir/$subject") unless ( -d "$prokdir/$subject");


 # do the blasting, read from STDIN to avoid writing xml file
my $fh;
my $cutoff = 0.001;
my $blastp   = "/home/rec3141/program/bin/blastp";
my $options = "-db $workdir/$sub_dir/$subject.faa -query $workdir/$query_dir/$query.faa -outfmt 5 -comp_based_stats 2 -evalue $cutoff -parse_deflines";

my $command = "$blastp $options";
 
print "blast $subject $query begin: $command\n";

open $fh,"$command |" || die("cannot run $command: $!\n");

# set up the Search handle and supply the blast report

my $searchio  = Bio::SearchIO->new(-format => 'blastxml', -fh => $fh);

# this will cause the parser to avoid including the additional
# HSP memory structures, should be marginally faster

$searchio->attach_EventHandler(Bio::SearchIO::FastHitEventBuilder->new);

# open a new file, to be overwritten (not appended), to store the
# hits in the following format (one per line):
#
# query_refseq_id hit_refseq_id significance match_length/length_shorter_gene 
# my $outf = "$prokdir/$subject/$subject.$query.hits.tmp.$$";
my $outf = "/scratch/rec3141/$subject.$query.hits.tmp.$$";

open(OF, ">$outf") or die "couldn't open outfile $outf"; 

# loop over queries in the report 

while( my $result = $searchio->next_result ) {
# loop over hits for each query

  while( my $hit = $result->next_hit ) {

# only take hits that have an E value of less than the cutoff

    if ($hit->significance <= $cutoff)  {

# calculating the match_length / gene_length(length=shorter,longer)
      my $fml_long = 0.0;
      my $fml_short = 0.0;
      my $hsp = $hit->next_hsp;
      if ($result->query_length < $hit->length) {
         $fml_short = ($hsp->hsp_length - $hsp->gaps)/ $result->query_length;
         $fml_long = ($hsp->hsp_length - $hsp->gaps) / $hit->length;
      } else {
         $fml_long = ($hsp->hsp_length - $hsp->gaps) / $result->query_length;
         $fml_short = ($hsp->hsp_length - $hsp->gaps) / $hit->length;
      }

# make sure the hit name is of the right format (not sure if this
# is necessary) and extract the refseq id for the hit 

      if ($hit->name =~ m/^gi\|\d+\|(?:ref|gb)\|(.+)\|$/) {
        print OF $result->query_accession, " ", $1, " ", $hit->significance, " ", sprintf("%.3f", $fml_short), " " , sprintf("%.3f", $fml_long), "\n";
      } else {
	print OF $result->query_accession, " ", $hit->name, " ", $hit->significance, " ", sprintf("%.3f", $fml_short), " ", sprintf("%.3f", $fml_long), "\n";
#             warn("whoops! can't parse out refseq id!!", $hit->name, "\n");
      }
    }
  }
}
print "blast $subject $query done\n";
close(OF);

# move("$prokdir/$subject/$subject.$query.hits.tmp.$$","$prokdir/$subject/$subject.$query.hits");
move("/scratch/rec3141/$subject.$query.hits.tmp.$$","/scratch/rec3141/$subject.$query.hits");

exit 0;
