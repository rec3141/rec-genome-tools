#!/usr/bin/perl
use strict;

# Summary:
# this script will extract all hits exceeding a particular
# threshold from a blast xml output file and store them 
# in a new file.  It includes both the refseq ID for the
# query and hit gene, as well as the E value and the length
# of the matching fraction for both the shorter and longer
# gene

# the BioPerl modules should be installed in a system location
# but we can just include them from any location they are saved.
# in the case of the former this line can be commented out.
use lib "/opt/sharcnet/bioperl/current/lib/perl5";

# use lib "/home/merz/gene_clustering/parser/BioPerl-1.6.0/";

# use the SearchIO module 

use Bio::SearchIO;

# use this listener extension to avoid allocating memory for 
# high scoring segment pairs

use Bio::SearchIO::FastHitEventBuilder;

# check to make sure we've been supplied two input arguments  

# in practice we will not cutoff hits here, so the second 
# argument can be set to the blast default of 10.0
# [REC] except the files are way too big, so cutting off at 0.001

if (scalar(@ARGV) != 3) {
  die("script takes 3 input arguments - 2 refseq identifiers and an expectation cutoff\n");
}

# check to make sure that we've been given a file with the correct
# suffix, indicating that it is a blast xml file generated as 
# described by the workflow.

my $prokdir = "/work/rec3141/gene_clustering/prok";
my $fn;
my $outf;
my $ref1 = shift(@ARGV);
my $ref2 = shift(@ARGV);
if (-e  "$prokdir/$ref1/$ref1.$ref2.xml.blastout") {
	$fn =   "$prokdir/$ref1/$ref1.$ref2.xml.blastout";
	$outf = "$prokdir/$ref1/$ref1.$ref2.hits";
} else {
        die("xml file does not exist: $prokdir/$ref1/$ref1.$ref2.xml.blastout\n");
}

# check to make sure that we've been given a realistic expectation
# cutoff (always positive, optional integer exponent) 

my $cutoff = shift(@ARGV); 
unless ($cutoff =~ /^(\d+)(\.\d*)?e?([-+]?\d+)?$/) {
        die("failed to understand cutoff: $cutoff\n");
}

# open a new file, to be overwritten (not appended), to store the
# hits in the following format (one per line):
#
# query_refseq_id hit_refseq_id significance match_length/length_shorter_gene 

open(OF, ">$outf") or die "couldn't open outfile $outf"; 

# set up the Search handle and supply the blast report

my $in = new Bio::SearchIO(-format => 'blastxml',
                           -file   => $fn);

# this will cause the parser to avoid including the additional
# HSP memory structures, should be marginally faster

$in->attach_EventHandler(Bio::SearchIO::FastHitEventBuilder->new);

# loop over queries in the report 

while( my $result = $in->next_result ) {

# loop over hits for each query

  while( my $hit = $result->next_hit ) {

# only take hits that have an E value of less than the cutoff

    if ($hit->significance <= $cutoff)  {

# check each of the high scoring pairs
# this isn't necessary, the hsp with rank 1 is always used
# as the "hit"
#      if ($hit->num_hsps  != 1) {
#         print $result->query_accession, "::", $hit->name, "::", $hit->significance, "::", $hit->length, "::", $result->query_length, "::", $hit->num_hsps, "\n"; 
#         while( my $hsp = $hit->next_hsp ) {
#            print "----", $hsp->rank, "::", $hsp->expect, "::", $hsp->hsp_length, "\n";
#         }
#      }

# calculating the match_length / gene_length(length=shorter,longer)
      my $fml_long = 0.0;
      my $fml_short = 0.0;
      my $hsp = $hit->next_hsp;
      if ($result->query_length < $hit->length) {
         $fml_short = $hsp->hsp_length / $result->query_length;
         $fml_long = $hsp->hsp_length / $hit->length;
      } else {
         $fml_long = $hsp->hsp_length / $result->query_length;
         $fml_short = $hsp->hsp_length / $hit->length;
      }

# make sure the hit name is of the right format (not sure if this
# is necessary) and extract the refseq id for the hit 

# $result->query_length  (length of query)
# $hit->length  (length of hit)

      if ($hit->name =~ m/^gi\|\d+\|(?:ref|gb)\|(.+)\|$/) {
        print OF $result->query_accession, " ", $1, " ", $hit->significance, " ", sprintf("%.3f", $fml_short), " " , sprintf("%.3f", $fml_long), "\n";
      } else {
	print OF $result->query_accession, " ", $hit->name, " ", $hit->significance, " ", sprintf("%.3f", $fml_short), " ", sprintf("%.3f", $fml_long), "\n";
#             warn("whoops! can't parse out refseq id!!", $hit->name, "\n");
      }
    }
  }
}
close(OF);
