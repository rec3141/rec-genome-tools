#!/usr/bin/perl -w
use strict;
use diagnostics;
use Archive::Tar;
use File::Slurp;
print "this program takes one or two blastlists and sets up
all-against-all blast jobs for all of the taxa
in the lists, on the cluster, e.g.

~/repo/blastlist_process.pl blastlist-archaea

creates .run files containg the blast commands, 
which can be run by ./newblasts.run
or in parallel using quickblast.sh

change: if 2 blastlists provided, will force re-blast of all of first against all of second
change: now puts all files in /prok/ and softlinks to them
change: now extracts files from archive tar.gz files when they exist
        the code is: (e)xists on filesystem, archived in tg(z), or (n)ew

\n";

# change: now splits blasts into batches of 50 per file
#to change: when finding files, only identify if it is archived, don't extract it (let signifier do that)
#           (in signifier) when making new blast files, add them to archive immediately, don't keep them on disk

my $prok = '/work/rec3141/gene_clustering/prok';
my $trytgz = 1;
my @taxa1;
my @taxa2;
my $forced;

if (scalar(@ARGV) == 0) {die "include at least one blastfile\n";}
elsif (scalar(@ARGV) == 1) {
  @taxa1 = read_file($ARGV[0]) or die ("couldn't read file\n");
  @taxa2 = read_file($ARGV[0]) or die ("couldn't read file\n");
  $forced = 0;
} elsif (scalar(@ARGV) == 2) {
  @taxa1 = read_file($ARGV[0]) or die ("couldn't read file\n");
  @taxa2 = read_file($ARGV[1]) or die ("couldn't read file\n");
  $forced = 1;
} else {die "inputs not correct\n";}
my %hash;

# print "read files\n";
#1 get a list of existing .hits files on the filesystem
#2 get a list of existing .hits files in archives
#3 write instructions for making new .hits files if not found above

my %seen = ();
my @utaxa = grep { ! $seen{$_} ++ } (@taxa1,@taxa2);

for (my $i=0;$i<scalar(@utaxa);$i++) {
 my $taxon1 = $utaxa[$i];
 chomp $taxon1;
 mkdir("$prok/$taxon1") unless (-d "$prok/$taxon1");
 symlink("$prok/$taxon1",$taxon1) unless (defined(readlink($taxon1)));

 opendir(DIR, "$prok/$taxon1");
 my @taxafiles = readdir(DIR);
 close(DIR);

  foreach my $hits (@taxafiles) {
    next unless $hits =~ m/hits/;
    my ($tax1,$tax2) = split(/\./,$hits);
    $hash{$tax1}{$tax2} ='fs';
  }
}

open(my $fhrun,">","newblasts.run") or die "couldn't open newblasts.run";
chmod(0755,$fhrun);
my $newjobs=0;
my $i=0;
my $j=0;

for ($i=0;$i<scalar(@taxa1);$i++) {
  my $taxon1 = $taxa1[$i];
  chomp $taxon1;
  print $taxon1,"\n";
  my $printcount = 0;
  my $fhfile1;
  my $filename;

  my $tar;
  my $tarlink;
  my @tarfiles;
  if ($trytgz == 1) {
    $tarlink = readlink("$prok/$taxon1/$taxon1.hits.tgz");
    if ( -e $tar && defined($tarlink)) {
      $tar = Archive::Tar->new($tarlink,COMPRESS_GZIP);
      @tarfiles = $tar->list_files;  # returns list of Archive::Tar::File obj
      foreach my $file (@tarfiles) {
	my ($t1,$t2) = split(/\./,$file);
	$hash{$t1}{$t2}='tgz' unless $hash{$t1}{$t2};
      }
    }
  }


  TAXON2: for ($j=0;$j<scalar(@taxa2);$j++) {
  if ($printcount == 0) {
    $filename = $taxon1 . "_" . $printcount . ".run";
    open($fhfile1,">>",$filename) or die "couldn't open $fhfile1";
  } 
#   elsif ($printcount % 50 == 0) {
#     print $fhrun "sqsub -q serial -r 7d -o ./$taxon1/$taxon1.jobout ./$filename >> submission_record\n";
#     chmod(0744,$fhfile1);
#     close($fhfile1);
#     $filename = $taxon1 . "_" . $printcount . ".run";
#     open($fhfile1,">>",$filename) or die "couldn't open $fhfile1.run";
#   }

    my $taxon2 = $taxa2[$j];
    chomp $taxon2;

    if ($forced ==1) {
      print $fhfile1 "/home/rec3141/repo/blast_extract.pl $taxon1 $taxon2\n";
      print $fhfile1 "/home/rec3141/repo/blast_extract.pl $taxon2 $taxon1\n";
      $printcount++;
      $printcount++;
      print "n"; $newjobs++;
    } else {
      if ($hash{$taxon1}{$taxon2}) {
	if ($hash{$taxon1}{$taxon2} eq 'fs') {
	  print "e";
	} elsif ($hash{$taxon1}{$taxon2} eq 'tgz') {
	  print $fhfile1 "tar xvzf $prok/$taxon1/$taxon1.hits.tgz $taxon1.$taxon2.hits\n";
	  print "z";
	} else {print "confused...\n\n"}
      } else {
	print $fhfile1 "if [ -L \"$prok/$taxon1/$taxon1.$taxon2.hits\" ]; then echo \"$taxon1.$taxon2.hits.tmp already in progress\"; else ln -s /scratch/rec3141/$taxon1.$taxon2.hits $prok/$taxon1/$taxon1.$taxon2.hits; /home/rec3141/repo/blast_extract.pl $taxon1 $taxon2;fi;\n";
	print "n"; $newjobs++;
	$printcount++;
      }
    }
  }

  print "\n";
  chmod(0755,$fhfile1);
  close($fhfile1);
  print $fhrun "sqsub -q serial -r 7d -o ./$taxon1/$taxon1.jobout ./$filename >> submission_record\n" if ($newjobs>0);
  $newjobs=0;
}

close($fhrun);
