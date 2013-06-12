#!/usr/bin/perl
use Bio::DB::Taxonomy;
use Data::Dumper;
use strict;

use warnings;
# use diagnostics;

die
"\nthis program takes 2 inputs:
1) a taxon name  -- e.g. 'Deltaproteobacteria' will find all Deltaproteobacteria
  OR a filename containing a list of Taxon Ids or Taxon Names or a mixture of both
  (taxa to exclude may be prefixed by a '-')
2) a prefix to use for new identifiers

e.g. gettaxa.pl NCBI-summary.txt Deltaproteobacteria DP


and outputs:
1) a blastlist   - RefSeq ids for blasting
2) a taxalist    - indices for groups for clustering
3) a key list    - key for phylogenetics 
4) a taxid list  - taxon id, project id, and taxon name suitable for input into this program
\n\n" if ( scalar(@ARGV) != 2 );


### SET UP TAXONOMIC DATABASE

my $idx_dir = '/work/rec3141/gene_clustering/all_genomes/idx/';

my ( $nodefile, $namesfile ) = ( 'nodes.dmp', 'names.dmp' );

my $db = new Bio::DB::Taxonomy( -source    => 'flatfile',
                                -nodesfile => $nodefile,
                                -namesfile => $namesfile,
                                -directory => $idx_dir
);

die "Couldn't load taxonomic database\n\n" unless $db;



### FIND WHICH TAXONOMIC GROUPS WERE REQUESTED

my $go_up_tree = 0; #taxonomic levels to travel up tree
my $printname = $ARGV[0]; $printname =~ s#[ /']#-#g;
my $prefix   = $ARGV[1];
my @requestedids; #taxon ids
my @excludeids; #taxon ids, to skip

# IS IT A FILE?
if ( open( IF, "<$ARGV[0]" ) ) {
    while (<IF>) {
		chomp;
    	my @line = split(/\t/,$_);
    	# list can be taxonids or taxa names or mixture
    	if ($line[0] =~ m/^(-)?(\d+)$/) {
	  $1 eq '-' ? push(@excludeids,$2) : push(@requestedids,$line[0]);
	} elsif ($line[0] =~ m/^(-)?(\w+.*)$/) {
	  $1 eq '-' ? push(@excludeids,$db->get_taxonids($2)) : push(@requestedids,$db->get_taxonids($line[0]));
	}
    }
    close(IF);
} else {
    push( @requestedids, $db->get_taxonids($ARGV[0]) );
}

die "couldn't find requested taxonid\n" unless (@requestedids);


### FIND WHICH GENOMES BELONG TO REQUESTED GROUPS

my $draftfile = "/work/rec3141/gene_clustering/all_genomes/prokaryotes.txt";

my %taxahash;
my %hasgenome;

## MAKE AN INDEX OF KNOWN MICROBIAL GENOMES
## indexed by BioProject ID

open (IF, "<$draftfile") or die "couldn't open $draftfile file\n";
while (<IF>) {
	next if $.==1; #skip header line

	# 11/2/2012 LPROKS OBSOLETE, USING NEW PROKARYOTES.TXT
	#Organism/Name	BioProject Accession	BioProject ID	Group	SubGroup	Size (Mb)	GC%	Chromosomes/RefSeq	Chromosomes/INSDC	Plasmids/RefSeq	Plasmids/INSDC	WGS	Scaffolds	Genes	Proteins	Release Date	Modify Date	Status	Center
	# Colwellia psychrerythraea 34H	PRJNA57855	57855	Proteobacteria	Gammaproteobacteria	5.37318	38	NC_003910.7	CP000083.1	-	-	-	1	5054	4910	2002/05/16	2012/01/20	Complete	TIGR
	# Abiotrophia defectiva ATCC 49176	PRJNA55729	55729	Firmicutes	Bacilli	3.4774	37.1	-	-	-	-	ACIN02	8	3346	3291	2009/03/17	2010/02/03	Scaffolds or contigs	Washington University Genome Sequencing Center

	my ($taxaname, undef , $projectid, undef ,undef,undef,undef,$chromosomes,undef,$plasmids,undef,$wgs) = split( "\t", $_ );

	## get Taxonomic ID
	## used to be from summary.txt and lproks_2.txt but no more
	## now need to use Bio::DB::Taxonomy
	
    my @taxonids = $db->get_taxonids($taxaname);    #just taxonids
    if ( @taxonids < 1 ) {
    	# these species are not found in the current taxonomic database
    	# haven't looked into why
		print "skipping. NOT FOUND IN TAXONOMY DB: $taxaname\n";
		next;
	} elsif ( @taxonids >1 ) { 
		# get_taxonids does case-insensitive matching so....
		# 	Bacillus sp. M2-9 <-- has genome record
		# 	Bacillus sp. m2-9 <-- does not have genome record
		# but they are both found
		# only 1 known case... for now, skipping with error message
		print "skipping. MULTIPLE MATCHES IN TAXONOMY DB: $taxaname\n"if scalar(@taxonids)>1;
		next;
	}

	my $taxid = shift(@taxonids);


	## hash information about sequence accessions
	# some strains have multiple genome sequences or have them in progress
	# for now we don't care, one is enough
	# further below we'll choose the complete sequence if it is available
	# else choose the draft
	
    if ($wgs !~ m/-/) {
    	$wgs =~ s/^(\w{4})\d+/$1/;
		$taxahash{$taxid}{'scaffolds'} = "NZ_".$wgs."00000000";
    } elsif ($chromosomes !~ m/-/) {
		$taxahash{$taxid}{'chromosomes'} = $chromosomes;
		$taxahash{$taxid}{'plasmids'} = $plasmids unless ($plasmids =~ m/-/);
	} elsif ($plasmids !~ m/-/) {
		$taxahash{$taxid}{'plasmids'} = $plasmids;
    } else {
    	#no data, skip
    	next;
    }

	$taxahash{$taxid}{'taxaname'} = $taxaname;
	$taxahash{$taxid}{'taxid'} = $taxid;
	$taxahash{$taxid}{'projectid'} = $projectid;
	
	$hasgenome{$taxid}++; #set up index for taxa with genomes

}
close IF;


### Get all descendants of the requested taxa
my @ancestors; #Bio::Taxon objects
@ancestors = map { $db->get_taxon(-taxonid => $_) } @requestedids;

## undocumented: Get ancestors of the requested taxa
while ($go_up_tree > 0) {
  my @newancestors = map { $db->ancestor($_) } @ancestors;
  for(my $i=0;$i<scalar(@ancestors);$i++) {
    print $ancestors[$i]->rank,"\t",$ancestors[$i]->scientific_name, "\t-->\t", $newancestors[$i]->rank,"\t",$newancestors[$i]->scientific_name,"\n";
  }
  @ancestors = @newancestors;
  $go_up_tree--;
}

my @descendants = @ancestors; #Bio::Taxon objects
foreach my $taxon (@ancestors) {
    push(@descendants,$db->get_all_Descendents($taxon));
}

my %req_genomes; #hash keys->taxonids, values->counters
{
  no warnings;
  #keep only descendants with genomes
  @descendants = grep { $hasgenome{$_->id} } @descendants; 
  # print "\n $printname with genomes:\n";
  # print join("\n",map {$_->scientific_name} @descendants);
  map { $req_genomes{$_->id}++ } @descendants;
  undef @descendants;
}

#process excluded genomes
my @excluded; #Bio::Taxon objects
@excluded = map { $db->get_taxon(-taxonid => $_) } @excludeids;
foreach my $taxon (@excluded) {
    push(@excluded,$db->get_all_Descendents($taxon));
}
map { print $_->scientific_name," excluded\n" } @excluded;
map { delete($req_genomes{$_->id}) } @excluded;

### OUTPUT NECESSARY FILES

my $newaccession = 1;
my $pad_len      = 6;

my @allmols;

open( OUTTAXA,  ">taxalist-$printname" );
open( OUTBLAST, ">blastlist-$printname" );
open( OUTKEY, ">key-$printname");
open( OUTIDS, ">ids-$printname");

foreach my $k1 (sort { $taxahash{$a}{'taxaname'} cmp $taxahash{$b}{'taxaname'} } keys %taxahash) {

	#skip genome unless it has been requested
	next unless $req_genomes{$k1};

	if ($taxahash{$k1}{'chromosomes'}) {
	  my $chr = $taxahash{$k1}{'chromosomes'};
	  $chr =~ s/\.\d//g;
	  my @chr = split(',',$chr);
	  my @sortchr = sort(@chr);
	  my $pls;
	  if ($taxahash{$k1}{'plasmids'}) {
	    $pls = $taxahash{$k1}{'plasmids'};
	    $pls =~ s/\.\d//g;
	  } else {
	    $pls = '';
	  }
	  my @pls = split(',',$pls);
	
# 		print "COMPLETE\t",$k1,"\t",$taxahash{$k1}{'taxaname'},$taxahash{$k1}{'chromosomes'},"\t","\n";
	
	  my $padded = sprintf( "$prefix\_%0${pad_len}d", $newaccession );
	
	  push(@allmols,$padded);

	  print OUTTAXA $padded," ",join(' ',(@chr,@pls)),"\n";
	  print OUTBLAST join("\n",(@chr,@pls)),"\n";
	  print OUTKEY $sortchr[0], "\t",$taxahash{$k1}{'taxaname'},"\n";
	  print OUTIDS $taxahash{$k1}{'taxid'},"\t",$taxahash{$k1}{'projectid'},"\t",$taxahash{$k1}{'taxaname'},"\n";

	} elsif ($taxahash{$k1}{'scaffolds'}) {
# 		print "SCAFFOLD\t",$k1,"\t",$taxahash{$k1}{'taxaname'},$taxahash{$k1}{'scaffolds'},"\t","\n";

	  my $padded = sprintf( "$prefix\_%0${pad_len}d", $newaccession );
	  push(@allmols,$padded);

	  print OUTTAXA $padded," ",$taxahash{$k1}{'scaffolds'},"\n";
	  print OUTBLAST $taxahash{$k1}{'scaffolds'},"\n";
	  print OUTKEY $taxahash{$k1}{'scaffolds'}, "\t",$taxahash{$k1}{'taxaname'},"\n";
	  print OUTIDS $taxahash{$k1}{'taxid'},"\t",$taxahash{$k1}{'projectid'},"\t",$taxahash{$k1}{'taxaname'},"\n";

	} else {
		# this could be builds that have only plasmids?
	  print "error... no valid molecules for $taxahash{$k1}{'taxaname'}?\n";
	  next;
	}
	$newaccession++;
}

print OUTTAXA "$prefix\_110000 ", join(' ',@allmols), "\n";

close(OUTTAXA);
close(OUTBLAST);
close(OUTKEY);
close(OUTIDS);

print "DONE\n";
exit 0;
