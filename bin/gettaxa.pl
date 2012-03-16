#!/usr/bin/perl -w

use lib '/home/rec3141/parser/BioPerl-1.6.0';

use Bio::Perl;
use Bio::DB::Taxonomy;
use Data::Dumper;
use strict;

use warnings;
use diagnostics;

# 1) an NCBI genome summary file
die
"\nthis program takes 2 inputs:
1) a taxon name -- use 'cellular organisms' with quotes for all
  OR a filename containing a tab-delimited list of Taxon Names and Taxon Ids
2) a prefix to use for new identifiers

e.g. gettaxa.pl NCBI-summary.txt Archaea AR

and outputs:
1) a blastlist- file containing RefSeq ids for blasting, with and without plasmids
2) a taxalist- file containing indices for groups, with and without plasmids\n\n"

if ( scalar(@ARGV) != 2 );

my $idx_dir = '/work/rec3141/gene_clustering/all_genomes/idx/';

my ( $nodefile, $namesfile ) = ( 'nodes.dmp', 'names.dmp' );

#my $db = Bio::DB::Taxonomy->new(-source => 'entrez');
my $db = new Bio::DB::Taxonomy( -source    => 'flatfile',
                                -nodesfile => $nodefile,
                                -namesfile => $namesfile,
                                -directory => $idx_dir
);

my $default = 'cellular organisms';
my @requestednames;
push( @requestednames, $ARGV[0] );
my @requestedids;
my $printname = $ARGV[0];
$printname =~ s#[ /']#-#g;
my $prefix   = $ARGV[1];
my $ncbifile = "/work/rec3141/gene_clustering/all_genomes/NCBI-genome-summary.newest"; #$ARGV[0];
my $draftfile = "/work/rec3141/gene_clustering/all_genomes/lproks_2.txt";

my %taxahash;
my %multiples;
my @taxaorder;

open (IF, "<$draftfile") or die "couldn't open $draftfile file\n";
while (<IF>) {
## Columns:     "RefSeq project ID"     "Project ID"    "Taxonomy ID"   "Organism Name" "Super Kingdom" "Group" "Sequence availability" "RefSeq Accession"     "Genbank Accession"     "Number of Contigs"     "Number of proteins on contigs" "Genome Size"   "GC Content"    "Released date""Center Name"   "Ceneter URL"
# 55729   33011   592010  Abiotrophia defectiva ATCC 49176        Bacteria        Firmicutes      wgs assembly    NZ_ACIN00000000 ACIN00000000  26       3291    3.4774  37.1    03/17/09        Genome Sequencing Center (GSC) at Washington University (WashU) School of Medicine      http://www.genome.wustl.edu/sub_genome_group.cgi?GROUP=3&SUB_GROUP=4
my ($refseq,$projectid,$taxid,$taxaname,$kingdom,$group,$avail,$accession) = split( "\t", $_ );
    next if $refseq =~ m/\#\#/;
    $taxahash{$accession}{'taxaname'} = $taxaname;
    $taxahash{$accession}{'replicon'} = 'scaffold';
    $taxahash{$accession}{'taxid'}    = $taxid;
    $taxahash{$accession}{'projectid'} = $projectid;
    push( @taxaorder, $accession );
}
close IF;

open( IF, "<$ncbifile" ) or die "couldn't open $ncbifile file\n";
while (<IF>) {
# Accession       GenbankAcc      Length  Taxid   ProjectID       TaxName Replicon        Create Date     Update Date
# NC_000907.1     L42023.1        1830138 71421   57771   Haemophilus influenzae Rd KW20  chromosome      Oct 19 2001     Mar 17 2011  7:06:34:310PM

    my ( $accession, $genbankacc, $length, $taxid, $projectid, $taxaname, $replicon ) = split( "\t", $_ );
    #accessions are unique per molecule (RefSeq numbers)
    #taxids are shared by molecules of the same strain in different projects
    #projectids are shared by all molecules of the same project
    next if $accession =~ m/Accession/;
    $taxahash{$accession}{'taxaname'} = $taxaname;
    $taxahash{$accession}{'replicon'} = $replicon;
    $taxahash{$accession}{'taxid'}    = $taxid;
    $taxahash{$accession}{'projectid'} = $projectid;
    push( @taxaorder, $accession );
}
close IF;

my %seentaxa;
my %usedtaxa;
my %roottaxa;
my %rootplasmids;
my %rootscaffolds;

foreach my $k1 ( sort keys %taxahash ) {
    my $node = $db->get_taxon( -taxonid => $taxahash{$k1}{'taxid'} );    #a taxon object
    unless ( defined( $node ) ) {
	print $taxahash{$k1}{'replicon'} ." of $k1 not found in taxahash\n";
	next;
}

    #follow up the tree to get ancestors
    push( @{ $taxahash{$k1}{'ancestors'} }, $node );
    my $ancestor;
    while ( defined( @{ $taxahash{$k1}{'ancestors'} }[0]->ancestor ) ) {
        $ancestor = @{ $taxahash{$k1}{'ancestors'} }[0]->ancestor;
        unshift( @{ $taxahash{$k1}{'ancestors'} }, $ancestor );
        unless ( defined( @{ $usedtaxa{ @{ $taxahash{$k1}{'ancestors'} }[0]->scientific_name } } ) ) {
            push( @{ $usedtaxa{ @{ $taxahash{$k1}{'ancestors'} }[0]->scientific_name } },
                  ( @{ $taxahash{$k1}{'ancestors'} }[0]->id, @{ $taxahash{$k1}{'ancestors'} }[0]->rank )
            );
        }
    }
}

my $baseid = @{ $usedtaxa{$default} }[0];

if ( open( IF, "<$requestednames[0]" ) ) {
    shift @requestednames;
    while (<IF>) {
        chomp;
	my @extras = split('\t',$_);
	my @extrataxid;
	if ($extras[0] =~ m/^\d+/) {
		$extrataxid[0] = $extras[0];
	} else {
	@extrataxid = split('\s',$extras[1]);
        push( @requestednames, $extras[0] );
	}
        push( @requestedids, @extrataxid ) if @extrataxid;
    }
    close(IF);
} else {
    push( @requestedids, @{ $usedtaxa{ $requestednames[0] } }[0] );
}
# print join("\n",@requestedids) , "\n";

die "couldn't find requested taxon: $requestednames[0]\n" unless defined( $requestedids[0] );

my @reqarray;

foreach my $k1 ( sort keys %taxahash ) {
    my $node = $db->get_taxon( -taxonid => $taxahash{$k1}{'taxid'} );    #a taxon object
    unless ( defined( $node ) ) { 
# 	print "$k1 not found\t" . $taxahash{$k1}{'replicon'} . "\n";
	next;
}

    if ( grep { $_->id == $baseid } @{ $taxahash{$k1}{'ancestors'} } ) {
        my $brieftaxa = $k1;
        $brieftaxa =~ s/(\D{2}\_\d+)\.\d+/$1/;
#         push( @{ $roottaxa{ $taxahash{$k1}{'taxaname'} } }, $brieftaxa )
        push( @{ $roottaxa{ $taxahash{$k1}{'projectid'} } }, $brieftaxa )
          if ( $taxahash{$k1}{'replicon'} =~ m/chromosome/ );
          #multiple sequences of the same organism (species? strain?) can show up here
#         push( @{ $rootplasmids{ $taxahash{$k1}{'taxaname'} } }, $brieftaxa )
        push( @{ $rootplasmids{ $taxahash{$k1}{'projectid'} } }, $brieftaxa )
        if ( $taxahash{$k1}{'replicon'} =~ m/plasmid/ );
#?	if ( $taxahash{$k1}{'replicon'} =~ m/circular/ );
	if ( $taxahash{$k1}{'replicon'} =~ m/scaffold/ ) {
	  push( @{ $rootscaffolds{ $taxahash{$k1}{'projectid'} } }, $brieftaxa );
          push( @{ $roottaxa{ $taxahash{$k1}{'projectid'} } }, $brieftaxa );
	}
    }
    foreach my $requestedid (@requestedids) {
        if ( grep { $_->id eq $requestedid } @{ $taxahash{$k1}{'ancestors'} } ) {
#             push( @reqarray, $taxahash{$k1}{'taxaname'} );
            push( @reqarray, $taxahash{$k1}{'projectid'} );
        }

    }
}

#print(Dumper(%rootscaffolds));

open( OUT, ">allseen" );
foreach my $usedk ( sort keys %usedtaxa ) { print OUT $usedk, "\t", join( ' ', @{ $usedtaxa{$usedk} } ), "\n"; }
close(OUT);

open( OUT,  ">taxalist-$printname" );
open( OUT3, ">taxalist-$printname\-with-plasmids");
open( OUT9, ">taxalist-$printname\-with-scaffolds" );
open( OUT2, ">blastlist-$printname" );
open( OUT5, ">blastlist-$printname\-plasmids" );
open( OUT10, ">blastlist-$printname\-scaffolds" );
open( OUT6, ">key-$printname");
open( OUT7, ">key-$printname\-plasmids" );
open( OUT8, ">key-$printname\-scaffolds" );
my $newaccession = 1;

my $newplasmidacc = 1;
my $withplasmid = 220001;

my $newscaffoldacc =1;
my $withscaffold = 330001;

my $pad_len      = 6;

my @justplasmids;
my @withplasmids;
my @noplasmids;
my @justscaffolds;
my @withscaffolds;

foreach my $key1 (@taxaorder) {
#     next unless defined( $taxahash{$key1}{'taxaname'} );
#     my $k1 = $taxahash{$key1}{'taxaname'};
    next unless defined( $taxahash{$key1}{'projectid'} );
    my $k1 = $taxahash{$key1}{'projectid'};
    
    next unless defined( $roottaxa{$k1} );
    if ( scalar( @{ $roottaxa{$k1} } ) > 1 ) {
        if ( grep { $_ eq $k1 } @reqarray ) {
            my $padded = sprintf( "ZZ_%0${pad_len}d", $newaccession );
            my $paddedplasmids = sprintf( "$prefix\_%0${pad_len}d", $newplasmidacc);
	    my $paddedscaffolds = sprintf( "$prefix\_%0${pad_len}d", $newscaffoldacc);
            print OUT $padded, " ", join( ' ', @{ $roottaxa{$k1} } ), "\n";
            print OUT2 join( "\n", @{ $roottaxa{$k1} } ), "\n";
            print OUT5 join( "\n", @{ $rootplasmids{$k1} } ), "\n" if defined( @{ $rootplasmids{$k1} } );
            print OUT10 join( "\n", @{ $rootscaffolds{$k1} } ), "\n" if defined( @{ $rootscaffolds{$k1} } );
	    print OUT3 $paddedplasmids, " ", join( ' ', (@{ $roottaxa{$k1} } , @{ $rootplasmids{$k1} } )), "\n";
	    print OUT9 $paddedscaffolds, " ", join( ' ', (@{ $roottaxa{$k1} } , @{ $rootscaffolds{$k1} } )), "\n";
#             push( @noplasmids, @{ $roottaxa{$k1} } ); #to include these ids and not the new ids in the blastlist-* file
            print OUT6 join( "\n", @{ $roottaxa{$k1} } ), "\t",$taxahash{$key1}{'taxaname'},"\n";
	    print OUT7 join( "\n", @{ $rootplasmids{$k1} } ), "\t", $taxahash{$key1}{'taxaname'},"\t",$taxahash{$key1}{'replicon'},"\n" if defined( @{ $rootplasmids{$k1} } );
	    print OUT8 join( "\n", @{ $rootscaffolds{$k1} } ), "\t", $taxahash{$key1}{'taxaname'},"\t",$taxahash{$key1}{'replicon'},"\n" if defined( @{ $rootscaffolds{$k1} } );

            push( @noplasmids, $padded );
            push( @justplasmids, @{ $rootplasmids{$k1} } ) if defined( $rootplasmids{$k1} );
            push( @withplasmids, $paddedplasmids);
	    push( @justscaffolds, @{ $rootscaffolds{$k1} } ) if defined( $rootscaffolds{$k1} );
            push( @withscaffolds, $paddedscaffolds);

            $newplasmidacc++;
	    $newscaffoldacc++;
        }
        $newaccession++;
    } else {
        if ( grep { $_ eq $k1 } @reqarray ) {
            my $padded = sprintf( "$prefix\_%0${pad_len}d", $withplasmid );
            my $paddedplasmids = sprintf( "$prefix\_%0${pad_len}d", $newplasmidacc);
            my $paddedscaffolds = sprintf( "$prefix\_%0${pad_len}d", $newscaffoldacc);
            print OUT join( ' ', @{ $roottaxa{$k1} } ), "\n";
            print OUT2 join( "\n", @{ $roottaxa{$k1} } ), "\n";
            print OUT5 join( "\n", @{ $rootplasmids{$k1} } ), "\n" if defined( @{ $rootplasmids{$k1} } );
            print OUT10 join( "\n", @{ $rootscaffolds{$k1} } ), "\n" if defined( @{ $rootscaffolds{$k1} } );
	    print OUT3 $paddedplasmids, " ", join( ' ', ( @{ $roottaxa{$k1} } , @{ $rootplasmids{$k1} } )), "\n";
	    print OUT9 $paddedscaffolds, " ", join( ' ', ( @{ $roottaxa{$k1} } , @{ $rootscaffolds{$k1} } )), "\n";
            #		print OUT3 " ", join(" ",@{$roottaxa{$k1}});
            print OUT6 join( "\n", @{ $roottaxa{$k1} } ), "\t",$taxahash{$key1}{'taxaname'},"\n";
	    print OUT7 join( "\n", @{ $rootplasmids{$k1} } ), "\t", $taxahash{$key1}{'taxaname'},"\t",$taxahash{$key1}{'replicon'},"\n" if defined( @{ $rootplasmids{$k1} } );
	    print OUT8 join( "\n", @{ $rootscaffolds{$k1} } ), "\t", $taxahash{$key1}{'taxaname'},"\t",$taxahash{$key1}{'replicon'},"\n" if defined( @{ $rootscaffolds{$k1} } );

            push( @noplasmids, @{ $roottaxa{$k1} } );
            push( @justplasmids, @{ $rootplasmids{$k1} } ) if defined( $rootplasmids{$k1} );
            push( @justscaffolds, @{ $rootscaffolds{$k1} } ) if defined( $rootscaffolds{$k1} );
            push( @withscaffolds, $paddedscaffolds);


        $withplasmid++;
        $newplasmidacc++;
	$newscaffoldacc++;
	$withscaffold++;
        }

    }
    undef( $roottaxa{$k1} );
}
print OUT "$prefix\_110000 ", join( ' ', ( @noplasmids, @justplasmids, @justscaffolds ) ), "\n";
print OUT "$prefix\_100000 ", join( ' ', @noplasmids ),   "\n";
print OUT "$prefix\_010000 ", join( ' ', @justplasmids ), "\n";
print OUT "$prefix\_030000 ", join( ' ', @justscaffolds ), "\n";
print OUT3 "$prefix\_110000 ", join( ' ', ( @noplasmids, @justplasmids, @justscaffolds ) ), "\n";
print OUT3 "$prefix\_110001 ", join( ' ', @withplasmids ), "\n";
print OUT3 "$prefix\_130001 ", join( ' ', @withscaffolds ), "\n";

close(OUT);
close(OUT2);
close(OUT3);
close(OUT5);
close(OUT6);
close(OUT7);
close(OUT8);
close(OUT9);
close(OUT10);
