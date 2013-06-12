#!/usr/bin/perl
use strict;
use warnings;
use diagnostics;

use RECGenomeTools qw(:All);

use Bio::Seq;
use Bio::SeqIO;
use Bio::AlignIO;
use Data::Dumper;
#use Env qw(HOME);

sub parse_outputs;

$|++; #flush buffer 

# Summary:

# This script can be called to retrieve gene sequences and other characteristics belonging to any 
# specific cluster(s)

# clusters are currently specified by their index (line number) in the 
# cluster list file, separated by commas

# a cluster number 0 indicates all clusters in file

# to do command-line arguments properly:
# my $debug = 0;
# my $dbname;
# my $output;
# 
# GetOptions(
# 	   'v|verbose!'  => \$opt_debug,
# 	   'o|output:s'  => \$opt_output,
#		'm|maxvals:i' => \$opt_maxvals,
#		'c|cindex:i' => \$opt_cluster_index,
#		'f|cfile:i' => \$opt_cluster_file,
# 	   );
# 
# unless(  defined $dbname ) {
#     die("no dbname provided\n");
# }
# if any {!defined $_} $text, $cols, $filler;


if (scalar(@ARGV) != 6) { 
  die("\nscript takes 6 inputs:\n
    1) a parameter to output for each protein, comma-separated list of 'fastant', 'fastaaa', 'pid', 'cog', 'gene', 'cogfunc', 'projectid', 'taxaname', 'id', 'replicon', 'group', 'subgroup', 'loc', 'desc', or 'evalue', or 'locus2refseq', or 'singlecopy'
       which require the following files
		 prokaryotes.txt	:	projectid, taxaname, id, replicon, taxonomic group, taxonomic subgroup
		 			.hits	:	evalue
					.faa 	:	fastaaa, pid; also required for .ffn and .ptt parameters
					.ffn	:	fastant, loc; also required to find .ptt parameters
		           	 WARNING: 'fastant' sequences from DRAFT genomes with many contigs/scaffolds may be incorrect!
		           	.ptt 	:	cog, gene, cogfunc, desc
    2) maximum number of values to output per cluster; 0 for all
    3) a filename containing clusters of protein ids
    4) a cluster index in the file, 0 for all, comma-separated for multiple
    5) a taxalist file
    6) a taxalist id to limit output clusters to

    eg. ./retrieve_sequences.pl cog 42 NC_000123.cluster_list taxalist-Archaea AR_110001
	    
	script outputs:\n
	1) requested parameters, one per line (if 1 cluster) or summarized (if >1 cluster) on one line
	    \n");
}

### OUTPUT
### a parameter to output for each pid
my $opt_output = shift(@ARGV);
my %inputs= ('opt_output'=>$opt_output);

my @outputs = parse_outputs($opt_output);
my $parse_locus++ if (grep {$_ eq 'locus2refseq'} @outputs);
$inputs{'outputs'} = \@outputs;
$inputs{'parse_locus'} = $parse_locus;

my $singlecopy++ if (grep { $_ =~ s/singlecopy/fastaaa/ } @outputs);


### MAXVALS
### maximum number of values to output per genome (for sequences) or cluster (for others); 0 for all (can be slow)
my $opt_maxvals = shift(@ARGV);
unless ($opt_maxvals =~ m/^\d+$/) { die("$0 failed to understand input argument 2: ", $opt_maxvals, "\n") }
$inputs{'opt_maxvals'}=$opt_maxvals;

### CLUSTER FILE
### a filename containing clusters of protein ids
my $opt_cluster_file = shift(@ARGV);
my $cluster_file;
if ( -e $opt_cluster_file) {
	$cluster_file = $opt_cluster_file;
} else {
	die("$0 failed to understand input argument 3: ", $opt_cluster_file, "\n")
}

$inputs{'cluster_file'}=$cluster_file;

### CLUSTERS
### a cluster index in the file, 0 for all, comma-separated for multiple
### read in cluster index (just an integer or CSV list of integers)
my $opt_cluster_index = shift(@ARGV);

$inputs{'opt_cluster_index'}=$opt_cluster_index;
my $wanted_clusters_ref = parse_cluster_index(\%inputs);
my %wanted_clusters = %{$wanted_clusters_ref}; #requested clusters
$inputs{'wanted_clusters'}=\%wanted_clusters;

my $prot_seen_ref = check_for_clusters(\%inputs);
my %prot_seen = %{$prot_seen_ref}; #proteins in requested clusters
$inputs{'prot_seen'} = \%prot_seen;

### TAXALIST -- a taxalist file
my $taxalist = shift(@ARGV);
$inputs{'taxalist'}=$taxalist;
my $taxahash_ref = parse_taxalist(\%inputs);
my %taxahash = %{$taxahash_ref};
$inputs{'taxahash'}=\%taxahash;

### TAXALIST ID
### a taxalist id to limit clusters to
my $id = shift(@ARGV);
$inputs{'id'}=$id;
my ($ids_ref,$all_ids_ref) = parse_taxaid(\%inputs);
my @ids = @{ $ids_ref }; #list of lists of ids belonging to each genome
my %all_ids = %{ $all_ids_ref }; # $all_ids{'NC_123456'}

$inputs{'ids'}=\@ids;
$inputs{'all_ids'}=\%all_ids;

my ($accgenomes_ref,$genes_ref) = parse_accession_index(\%inputs);
# %accgenomes is a hash to hold id of each protein from requested clusters by genome
# %genes is a hash to hold genome from requested clusters by protein
my %genes = %{ $genes_ref }; # %accgenomes{'YP_006429118.1'} = 'NC_010814'

$inputs{'accgenomes'}=$accgenomes_ref;
$inputs{'genes'}=$genes_ref;

my $ncbi_ref = parse_proks(); #parse prokaryotes.txt into %ncbi
my %ncbi = %{ $ncbi_ref };

$inputs{'ncbi'}=\%ncbi;


# process additional files if necessary

my %pids; # %pids{'1234098'} = 'YP_10901901'
#gets $genes{$refseq}{'faa_obj'} and $genes{$refseq}{'pid'} for all requested genomes
if ( grep {$_ =~ m/(fastaaa|fastant|loc|desc|cog|cogfunc|gene|pid|locus2refseq)/ } @outputs ) {
	my $pids_ref = parse_faa( \%inputs ); 
	%pids = %{ $pids_ref };
}
$inputs{'pids'}=\%pids;

my %locs; # %locs{'NC_012345'}{'123124'}{'234980'} = @[NP_0911223.1 NP_0917541.1 ]
my %loci; # %loci{'DVU_12345'} = 'NP_12345'
#gets ... for all requested genomes
if ( grep {$_ =~ m/(fastaaa|fastant|loc|desc|cogfunc|cog|gene|locus2refseq)/ } @outputs ) {
	my ($locs_ref,$loci_ref) = parse_ptt( \%inputs ); 
	%locs = %{ $locs_ref };
	%loci = %{ $loci_ref };
}
$inputs{'locs'}=\%locs;

#gets ... for all requested genomes
if ( grep {$_ =~ m/(fastant)/} @outputs ) {
	parse_ffn( \%inputs );
}



### populate %out for subsequent output
open(IF1, "< $cluster_file")  || die("couldn't open file $cluster_file: $!\n");
while (<IF1>) {
	next unless $wanted_clusters{$.};
	chomp;
	my $cluster_line = $_;
	my $cluster = $.;
	
	my @cluster_list = split(/\s+/,$cluster_line);
	my @filtered = ();
 	my %taxonids = ();

	if ( $parse_locus ) {
		@filtered = map {lc} @cluster_list;
	} else {
		foreach my $refseq (@cluster_list) {
			next unless exists($genes{$refseq});
			$taxonids{ $genes{$refseq}{'taxonid'} } = undef;
			push(@filtered,$refseq);
		}
	}

	#if genome subset requested, not all clusters will be represented
	unless (scalar(@filtered)) {
		print "\n";
		next;
	}
	
	my $counter = scalar(@outputs);
	foreach my $output (@outputs) {
		$counter--; #decrement
		
		### output information about genomes
		if (grep { $_ eq $output } qw/projectid taxaname id replicon group subgroup/) {
			my @out=();
			foreach my $refseq (@filtered) {
 				push(@out, $ncbi{ $genes{$refseq}{'taxonid'} }{$output});
			}
				print join("\t", grep {defined} @out);
				print "\x01" if $counter; #don't print separator on last $output
				next;
		}

		### output information about genes		
		if (grep { $_ eq $output } qw/loc locs loce strand locus desc cogfunc cog gene pid/) {
			my @out = ();
			foreach my $refseq (@filtered) {
				push(@out, $genes{$refseq}{$output} );
			}
			print join("\t", grep {defined} @out);
			print "\x01" if $counter; #don't print separator on last $output
			next;
		}
		
		### output sequences in fasta format
		if ($output =~ m/(fastant|fastaaa)/) {
			foreach my $refseq (@filtered) {
				my $sequence = Bio::SeqIO->new(-fh => \*STDOUT, '-format' => 'fasta' );
				if (exists($genes{$refseq}{$output})) {
					foreach( sort { $a->length <=> $b->length } @{ $genes{$refseq}{$output} }) {
						$sequence->write_seq($_);
						last if $singlecopy;
					}
				} else {
					print ">$refseq\n";
				}
			}
			next;
		}
		
		### output refseq from loci
		if ($output eq 'locus2refseq') {
			my @out=();
			foreach my $locus (@filtered) {
				push(@out, $loci{$locus});
			}
			print join("\t", grep {defined} @out);
			print "\x01" if $counter; #don't print separator on last $output
		}

		### output evalues
		if ($output eq 'evalue') {
 			get_hits( (\%taxonids,\@filtered, \%genes ) );
 			next;
		}
		
	} #end foreach my $output
	print "\n";
	
} #end while
close(IF1);

exit 0;




#------------------------
# SUB PARSE_OUTPUTS
#------------------------
sub parse_outputs {
	my $opt_output = shift;
	my @output_list = split(',',$opt_output);

	my @option_list = qw|fastant fastaaa cog gene cogfunc taxid projectid taxaname id replicon loc desc evalue group subgroup location pid locus2refseq singlecopy|;
	my @_outputs = ();
	
	foreach my $output (@output_list) {
		if (grep { $_ eq $output } @option_list) {
			push(@_outputs,$output);
		} else {
			die("these outputs are available: ", join(' ',@option_list), "\n");
		}
	}
	
	return(@_outputs);
}





#####################################
# Here is a process to get the nucleotide sequence given the RefSeq ID:

#RefSeq->PID(.faa)->Loc(.ptt)->Seq(.ffn)

#[merz@wha780 program]$ head -5 NC_000964.faa
#>gi|16077069|ref|NP_387882.1| chromosomal replication initiation protein [Bacillus subtilis subsp. subtilis str. 168]
#MENILDLWNQALAQIEKKLSKPSFETWMKSTKAHSLQGDTLTITAPNEFARDWLESRYLHLIADTIYELT
#GEELSIKFVIPQNQDVEDFMPKPQVKKAVKEDTSDFPQNMLNPKYTFDTFVIGSGNRFAHAASLAVAEAP
#AKAYNPLFIYGGVGLGKTHLMHAIGHYVIDHNPSAKVVYLSSEKFTNEFINSIRDNKAVDFRNRYRNVDV
#LLIDDIQFLAGKEQTQEEFFHTFNTLHEESKQIVISSDRPPKEIPTLEDRLRSRFEWGLITDITPPDLET

#[merz@wha780 program]$ head -5 NC_000964.ptt 
#Bacillus subtilis subsp. subtilis str. 168, complete genome - 1..4215606
#4176 proteins
#Location        Strand  Length  PID     Gene    Synonym Code    COG     Product
#410..1750       +       446     16077069        dnaA    BSU00010        -       COG0593L        chromosomal replication initiation protein
#1939..3075      +       378     16077070        dnaN    BSU00020        -       COG0592L        DNA polymerase III subunit beta

#[merz@wha780 program]$ head -5 NC_000964.ffn
#>gi|932480293|ref|NC_000964.3|:410-1750
#ATGGAAAATATATTAGACCTGTGGAACCAAGCCCTTGCTCAAATCGAAAAAAAGTTGAGCAAACCGAGTT
#TTGAGACTTGGATGAAGTCAACCAAAGCCCACTCACTGCAAGGCGATACATTAACAATCACGGCTCCCAA
