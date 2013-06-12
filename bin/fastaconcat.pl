#!/usr/bin/perl
use strict;
use warnings;
use diagnostics;

use RECGenomeTools qw(:All);
use Data::Dumper;
use Bio::SeqIO;

sub get_evalues;
sub average;

die
"\nthis program takes 3 inputs:
1) a taxalist
2) a taxalist id
3) the path to a directory of .fasta alignments to concatenate 
4) a method to select amongst multiple homologs:
	'lowest' (evalue to single-copy genes),
	'average' (of evalues to single-copy genes),
	'longest' (longest protein)
	'lengthdiff' (protein with closest length to single-copy proteins
e.g. fastaconcat.pl taxalist-Bacilli-with-plasmids BL_110001 ./trees average

and outputs:
1) a concatenated fasta file named with RefSeq IDs" unless ( scalar(@ARGV) == 4 );

my %inputs;
my $forceall = 1; #for when I know all the genes are OK

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

my %accgenomes; # %accgenomes{'NC_010814'}{'YP_006429118.1'} = undef
my %genes; # %accgenomes{'YP_006429118.1'} = 'NC_010814'

my @faadirs = qw|/work/rec3141/gene_clustering/all_genomes/faa /work/rec3141/gene_clustering/all_genomes/scaffold/faa /work/rec3141/gene_clustering/all_genomes/contig/faa|;
foreach my $faadir (@faadirs) {
  open(IF,"<$faadir/accession_index") or die "couldn't open file $faadir/accession_index\n";
  while (<IF>) {
	chomp;
	my @acc = split(/\s+/,$_);
	my $genome = shift(@acc);
	next unless exists($all_ids{$genome});
	map { $accgenomes{$genome}{$_} = undef } @acc;
 	map { $genes{$_} = $genome } @acc;
  }
  close(IF);
}

my $dir = shift(@ARGV);
my @files = <$dir/*.fasta>; #SF_110010.404.aa.fasta
my %fasta;
my %lengths;
foreach my $file (@files) {
	my $faaio = Bio::SeqIO->new(-file => $file, -format=>"fasta");
	while (my $faa_obj = $faaio->next_seq) {
		my @display_ids = split(/\|/,$faa_obj->display_id);
		chomp @display_ids;
		my $refseq = $display_ids[3];
		my $genome = $display_ids[5];
		unless (exists($genes{$refseq})) {
			$genes{$refseq} = $genome;
		}
		next unless $display_ids[4] eq 'mol';
		$lengths{$file} = $faa_obj->length;
 		$fasta{$file}{$genome}{$refseq} = $faa_obj;
		
	} #end while sequence
}

my $method = shift(@ARGV); #average, lowest, longest, lengthdiff

my %sequences; # $sequences{$tag}{'seq'} .= $newclusterseq $sequences{$tag}{'refseqs'}
#for each cluster, find best sequences for each tag

	#$file =  SF_110010.405.aa.fasta
	#$id = SF_110010
	#$tag = SF_000001
	#$genomes = NC_012345

my @gen;
my @pro;
foreach my $file (keys %fasta) {
	foreach my $genome (keys %{$fasta{$file}}) {
#		print "$genome\n";
		push @gen,$genome;
		push @pro, keys %{$fasta{$file}{$genome}};
	}
}

my %evalues;
if (grep {$_ eq $method } qw/lowest average/) {
	%evalues = %{ get_evalues( (\@gen,\@pro) ) };
}

print "$id\n";

foreach my $file (sort keys %fasta) {
 	print "$file\n";
	my %choices = ();
	foreach my $tag ( sort grep { defined } @{$taxahash{$id}} ) {
 		print ":$tag\n";
		foreach my $genome ( sort grep { defined } @{$taxahash{$tag}} ) {
 			print "::$genome\n";
			foreach my $refseq (keys %{ $fasta{$file}{$genome} }) {
				if (exists($accgenomes{$genome}{$refseq})) {
  					print ":::$refseq\n";
					push(@{ $choices{$tag} },$refseq);
				} elsif ($forceall) {
  					print ":::forced $refseq\n";
					push(@{ $choices{$tag} },$refseq);
				}
			} #end foreach my $refseq
		} #end foreach my $genome
	} #end foreach my $tag
	
	#at this point we've gone through all the genomes in the cluster
	#if multiple potential matches, for each potential match find evalues to other members of cluster
	#select potential match with best hits to single-member matches
	my @nocopy =     grep { !exists $choices{$_} }                                 @{$taxahash{$id}}; #[SF_000001 SF_000002]
	my @singlecopy = grep { exists $choices{$_} && scalar(@{ $choices{$_} }) == 1 } @{$taxahash{$id}}; #[SF_000001 SF_000002]
	my @multicopy =  grep { exists $choices{$_} && scalar(@{ $choices{$_} })  > 1 } @{$taxahash{$id}}; #[SF_000001 SF_000002]
	
# 	print "$file nocopy: @nocopy\n";
# 	print "$file single: @singlecopy\n";
# 	print "$file multi : @multicopy\n";
# 	print Dumper %choices;


	foreach my $tag ( sort grep { defined } @{$taxahash{$id}} ) {
		if (grep { $_ eq $tag } @singlecopy) {
			my $refseq = @{ $choices{$tag} }[0];
			# $sequences{$tag}{'seq'} .= sprintf('%14.14s',$refseq);
			$sequences{$tag}{'seq'} .= 'N' x 12;
			$sequences{$tag}{'seq'} .= $fasta{$file}{ $genes{$refseq} }{$refseq}->seq;
			push(@{$sequences{$tag}{'refseq'}},$refseq);
		} elsif (grep { $_ eq $tag } @nocopy) {
# 			$sequences{$tag}{'seq'} .= sprintf('%14.14s',' ');
			$sequences{$tag}{'seq'} .= 'N' x 12;
			$sequences{$tag}{'seq'} .=	('-' x $lengths{$file});
			push(@{$sequences{$tag}{'refseq'}},"none");
		} elsif ( grep { $_ eq $tag } @multicopy ) {
			my %best;
			foreach my $prot ( @{ $choices{$tag} } ) {
				my @prots = ($prot, map { @{ $choices{$_} }[0] } @singlecopy);
				my @genomes = map { $genes{$_} } @prots;
				my %evals;
				foreach my $k1 (@prots) {
					foreach my $k2 (@prots) {
						$evals{$k1}{$k2} = $evalues{$k1}{$k2};						
					}
				}

				# if lowest: to pick the protein with the lowest hit to any other single copy protein
				# elsif average: to pick the protein with the lowest average hit to the other single copy proteins
				# elsif longest: to pick the protein with the longest length
				if ($method eq 'lowest') {
# 					%evals = %{ get_evalues( (\@genomes,\@prots) ) };
					foreach my $single (sort { $evals{$prot}{$a} <=> $evals{$prot}{$b} } keys %{ $evals{$prot} } ) {
						$best{$prot} = $evals{$prot}{$single};
						last;
					}
				} elsif ($method eq 'average') {
# 					%evals = %{ get_evalues( (\@genomes,\@prots) ) };
					my @meanarray;
					foreach my $single (sort { $evals{$prot}{$a} <=> $evals{$prot}{$b} } keys %{ $evals{$prot} } ) {
						push(@meanarray,$evals{$prot}{$single});
					}
					$best{$prot} = average(\@meanarray);
				} elsif ($method eq 'longest') {
					$best{$prot} = -1 * ($fasta{$file}{ $genes{$prot} }{$prot}->seq =~ tr/[A-Z]// );
				} elsif ($method eq 'lengthdiff') {
					my $protlength = $fasta{$file}{ $genes{$prot} }{$prot}->seq =~ tr/[A-Z]//;
					my @difflengths = map { $protlength - ($fasta{$file}{ $genes{$_} }{$_}->seq =~ tr/[A-Z]//) } (map { @{ $choices{$_} }[0] } @singlecopy);
					$best{$prot} = average(\@difflengths);
				} else {die "select a method to choose from amongst homologs: lowest, average, longest, lengthdiff\n";}
			}

			foreach my $best_prot ( sort { $best{$a} <=> $best{$b} } keys %best ) {
# 				$sequences{$tag}{'seq'} .= sprintf('%14.14s',$best_prot);
				$sequences{$tag}{'seq'} .= 'N' x 12;
				$sequences{$tag}{'seq'} .= $fasta{$file}{ $genes{$best_prot} }{$best_prot}->seq;
				push(@{$sequences{$tag}{'refseq'}},$best_prot);

				last;
			}
			
		} else {
 			# print "$tag not found\n";
		}
	} #end foreach my $tag
	
} #end foreach my $file



#>TAG1
#seq1seq2seq3seq4seq5
#>TAG2
#seq1seq2seq3seq4seq5

mkdir("trees") unless ( -d './trees' );
open(OUT,">","./trees/$id.concat.$method.faa");
foreach my $tag (sort keys %sequences) {
	print OUT ">$tag ",join(' ',@{$sequences{$tag}{'refseq'}}),"\n";
	print OUT $sequences{$tag}{'seq'},"\n";
}
close(OUT);


open(OUT,">","./trees/$id.$method.raxml.q.txt");
my $length=1;
foreach my $file (sort keys %lengths) {
	my $newlength = $length + $lengths{$file} + 12;
	my (undef,$clusternum) = split(/\./,$file);	
	print OUT "WAG, clust$clusternum = ", $length, "-", $newlength,"\n";
	$length = $newlength;
}
close(OUT);

exit 0;


sub get_evalues {
# honestly i'm not sure what this is supposed to be useful for....
# maybe for subsetting a cluster into MCL??
# ps not sure it works properly

# 	my $taxonids_ref = shift;
# 	my @ids = keys %{$taxonids_ref};
# 	my $cluster_list_ref = shift;
# 	my @_proteins = @{$cluster_list_ref};
# 	my $genes_ref = shift;
# 	my %genes = %{ $genes_ref };
	my $cutoff = 0.6;
	my @ids = @{ shift(@_) };
	my @_proteins = @{ shift(@_) };
	my $prokdir = '/work/rec3141/gene_clustering/prok';
	my %eval_out;
	foreach my $subject (sort @ids) {
		foreach my $query (sort @ids) {
			my %evalues = ();
			open(IF,"<$prokdir/$subject/$subject.$query.hits") or die("couldn't open hits file $prokdir/$subject/$subject.$query.hits\n");
			while(<IF>) {
				chomp;
				my @line = split(/\s+/,$_);
				$evalues{$line[0]}{$line[1]} = $line[2] if ( $line[3]>$cutoff || $line[4]>$cutoff);
			}
			close(IF);
		
			my @w = keys(%evalues);
			foreach my $e (keys %evalues) {
				push (@w,keys %{$evalues{$e}});
			}
			
			my %union = ();
			my %isect = ();
			foreach my $f (@w) { $union{$f} = 1 }
			foreach my $g (@_proteins) {
				if ( $union{$g} ) { $isect{$g} = 1 }
				$union{$g} = 1;
			}		

			foreach my $i (keys %isect) {
				foreach my $j (keys %isect) {
					if (exists($evalues{$i}{$j})) {
						$eval_out{$i}{$j} = $evalues{$i}{$j};
					}
				}
			}
		} #end foreach my query
	} #end foreach my subject
	
	return \%eval_out;
} #end sub


sub average{
        my($data) = @_;
        if (not @$data) {
                die("Empty array\n");
        }
        my $total = 0;
        foreach (@$data) {
                $total += $_;
        }
        my $average = $total / @$data;
        return $average;
}
