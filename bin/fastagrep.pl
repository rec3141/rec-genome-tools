#!/usr/bin/perl

# FASTA grep tool by Tim Booth - tbooth@ceh.ac.uk, http://nebc.nerc.ac.uk
#
# This script may be redistributed under the same terms as Perl itself.

use strict; use warnings;
use Getopt::Std;
use IO::File;
use Data::Dumper;

#Not using BioPerl at all - reading FASTA is too easy

our $DEFAULT_WIDTH = 80;
our $BLURB = " FASTAgrep - extract sequences from a multi-FASTA file by regex.

  Usage: fastagrep [opts] <pattern> infile > outfile

  Options: grep-style options
        -F = fixed mode, no regex
        -x = whole line match
        -X = match whole of first word (ie. the ID)
        -v = invert match
        -i = ignore case 
        -f FILE = take patterns from file, one per line

      FASTA options
        -w COLS = re-wrap to new width, 0 for no wrap.  Will also remove blank lines
        -s = search in the actual sequence, not the header, implies -w$DEFAULT_WIDTH
		-C = force use of DOS/Windows CRLF line terminators
";
our $VERSION = 0.3;
sub HELP_MESSAGE { print $BLURB };

our %opts = ();
$Getopt::Std::STANDARD_HELP_VERSION = 1;
getopts('FxXvif:w:sC', \%opts) || HELP_MESSAGE();

our @patterns;
our %patternhash;
our $rewrap;
our $invert = ($opts{v} xor 0);

#What to cut off the info line to get the ID.
our $id_chopper = qr/[\s|,:;\/].*/;

#Options -X and -s are meaningless but will work
#Options -i and -F are now supported together, with speedup for FX and large ID list
($opts{x} && $opts{X}) and die "Bad options - X and x flags are mutually exclusive.\n";
$opts{C} and $/ = "\015\012";

#Things to do if naughty carriage returns are found
sub checkcr1
{
	$_[0] =~ /\015/ and die ( $/ eq "\015\012" ? 
		"The first search pattern contains embedded CR line terminators despite CRLF line temrinators being set.\n" . 
		"This input file really needs cleaning up.\n" :
		"Your pattern file contains CR line terminators.  You should run 'sed -i \"s/\\r//\" $opts{f}' to purge them.\n" );
	#Note that if a Unix file is read in DOS mode it ends up all on one line.
 	$_[0] =~ /\012/ and die ( 
 		"The first search pattern contains embedded LF line terminators, so something went wrong when " .
 		"reading the file.\nMaybe you are trying to read a Unix-type file in Windows or with -C in effect?\n" );
					
}
sub checkcr2
{
	$_[0] =~ /\015/ and die ( $/ eq "\015\012" ? 
		"The first line of your FASTA file contains embedded CR line terminators despite CRLF line terminators being set.\n" . 
		"This input file really needs cleaning up.\n" :
		"Your FASTA file contains CR line terminators.  You should run 'sed -i \"s/\\r//\" <filename>' to purge them,\n" .
        "or use the -C flag if you are sure you want to keep all files in this format.\n");
}

#Option f means we have to go and find the patterns in a file.
if($opts{f})
{
	my $patternfile = new IO::File("<$opts{f}") or die "Unable to open pattern file $opts{f}\n";
	if($opts{F})
	{
		if($opts{i})
		{
			@patterns = map {chomp; lc($_)} <$patternfile>;
		}
		else
		{
			@patterns = map {chomp; $_} <$patternfile>;
		}
	}
	elsif($opts{x} || $opts{X})
	{
		#For patterns deal with -x up-front, for fixed match deal with it later.
		if($opts{i})
		{
			@patterns = map {chomp;	qr{^$_$}i} <$patternfile>;
		}
		else
		{
			@patterns = map {chomp; qr{^$_$}} <$patternfile>;
		}
	}
	else
	{
		if($opts{i})
		{
			#Have to replace empty patterns due to Perl silliness.
			@patterns = map {chomp; ($_ eq '' ? qr{|} : qr{$_}i)} <$patternfile>;
		}
		else
		{
			@patterns = map {chomp; ($_ eq '' ? qr{|} : qr{$_})} <$patternfile>;
		}
	}
	#Finally, check to see if any unwanted carriage returns crept in
	checkcr1($patterns[0]) if @patterns;
}
else
{
	#Similar to above
	my $patstring = shift;
	defined $patstring or die $BLURB;
	if($opts{F})
	{
		if($opts{i})
		{
			@patterns = split("\n", lc($patstring));
		}
		else
		{
			@patterns = split("\n", $patstring);
		}

		@patterns = ("") unless @patterns;
	}
	elsif($opts{x} || $opts{X})
	{
		if($opts{i})
		{
			@patterns = map {qr{^$_$}i} split(/\|/, $patstring);
		}
		else
		{
			@patterns = map {qr{^$_$}} split(/\|/, $patstring);
		}
		@patterns = (qr{^$}) unless @patterns;
	}
	else
	{
		if($patstring eq '')
		{
			@patterns = qr{|};
		}
		elsif($opts{i})
		{
			@patterns = qr{$patstring}i;
		}
		else
		{
			@patterns = qr{$patstring};
		}
	}
}

#If there is nothing to match then we are done.
#Nope - still want to reformat
#(!@patterns && !$opts{v}) and exit 0;
#(@patterns[0] eq '' && $opts{v}) and exit 0;

#die Dumper \@patterns;

#Now if I have more than about 5 patterns, and mode is Fx or FX, I can make things
#faster (or should that be fasta?) with a hash.
if(@patterns >= 5  && $opts{F} && ($opts{x} || $opts{X}))
{
	$patternhash{pop(@patterns)}++ while @patterns;
}

# Re-formatting?  Internally 0 = no reformat, -1 = no wrap, N = N cols
$rewrap = 0;
if(defined($opts{w})) { $rewrap = ($opts{w} > 0 ? $opts{w} : -1) }
elsif($opts{s})       { $rewrap = $DEFAULT_WIDTH; }

#Decide how I want to match.  Sorry this got a bit messy, but I wanted to eliminate
#all possible conditionals from the matcher routine.
my $matcher;
if($opts{F}) { if($opts{x}) { if(%patternhash) { if($opts{i}) {
    #Fixed case insensitive match to whole info line based on hash
	$matcher = sub{ return $patternhash{lc($_[0])} }
} else {
	#Fixed match to whole info line based on hash
    $matcher = sub{ return $patternhash{$_[0]} }  
} } else { if($opts{i}) {
    #Fixed case insensitive match to whole info line based on array
	$matcher = sub{ my $in = lc($_[0]) ; for(@patterns){$in eq $_ && return 1} ; 0 }
} else {
	#Fixed match to whole info line based on array
	$matcher = sub{ for(@patterns){$_[0] eq $_ && return 1} ; 0 }
} } } elsif($opts{X}) { if(%patternhash) { if($opts{i}) {
	#Fixed case insensitive match to first word of info line based on hash
	$matcher = sub{ my $fw = $_[0]; $fw =~ s/^\s+//; $fw =~ s/$id_chopper//; $patternhash{lc($fw)} }
} else {
	#Fixed match to first word of info line based on hash
	$matcher = sub{ my $fw = $_[0]; $fw =~ s/^\s+//; $fw =~ s/$id_chopper//; $patternhash{$fw} }
} } else { if($opts{i}) {
	#Fixed case insensitive match to first word of info line based on array
	$matcher = sub{ my $fw = $_[0]; $fw =~ s/^\s+//; $fw =~ s/$id_chopper//; $fw = lc($fw); for(@patterns){$fw eq $_ && return 1} ; 0 }
} else {
	#Fixed match to first word of info line based on array
	$matcher = sub{ my $fw = $_[0]; $fw =~ s/^\s+//; $fw =~ s/$id_chopper//; for(@patterns){$fw eq $_ && return 1} ; 0 }
} } } else { if($opts{i}) {
	#Fixed case insensitive match to substring of info line based on array
	$matcher = sub{ for(@patterns){my $in = lc($_[0]); index($in, $_) >= 0 && return 1 } ; 0 }
} else {
    #Fixed match to substring of info line based on array
	$matcher = sub{ for(@patterns){index($_[0], $_) >= 0 && return 1 } ; 0 } 
} } } elsif($opts{X}) {
    #Pattern match to first word of info line (patterns will have been modified to match whole word)
	$matcher = sub{ my $fw = $_[0]; $fw =~ s/^\s+//; $fw =~ s/$id_chopper//; for(@patterns){$fw =~ $_ && return 1} ; 0 };
} else {
	#Pattern match within info line
	$matcher = sub{ for(@patterns){$_[0] =~ $_ && return 1} ; 0 };
}

#And how I want to print
#Note use of $/ (IRS) here rather then \n as this supports -C.
my $printer;
if($rewrap > 0)
{
	$printer = sub { print ">$_[0]$/";
					 while( $_[1] =~ /(.{1,$rewrap})/go ) { print "$1$/" };
				   };
}
else
{
	$printer = sub { print ">$_[0]$/$_[1]$/" };
}

#And thus how I want to process
my $process_seq;
if($opts{s})
{
	$process_seq = sub { if(&$matcher($_[1]) xor $invert) { &$printer(@_) } };
}
else
{
	$process_seq = sub { if(&$matcher($_[0]) xor $invert) { &$printer(@_) } };
}


#OK, let's go through the sequences.  Unless re-wrap is on I'm simply going to assume that a 
#header line begins with \s*> and everything until the next header is sequence.
#No checking line endings either.
if(!$rewrap)
{
	my $print = 1;
	for(<>)
	{
		/^\s*>(.*)/ and	$print = (&$matcher($1) xor $invert);
		print if $print;
	}
}
else
{
	#Here do a little more sanitization
	my $header = '';
	my $sequence = '';
	my $seqsdone = 0;

	#for lines
	#  chomp
	#  if line is blank ignore it
	#  if is_header
	#   if first_line check and save
	#	elsif seq_is_empty throw error
	#   else process_sequence, save
	#  else 
	#   if first_line throw error
	#   else add_to_seq
	#done
	#if seq_empty
	#   if !first_line throw error
	#else process_sequence

	for(<>)
	{
		chomp;
		next unless /\S/;
		if(/^\s*>(.*)/)
		{
			my $newheader = $1;
			if(!$seqsdone){ checkcr2($header = $newheader); $seqsdone++ }
			elsif(!$sequence){ die "Missing sequence data at record $seqsdone.\n"; }
			else{ &$process_seq($header, $sequence) ; $header = $newheader ; $sequence = '' ; $seqsdone++ }
		}
		else
		{
			if(!$seqsdone){ die "First line of file does not contain a FASTA header.\n" }
			s/^\s+//; s/\s+$//;
			$sequence .= $_;
		}
	}
	if(!$sequence)
	{
		$seqsdone && die "Missing sequence data at final record $seqsdone.\n";
	}
	else
	{
		&$process_seq($header, $sequence);
	}
}

#Et voila.
1;
