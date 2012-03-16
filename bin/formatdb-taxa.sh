#!/usr/bin/perl
use warnings;
use strict;

while(<>) {
chomp;
my($id,@seqs) = split(/\s/,$_);

print("formatdb -n $id -i \"". join(".faa ",@seqs) . ".faa\"\n");

}

