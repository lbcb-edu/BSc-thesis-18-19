#!/usr/bin/perl

use warnings;
use strict;
use feature 'say';

# location where fragmented chromosomes/plasmids are stored
my $dir = $ENV{"HOME"} . "/database/class-chromosome/";

my @files;
my $file;
my @lines;
my $n = 0;

opendir my $d, $dir
	or die "Cannot open directory: $dir";
@files = readdir $d;
closedir $d;

foreach $file (@files) {
	next if ($file eq "." || $file eq "..");
	
	open L, $dir . $file;
	chomp(@lines = <L>);
	close L;
	
	$n += scalar @lines / 2;
}

say "Number of fragments: $n";
