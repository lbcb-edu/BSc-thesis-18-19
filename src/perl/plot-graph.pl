#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long;
use feature 'say';

# default settings
my $help    = "";
my $amount  = "s";
my $type    = "p";
my $query   = "CP035315.1.fna";
my $hits    = 100;
my $count   = 50;
my $cutoff  = 100;
my $graph   = "s";

GetOptions (
  'help'      => \$help,
  'amount=s'  => \$amount,
  'type=s'    => \$type,
  'query=s'   => \$query,
  'hits=i'    => \$hits,
  'count=i'   => \$count, 
  'cutoff=f'  => \$cutoff,
  'graph=s'   => \$graph
);

if ($help) {
  print "Usage: $0 [options ...] \n",
		"Will generate a graph depending on the following options:\n",
		"    -help\n",
		"        prints the usage\n",
		"    -amount\n",
		"        input 's' for single query\n",
		"        input 'm' for multiple queries\n",
		"        default: 's'\n",
		"    -type\n",
		"        input 'p' for plasmid\n",
		"        input 'c' for chromosome\n",
		"        default: 'p'\n",
		"    -query\n",
		"        input .fna file which will be used as query\n",
		"        (only in case of single query mode)\n",
		"        default: plasmid 'CP035315.1.fna'\n",
		"    -hits\n",
		"        input number of hits used in blast method\n",
		"        (only in case of multiple query mode)\n",
		"        default: 100\n",
		"    -count\n",
		"        input number of wanted queries or 'a' for all queries\n",
		"        (only in case of multiple query mode)\n",
		"        default: 50\n",
		"    -cutoff\n",
		"        input percentage(%) to filter out the hits\n",
		"        default: 100\n",
		"    -graph\n",
		"        input 's' for scatter plot\n",
		"        input 'm' for marginal histogram\n",
		"        default: 's'\n";
  exit;
}

say "Please Wait Until Graph Is Generated...";

if ($amount eq "s") {
	$count = 1;
}
if ($count eq "a") {
	$count = 0;
}

# location where plasmids/chromosomes are stored
my $dir = ($type eq "p") ? "./database/plasmid" : "./database/chromosome";

opendir my $fh, $dir or die "Cannot open directory: $dir";
my @files = readdir $fh;
@files = @files[ 2 .. $#files+1 ];
closedir $fh;

my @blast;
my $temp_file = "./temp.txt";
open W, ">$temp_file";
say W "Type Identity Hits";

foreach my $f (@files) {
	if ($amount eq "s") {
		@blast = `slc-blastn -db "Gamma_plasmids.fna" "$dir/$query" -b 10000 -v 10 -e 1 | first-only.pl -c 1`;
	}
	else {
		@blast = `slc-blastn -db "Gamma_plasmids.fna" "$dir/$f" -b $hits -v 10 -e 1 | first-only.pl -c 1`;
	}
	$count--;

	foreach (@blast) {
		chomp;
		my ($q, $s, $i, $h) = split ' ', $_;
		
		next if ($i > $cutoff);
		next if ($q eq $s);
		
		# contains information about strains which consist of chromosome and plasmids
		open C, "Gamma_plasmids_meta.txt";
		while (<C>) {
			if (/^($s).*(Plasmid).*$/) {
				say W "Plasmid $i $h";
			}
			elsif (/^($s).*(Chromosome).*$/) {
				say W "Chromosome $i $h";
			}
		}
		close C;
	}
	
	last if ($count == 0);
}

close W;

$query = substr($query, 0, -4);
`Rscript ../R/doPlot.R $amount $type $query $graph`;
`rm -f $temp_file`;
