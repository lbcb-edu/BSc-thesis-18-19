#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long;
use List::Util qw/shuffle/;

# location where chromosomes/plasmids are stored
my $chromo_dir = $ENV{"HOME"} . "/database/chromosome/";
my $plasmid_dir = $ENV{"HOME"} . "/database/plasmid/";

# location where fragmented chromosomes/plasmids are going to be stored
# afterwards used for the classification
my $class_chromo = $ENV{"HOME"} . "/database/class-chromosome/";
my $class_plasmid = $ENV{"HOME"} . "/database/class-plasmid/";

my @chromo_files;
my @plasmid_files;
my $file;
my @lines;
my $header;
my $data;
my $data_len;
my $fragment;
my $i;

#default settings
my $help = 0;
my $fragment_len = 5000;
my $chromo_num = 100;

GetOptions (
  'help' => \$help,
  'n=i'  => \$fragment_len,
  'c=i'  => \$chromo_num
);

if ($help) {
  print "Usage: $0 [options ...] \n",
		"Will fragment chromosomes and plasmids into pieces of length 'n'.\n",
		"  options:\n",
		"      -help\n",
		"        prints the usage\n",
		"      -n\n",
		"        input length of fragmented pieces\n",
		"        default: 5000\n",
		"      -c\n",
		"        input number of chromosomes to fragment\n",
		"        default: 100\n";
  exit;
}

opendir my $dir_ch, $chromo_dir
	or die "Cannot open chromosome directory!";
@chromo_files = readdir $dir_ch;
closedir $dir_ch;

opendir my $dir_pl, $plasmid_dir
	or die "Cannot open plasmid directory!";
@plasmid_files = readdir $dir_pl;
closedir $dir_pl;

# optional - randomly chosen $chromo_num number of chromosomes to be fragmented
my @shuffled_indexes = shuffle(0 .. $#chromo_files);
my @pick_indexes = @shuffled_indexes[0 .. $chromo_num-1];
@chromo_files = @chromo_files[@pick_indexes];

foreach $file (@chromo_files) {
	next if ($file eq "." || $file eq "..");
	
	open L, $chromo_dir . $file;
	chomp(@lines = <L>);
	close L;
	
	$data = $lines[1];
	$data_len = length $data;
	
	for ($i=0; $i <= int($data_len/$fragment_len); $i++) {
		$fragment = substr($data, $fragment_len * $i, $fragment_len);
		$header = '>' . $i . '_' . substr($file, 0, -4) . '	' . 'length: ' . length($fragment);
		
		open W, '>>', $class_chromo . $file;
		print W "$header\n$fragment\n";
		close W;
	}
}

foreach $file (@plasmid_files) {
	next if ($file eq "." || $file eq "..");
	
	open L, $plasmid_dir . $file;
	chomp(@lines = <L>);
	close L;
	
	$data = $lines[1];
	$data_len = length $data;
	
	for ($i=0; $i <= int($data_len/$fragment_len); $i++) {
		$fragment = substr($data, $fragment_len * $i, $fragment_len);
		$header = '>' . $i . '_' . substr($file, 0, -4) . '	' . 'length: ' . length($fragment);
		
		open W, '>>', $class_plasmid . $file;
		print W "$header\n$fragment\n";
		close W;
	}
}
