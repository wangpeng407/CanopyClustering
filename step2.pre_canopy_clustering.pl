#!/usr/bin/perl -w
use strict;

use FindBin qw($Bin);
use Cwd qw(abs_path);
use Getopt::Long;
use lib $Bin;
use CanopyClustering;

my %opt = ('cluster'=>'pearson,0.9', 'outdir'=>'./');

GetOptions(\%opt, 'table:s', 'cluster:s', 'outdir:s', 'prefix:s');

($opt{table} && -s $opt{table}) || die &usage();

(-d $opt{outdir}) || mkdir $opt{outdir};

my $outdir = abs_path($opt{outdir});

my ($m1, $t1) = split /,/, $opt{cluster};

my ($id_num, $samples, $ids) = id_matrix_get($opt{table});

my %pre_canopy = canopies_initial_cluster($id_num, $m1, $t1, $t1, $opt{prefix}, 1);

my %args;
my $outfile = $outdir . "/$opt{prefix}" . '.pre_canopy.list';
%args = ('hash'=>\%pre_canopy, 'output'=>$outfile, 'sample'=>'', 'prefix'=>'PC', 'header'=>0) ;
write_files(%args);

sub usage{
	print "#"x108, "\n\n";
	print "Description: first step script for clustring pre-canopies\n\n";
	print "Usage: perl $0 --table <all.table.xls> [options]\n";
	print "\t--table *         <file>  input table\n";
	print "\t--cluster         <str>   canopy clustering parameters, default is pearson,0.9\n";
	print "\t--prefix *        <str>   prefix of files\n\n";
	print "Example: perl $0 --table all.table.xls --cluster pearson,0.90 --prefix split0\n\n";
	print "#"x108, "\n";
	exit;
}
