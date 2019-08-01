#!/usr/bin/perl -w
use strict;

use FindBin qw($Bin);
use Cwd qw(abs_path);
use Getopt::Long;
use lib $Bin;
use CanopyClustering;

my %opt = ('merge'=>'spearman,0.97', 'outdir'=>'./', 'prefix'=>'temp');
GetOptions(\%opt, 'table:s', 'list:s', 'merge:s', 'outdir:s', 'prefix:s');

my ($m2, $t2) = split /,/, $opt{merge};

($opt{table} && -s $opt{table} && $opt{list} && -s $opt{list}) || die &usage();

my $outdir = abs_path($opt{outdir});

(-d $outdir) || mkdir $outdir;

my ($all_id_num, $samples, $ids) = id_matrix_get($opt{table});

my ($all_pre_canopy, $allid_contained_nums, $all_canopy_profile) = get_hash($opt{list}, $all_id_num);

my %merged_canopy = canopies_initial_cluster($all_canopy_profile, $m2, $t2, $t2, 'Merged', 1, $allid_contained_nums);

my $outfile = $outdir . "/$opt{prefix}.merged.canopy.list";
my %args = ('hash'=>\%merged_canopy, 'output'=>$outfile, 'sample'=>'', 'header'=>1, 'prefix'=>$opt{prefix}) ;
write_files(%args);

my ($merge_extend_canopies, $merge_extend_canopies_profile) = final_canopy_profile(\%merged_canopy, $all_pre_canopy, $all_id_num);

my $outfile1 = $outdir . "/$opt{prefix}.extend.merged.canopy.list";
my %args1 = ('hash'=>$merge_extend_canopies, 'output'=>$outfile1, 'sample'=>'', 'header'=>0, 'prefix'=>$opt{prefix}) ;
write_files(%args1);

my $outfile2 = $outdir . "/$opt{prefix}.extend.merged.canopy.txt";
my %args2 = ('hash'=>$merge_extend_canopies_profile, 'output'=>$outfile2, 'sample'=>$samples, 'header'=>1, 'prefix'=>$opt{prefix}) ;
write_files(%args2);


sub usage{
	print "#"x108, "\n";
	print "Description: second step script for merging split canopy list!\n\n";
	print "Usage: perl $0 --table <input.table> --list <canopy.list> --merge method,cutoff \n";
	print "\t--table *         <file>  input table\n";
	print "\t--list *          <file>  input canopy list\n";
	print "\t--merge           <str>   parameters of merging canopies, default is spearman,0.97\n";
	print "\t--prefix          <str>   the prefix of output files and ids\n\n";
	print "Example: perl $0 --table all.table.xls --list split.pre.canopy.list --merge spearman,0.97 --outdir ./\n\n";
	print '#'x108, "\n";
	exit;
}
