#!/usr/bin/perl -w
use strict;

use FindBin qw($Bin);
use Cwd qw(abs_path);
use Getopt::Long;
use lib $Bin;
use CanopyClustering;

my %opt = ('min_num'=>50, 'top_share'=>0.9, 'correlation_cut'=>0.7, 'outdir'=>'./');
GetOptions(\%opt, 'table:s', "list:s", 'min_num:n', 'top_share:f', 
			'correlation_cut:f', 'outdir:s');

(($opt{table} && -s $opt{table}) || ($opt{list} && -s $opt{list}) )|| die &usage();

my $outdir = abs_path($opt{outdir});

(-d $outdir) || mkdir $outdir;

my ($all_id_num, $samples, $ids) = id_matrix_get($opt{table});

my ($all_extend_canopy, $extend_contained_nums) = get_hash($opt{list});

my ($filtered_canopies, $filterd_canopy_profile) = filter_genes_and_canopies($all_extend_canopy, $all_id_num, $opt{min_num}, $opt{top_share}, $opt{correlation_cut});

my $outfile1 = $outdir . '/good.canopy.list';
my %args1 = ('hash'=>$filtered_canopies, 'output'=>$outfile1, 'sample'=>'', 'prefix'=>'', 'header'=>1) ;
write_files(%args1);

my $outfile2 = $outdir . '/canopy.table.xls';
my %args2 = ('hash'=>$filterd_canopy_profile, 'output'=>$outfile2, 'sample'=>$samples, 'prefix'=>'', 'header'=>1) ;
write_files(%args2);

sub usage{
	print "#"x108, "\n\n";
	print "Description: forth step script for filtering low-quality canopies\n\n";
	print "Usage: perl $0 --table <input.table> --list <canopy.list> [options]\n";
	print "\t--table *         <file>  input table\n";
	print "\t--list *          <file>  input canopy list\n";
	print "\t--min_num         <num>   filter parameters, if one canopy contains < 50 genes, remove this canopy\n";
	print "\t--top_share       <float> filter parameters, if top 3 canopy profile >= 0.9, remove this canopy\n";
	print "\t--correlation_cut <float> filter parameters, if  gene-canopy correlation < 0.6, remove this gene\n\n";
	print "#"x108, "\n";
	exit;
}
