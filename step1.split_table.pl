#!/usr/bin/perl -w
use strict;
use Cwd qw(abs_path);

my ($table, $num_parts, $outdir) = @ARGV;

@ARGV == 3 || die &usage();

(-d $outdir) || mkdir $outdir;

$outdir = abs_path($outdir);

#my ($file, $num_part, $prefix, $outdir)
split_files($table, $num_parts, 'split', $outdir);

sub split_files{
	use POSIX qw(ceil);
	my ($file, $num_part, $prefix, $outdir) = @_;
	my ($count, %sign_data, @VAR);
	$count = 1;
	open(IN, "< $file") || die $!;
	my $h = <IN>; chomp $h;
	my @samples = split /\t+/, $h;
	shift @samples;
	while(<IN>){
		chomp;
		my @ll = split /\t/;
		$sign_data{$count} = \@ll;
		$count++;
	}
	close IN;
	
	my $chuncksize = ceil($count / $num_part);
	my @ids = 1..$count-1;
	push @VAR, [ splice @ids, 0, $chuncksize ] while @ids;
	for my $i (0..$#VAR){
		my %tmp = map{$_, $sign_data{$_}} @{$VAR[$i]};
		my $out = $outdir . '/' . $prefix . "" . $i . '.table.xls';
		open OUT, '>', $out || die $!;
		print OUT $h, "\n";
		for my $key (sort {$a <=> $b} keys %tmp){
			print OUT join("\t", @{$tmp{$key}}), "\n";
		}
		close OUT;
	}
}

sub usage{
	print "#"x80, "\n\n";
	print "Description: Script for spliting files into peices!\n\n";
	print "Usage: perl $0 input.file num_of_part output_directory\n\n";
	print "Example: perl $0 table.xls 10 ./\n\n";
	print "#"x80, "\n";
	exit;
}
