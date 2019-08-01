package CanopyClustering; 
require Exporter;
@ISA = qw(Exporter);

@EXPORT = qw(get_hash id_matrix_get canopies_initial_cluster median_all_canopy_number 
			write_files final_canopy_profile filter_genes_and_canopies);

BEGIN {push @INC, "/TJPROJ1/MICRO/wangpeng/software/perl_module/Statistics-RankCorrelation-0.1205/lib"};

### sub-routine 1 ###

sub get_hash{
	my ($list, $id_array) = @_;
	my (%id_hash, %id_nums, %canopy_profile);
	my $i = 1;
	my $t0 = time();
	open IN, $list;
#	<IN>;
	while(<IN>){
		chomp;
		my @ll = split /\t+/;
		my @genes = split /,/, $ll[2];
		$id_hash{$ll[0]} = \@genes;
		$id_nums{$ll[0]} = $ll[1];
		if($id_array){
			my @temp_profile = median_each_canopy_number(\@genes, $id_array);
			$canopy_profile{$ll[0]} = \@temp_profile;
		}
	}
	my $timeuse = time()-$t0;
	print "\nGet hash list finished ($timeuse s)!\n";
	if($id_array){
		return (\%id_hash, \%id_nums, \%canopy_profile);
	}else{
		return (\%id_hash, \%id_nums);
	}
	
}

sub id_matrix_get{
	$| = 1;
	my ($table) = @_;
	open IN, $table;
	my (@ids, %id_array);
	my $h = <IN>; chomp $h;
	my @samples = split /\t+/, $h;
	shift @samples;
	while(<IN>){
		chomp;
		my @ll = split /\t/;
		my $id = shift @ll;
		push @ids, $id;
		$id_array{$id} = \@ll;
	}
	close IN;
	print "Reading matrix finished!\n";
	return (\%id_array, \@samples, \@ids);
}

#sub routine 2###
sub canopies_initial_cluster{
	##########################################
	use Statistics::Basic qw(:all);
	use Statistics::RankCorrelation qw(:all);
	##########################################
	$| = 1;
	my ($idarr, $cor, $t1, $t2, $pre, $sign, $type) = @_;
#	\%id_array | spearman or pearson | threshold1 | threshold2 | prefix | 1
	print "\nStart $pre canopies:\n";
	my $t0 = time();
	my (%res_canopy);
	my %id_array = %{$idarr};
	my $num_ele = scalar(keys %id_array);
	my $all_num = $num_ele;
	while($num_ele > 0 && %id_array){
		my $key = $pre . '-' . $sign;
		my @allids = auto_sort(keys %id_array);
		srand(0);
		my $rand_integer = int(rand(scalar(@allids)));
		my $seed = $allids[$rand_integer];
		my $first_array = $id_array{$seed};
		delete $id_array{$seed};
		push @{$res_canopy{$key}}, $seed;
		$num_ele = scalar(@allids);
		for my $id (keys %id_array){
			my $temp_array = $id_array{$id};
			if($type){
				my %id_contain_nums = %{$type};
				if($id_contain_nums{$id} > 5000){
					push @{$res_canopy{$key}}, $id;
					delete $id_array{$id};
					next;
				}
			}
			my $correlation= $cor eq 'pearson' ? 
							 Statistics::Basic::Correlation->new($first_array, $temp_array) : 
							 Statistics::RankCorrelation->new($first_array, $temp_array)->spearman();
			if($correlation >= $t1 && $correlation < $t2){
				push @{$res_canopy{$key}}, $id;
			}elsif($correlation>=$t2){
				push @{$res_canopy{$key}}, $id;
				delete $id_array{$id};
			}
		}
		my $timeuse = sprintf("%.1f", (time()-$t0)/60);
		print "$sign: ";
		$sign++;
		print "using time $timeuse min.\r";
	}
	print "\n";
	return(%res_canopy);

}


###sub-routine 3 ###
sub median_all_canopy_number{
	my ($initial_canopies, $id_array) = @_;
	my %canopies_center;
	my %ini_canopy = %{$initial_canopies};
	my %idNums = %{$id_array};
	for my $i (keys %ini_canopy){
		my @ids = @{$ini_canopy{$i}};
		push @{$canopies_center{$i}}, median_each_canopy_number(\@ids, $id_array);
	}
	return %canopies_center;	
}

sub final_canopy_profile{
	my ($canopy1, $canopy0, $id_array) = @_;
	my %Canopy1 = %{$canopy1};
	my %Canopy0 = %{$canopy0};
	my (@allids, %final_canopies);
	for my $c (keys %Canopy1){
		my @every_canopy_contents = @{$Canopy1{$c}};
		@allids = ();
		for my $each (@every_canopy_contents){
			my @all_each_ids = @{$Canopy0{$each}};
			@allids = (@allids, @all_each_ids);
		}
		my @temp = uniq_array(@allids);
		$final_canopies{$c} = \@temp;
	}
	my %final_canopy_profile = median_all_canopy_number(\%final_canopies, $id_array);
	return (\%final_canopies, \%final_canopy_profile);
}

sub filter_genes_and_canopies{
	##########################################
	use List::Util qw(sum max);
	use Statistics::Basic qw(:all);
	use Statistics::RankCorrelation qw(:all);
	##########################################
	$| = 1;
	my $t0 = time(); 
	print "\nStart filtering:\n";
	my ($canopies, $id_array, $min_geneNum, $top_share, $cor_cut) = @_;
	my $sign = 1;
	my (%filter_canopies, %filter_canopy_profile, %canopy_id_cor, %ids_canopy);
	my %FinalCanopies = %{$canopies};
	my $total = scalar(keys %FinalCanopies);
	my %IDArray = %{$id_array};
	for my $each_canopy_id (auto_sort(keys %FinalCanopies)){
		my @temp_each_all_ids = my @each_all_ids = auto_sort(@{$FinalCanopies{$each_canopy_id}});
		scalar(@each_all_ids) < $min_geneNum && next; 
		my @temp_canopy_profile = my @canopy_profile = median_each_canopy_number(\@each_all_ids, $id_array);
		my @sort_canopy_profile = sort {$b <=> $a} (@canopy_profile);
		sum(@sort_canopy_profile[0..2])/sum(@sort_canopy_profile) >= $top_share && next;
		for my $single_id (@each_all_ids){
			my $corvalue = Statistics::RankCorrelation->new(\@temp_canopy_profile, $IDArray{$single_id})->spearman();
			if($corvalue >= $cor_cut){
				push @{$filter_canopies{$each_canopy_id}}, $single_id;
			}else{
				@temp_each_all_ids = delete_element($single_id, \@temp_each_all_ids);
				@temp_canopy_profile = median_each_canopy_number(\@temp_each_all_ids, $id_array);
			}
		}
		$filter_canopies{$each_canopy_id} || next;
		if(scalar(@{$filter_canopies{$each_canopy_id}}) < $min_geneNum || sum(sign(@temp_canopy_profile)) < 4){
			delete($filter_canopies{$each_canopy_id});
			next;
		}
		$filter_canopy_profile{$each_canopy_id} = \@temp_canopy_profile;
	#	print scalar(@{$filter_canopies{$each_canopy_id}}), "\t", join("\t", @{$filter_canopies{$each_canopy_id}}), "\n";
		my $timeuse = sprintf("%.1f", (time()-$t0)/60);
		print "$sign: ";
		$sign++;
		print "using time $timeuse min.\r";
	}
	print "\n";
	%filter_canopies = remove_null_value_hash(%filter_canopies);
	%filter_canopy_profile = remove_null_value_hash(%filter_canopy_profile);
	my $good_num = scalar(keys %filter_canopies);
	return(\%filter_canopies, \%filter_canopy_profile);
}

sub write_files{
	#%args
	my (%args) = @_;
	my ($hash, $output, $samples, $prefix, $h) = ($args{'hash'}, $args{'output'}, $args{'sample'}, $args{'prefix'}, $args{'header'});
	open OUT, ">$output";
	my $header = $samples ? "Canopy_id\t" . join("\t", @{$samples}) : 
				"Canopy_id\tGeneNum\tGene_ids";
	$h && print OUT $header, "\n";
	if(!$samples){
		for my $id (auto_sort(keys %{$hash})){
			my @ids = @{${$hash}{$id}};
			my $genenum = scalar(@ids);
			my $pre = $prefix ? "$prefix" . "_" . "$id\t$genenum\t" . join(",", @ids) : "$id\t$genenum\t" . join(",", @ids);
			print OUT $pre, "\n";
		}
	}else{
		for my $id (auto_sort(keys %{$hash})){
			my @num = @{${$hash}{$id}};
			my $pre = $prefix ? "$prefix" . "_" . "$id\t" . join("\t", @num) : "$id\t" . join("\t", @num);
			print OUT $pre, "\n";
		}	
	}
	close OUT;
}

sub datestr{
	use POSIX qw(strftime);
	return strftime "%Y-%m-%d %H:%M:%S", localtime;
}

###internal sub-routines###
sub delete_element{
	my ($ele, $arr) = @_;
	my %temp_exists; $temp_exists{$ele} = 1;
	my @arr2 = grep { !$temp_exists{$_}} @{$arr};
	return @arr2;
}

sub remove_null_value_hash{
	my (%hash) = @_;
	my %h = map {$_, $hash{$_}} grep { $hash{$_} } keys %hash;
	return %h;
}


sub mean_med{
    my ($nums, $type) = @_;
    $type ||= 'median';
    my @gnum = @{$nums};
    my ($stat);
    if($type eq 'mean'){
    	$stat = sum(@gnum)/($#gnum+1);
    }else{
    	my @sortnum = sort {$a <=> $b} @gnum;
    	$stat = $#gnum%2 ? sum(@sortnum[(($#gnum+1)/2-1)..($#gnum+1)/2])/2 : $sortnum[$#gnum/2];
    }
    return($stat);
}

sub median_each_canopy_number{
	my ($sel_ids, $id_array) = @_;
	my @canopies_center;
	my @canopies = @{$sel_ids};
	my %idNums = %{$id_array};
	my $num_length = scalar(@{$idNums{$canopies[0]}});
	for my $s (0 .. $num_length-1){
		my @temp;
		for my $id(@canopies){
			push @temp, ${$idNums{$id}}[$s]; 
		}
		push @canopies_center, mean_med(\@temp);
	}
	return @canopies_center;	
}

sub check_number{
    my ($scalar) = @_; 
    my $tf = $scalar =~ m/^-?\d+$|^-?[0-9]+?\.\d+?$|(^-?[1-9]\.?\d*?(E-|E\+)\d+?$)/i ? 1 : 0;
    return $tf;
}

sub auto_sort{
	my (@arr) = @_;
	my @sort_arr;
	my $ele = $arr[0];
	my $tof = check_number($ele);
	if($tof){
		@sort_arr = sort {$a <=> $b} @arr;
	}else{
		@sort_arr = sort {$a cmp $b} @arr;
	}
	return @sort_arr;
}

sub sign{
    my (@num) = @_; 
    my @sign;
    for(@num){
        my $s = $_ > 0 ? 1 : 0;
        push @sign, $s; 
    }   
    return @sign;
} 

sub uniq_array{
	my (@arr) = @_;
	my %count;
	my @uniq_arr = grep { ++$count{$_} < 2; } @arr;
	return @uniq_arr;
}

1;
