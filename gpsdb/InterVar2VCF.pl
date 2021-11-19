#!/usr/bin/perl -w
use warnings;
use strict;

my $i2v = "/export/home/huangdaoyi/GPS_test/InterVar2VCF.txt";
my $gene_omim_file = "/export/home/huangdaoyi/GPS_test/ALL_SOURCES_ALL_FREQUENCIES_diseases_to_genes_to_phenotypes.txt";
my $list_out = $ARGV[2];

my %pseudogene = ();
open IN,"</export/home/huangdaoyi/GPS_test/pseudogene.list";
while(my $line = <IN>)
{
	chomp $line;
	$pseudogene{$line} = 1;
}
close IN;

my %gene_omim = &LoadGeneOmim($gene_omim_file);
open TAB,"<$i2v";
my $vcf_out = $ARGV[1];
open OUT,">>$vcf_out";
my @tab = ();
print OUT '##fileformat=VCFv4.2',"\n";
while(my $line = <TAB>)
{
	chomp $line;
	my @arr = split/\t/,$line;
	push @tab,$arr[0];
	print OUT qq/##INFO=<ID=$arr[0],Number=$arr[1],Type=$arr[2],Description="$arr[3]">/,"\n";
}
my %gene_score = ();
my %inter_str = ();
my %sig = ();
print OUT "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO\n";
open IN,"<$ARGV[0]";
<IN>;
while(my $line = <IN>)
{
	chomp $line;
	my @arr = split/\t/,$line;
	print OUT "$arr[0]\t$arr[1]\t$arr[9]\t$arr[3]\t$arr[4]\t1\tPASS\t";
	for(my $i=0;$i<=$#tab;$i++)
	{
		if($i == 8)
		{
			$arr[5+$i] =~ /InterVar: (.+?) PVS1/;
			my $symbol = $1;
			print OUT "$tab[$i]=$1;";
			my $score = 0;
			my @res = &score($arr[5+$i]);
			if(not exists $pseudogene{$arr[5]})
			{
				$score = $res[0];
			}
			$i = $i+1;
			print OUT "$tab[$i]=$score;";
			if(exists $gene_score{$arr[5]})
			{
				if($gene_score{$arr[5]} < $score)
				{
					$gene_score{$arr[5]} = $score;
					$inter_str{$arr[5]} = $res[1];
					$sig{$arr[5]} = $symbol;
				}
			}
			else
			{
				$gene_score{$arr[5]} = $score;
				$inter_str{$arr[5]} = $res[1];
				$sig{$arr[5]} = $symbol;
			}
		}
		elsif($i == 7)
		{
			$arr[5+$i] =~ /clinvar: (.+) $/;
			print OUT "$tab[$i]=$1;";
		}
		else
		{
			if(not $arr[5+$i])
			{
				$arr[5+$i] = '';
			}
			$arr[5+$i] =~ s/\;/,/g;
			if($arr[5+$i] eq '')
			{
				$arr[5+$i] = '.';
			}
			print OUT "$tab[$i]=$arr[5+$i];";
		}
	}
	#print $arr[4];
	print OUT "\n";
}
open LIST,">>$list_out";
my %omim_list = ();
my %inter_list = ();
for my $k(keys %gene_omim)
{
	for my $k1(keys %{$gene_omim{$k}})
	{
		if(exists $omim_list{$k})
		{
			if(exists $gene_score{$k1})
			{
				if($gene_score{$k1} > $omim_list{$k})
				{
					$omim_list{$k} = $gene_score{$k1};
					$inter_list{$k}{'str'} = $inter_str{$k1};
					$inter_list{$k}{'gene'} = $k1;
					$inter_list{$k}{'sig'} = $sig{$k1};
				}
			}
		}
		else
		{
			if(exists $gene_score{$k1})
			{
				$omim_list{$k} = $gene_score{$k1};
				$inter_list{$k}{'str'} = $inter_str{$k1};
                                $inter_list{$k}{'gene'} = $k1;
				$inter_list{$k}{'sig'} = $sig{$k1};
			}
			else
			{
				$omim_list{$k} = 0;
				$inter_list{$k}{'str'} = '.';
				$inter_list{$k}{'gene'} = '.';
				$inter_list{$k}{'sig'} = '.';
			}
		}
	}
}
for my $k(sort{$omim_list{$b}<=>$omim_list{$a}} keys %omim_list)
{
	print LIST $k,"\t";
	print LIST $omim_list{$k};
	print LIST "\t",$inter_list{$k}{'str'};
	print LIST "\t",$inter_list{$k}{'gene'};
	print LIST "\t",$inter_list{$k}{'sig'};
	print LIST "\n";
}

sub LoadGeneOmim()
{
	my $file = shift;
	open IN,"<$file";
	<IN>;
	my %map = ();
	while(my $line = <IN>)
	{
		chomp $line;
		my @arr = split/\t/,$line;
		$map{$arr[0]}{$arr[1]} = 0;
	}
	return %map;
}

sub score()
{
	my $str = shift;
	my $item = $str;
	$str =~ s/\, /\,/g;
	$str =~ s/ InterVar: (.+) PVS1/PVS1/g;
	$str =~ s/\[|\]//g;
	#print $str,"\n";
        $item =~ /InterVar: (.+) PVS1/;
        my $symbol = $1;
        $item =~ /PVS1=(.+) PS=\[(.+)\] PM=\[(.+)\] PP=\[(.+)\] BA1=(.+) BS=\[(.+)\] BP=\[(.+)\]/;
        my $psv = $1;
        my @ps = split/\,/,$2;
        pop @ps;
        my @pm = split/\,/,$3;
        pop @pm;
        my @pp = split/\,/,$4;
        pop @pp;
        my $ba = $5;
        my @bs = split/\,/,$6;
        pop @bs;
        my @bp = split/\,/,$7;
        pop @bp;
	my $psv_s  = $psv;
	my $ps_s = &sumArr(@ps);
	my $pm_s = &sumArr(@pm);
	my $pp_s = &sumArr(@pp);
	my $ba_s = $ba;
	my $bs_s = &sumArr(@bs);
	my $bp_s = &sumArr(@bp);
	#print $psv_s,"|$ps_s|$pm_s|$pp_s|$ba_s|$bs_s|$bp_s\n";
	my $score = 0;
	if($ba_s != 1)
	{
		my $upNum = $pp_s/8 + $pm_s/4 + $ps_s/2 + $psv_s - $bp_s/8 - $bs_s/2;
		my $prior_p = 0.1;
		my $odd = 350;
		my $op = $odd**$upNum;
		$score = ($op*$prior_p)/(($op-1)*$prior_p+1);
		#print "$upNum|$prior_p|$odd|$op|$score\n";
	}
        #print $symbol,"\t",$psv,"\n",join("\t",@ps),"\n",join("\t",@pm),"\n",join("\t",@pp),"\n",$ba,"\n",join("\t",@bs),"\n",join("\t",@bp),"\n";
	#print $str;
	#print $symbol,"\t",$score,"\n";
	return ($score,$str);
	
}
sub sumArr()
{
	my $sum = 0;
	map{$sum+= $_;}@_;
	return $sum;
}
