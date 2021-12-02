#!/usr/bin/perl/

open IN, "<Nucleosome_GM12878_island.bed" or die;


while ($line = <IN>) {
	chomp($line);
	@eles = split(/\t/, $line);
	$chr = $eles[0];
	$nuStart = $eles[1];
	$nuEnd = $eles[2];
	$newline = "$chr\t$nuStart\t$nuEnd";
	push(@lines, $newline);
	
}
close IN;

@sorted_lines = sort {
	@a_fields = split /\t/, $a;
	@b_fields = split /\t/, $b;
	
	$a_fields[1] cmp $b_fields[1]  
	||
	$a_fields[2] <=> $b_fields[2]  
} @lines;

foreach $line (@sorted_lines) {
	chomp($line);
	@eles = split(/\t/, $line);
	$chr = $eles[0];
	push(@$chr, $line);
}

open OUT, ">InNucleosome_island_count.txt" or die;

open Slist, "<samplelist.txt" or die;
while ($samplename = <Slist>) {
	chomp ($samplename);
	print "$samplename\n";
	FIND4EACHCOR($samplename);
}
close Slist;

close OUT;

sub FIND4EACHCOR {
	$file2open = $_[0];
	$in_gene_count = 0;
	$cor_count =0;
	
	open IN, "<$file2open" or die;
	<IN>;
	
	while ($line = <IN>) {
		chomp($line);
		@eles = split(/\t/, $line);
		$chr = $eles[0];
		$cor = $eles[1];
		
		$first = 0;
		$last = @$chr;
		$midpoint = int($last/2);
		
		$in_gene=bcompare($chr, $cor, $first, $midpoint,$last );
		
		$cor_count++;
		$in_gene_count = $in_gene_count + $in_gene;
	}

	print OUT "$file2open:\t$in_gene_count\t out of \t$cor_count\t in hg19 nuleosome island\n";
	
	close IN;
}


sub bcompare {
	$chr_name = $_[0];
	$pos = $_[1];
	$first_e= $_[2];
	$mid_e = $_[3];
	$last_e = $_[4];
	$inside_gene = 0;
		
	$midpointcoordi = $$chr_name[$mid_e];
	@temp = split(/\t/, $midpointcoordi);
	$midstart = $temp[1];
	$midend = $temp[2];

	if ($pos < $midstart) {
		$first_c = $first_e;
		$last_c = $mid_e-1;
		$midpoint_c = int(($last_c -$first_c)/2) + $first_e +1;
		if($midpoint_c ne $mid_e) {
			bcompare($chr_name, $pos, $first_c,$midpoint_c, $last_c);
		}
		elsif(($midpoint_c-$first_c)==1) {
			$lastpointcoordi = $$chr_name[$first_c];
			@temp = split(/\t/, $lastpointcoordi);
			$midstart = $temp[1];
			$midend = $temp[2];
			if($pos >= $midstart and $pos <= $midend) {
				$inside_gene = 1;
			}
		}
	}
	
	elsif ($pos >= $midstart and $pos <= $midend) {
		$inside_gene = 1;
	}
	
	elsif ($pos > $midend) {
		$first_c = $mid_e+1;
		$last_c = $last_e;
		$midpoint_c = int(($last_c -$first_c)/2) + $mid_e +1;
		if($midpoint_c ne $mid_e) {
			bcompare($chr_name, $pos, $first_c,$midpoint_c, $last_c);
		}
	}
	return($inside_gene);
}
