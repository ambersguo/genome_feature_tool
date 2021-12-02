#!/usr/bin/perl/

####################################################
#read in gene table
open IN, "</Users/guos2/Genomes/hg19/refGenehg19_12122017nonredundant22209_w_RNAseq.txt" or die;
while ($line = <IN>) {
	push(@lines, $line);
}
close IN;

#sort refgene table by chr then by chr_cor
@sorted_lines = sort {
	@a_fields = split /\t/, $a;
	@b_fields = split /\t/, $b;
	
	$a_fields[2] cmp $b_fields[2]  # sort by chr
	||
	$a_fields[4] <=> $b_fields[4]  # then sort by coordinates
} @lines;

#create arrays for each chr
foreach $line (@sorted_lines) {
	chomp($line);
	@eles = split(/\t/, $line);
	$chr = $eles[2];
	push(@$chr, $line);
}

######################################################

#open OUTput file
open OUT, ">In1MbGeneDensity_countMonteCarlo.txt" or die;

#read in coordinates to check

for ($loopcount=1; $loopcount<1001; $loopcount++ ) {
	$samplename="MLV_HiSeq_3n5simpleCat.bed";
	RAND1k();
	print "$loopcount\t";
	FIND4EACHCOR($samplename);	
}

close OUT;


sub FIND4EACHCOR {
	$file2open = $_[0];
	@vals=();
	$med=0;
	$sum=0;
	$average=0;

	open IN, "<$file2open" or die;
	
	#check each line of integration sites, see how many genes fall within 500kb+/- region
	while ($line = <IN>) {
		#count only every 5th site
		for ($idle=1;$idle<4; $idle++) {<IN>;}
		
		chomp($line);
		@eles = split(/\t/, $line);
		$chr=$eles[0];
		$strand=$eles[1];
		$cor=$eles[2];
		$count_of_genes_in_1MB=0;
		$expval=0;
		
		$cor500kbup=$cor-500000;
		$cor500kbdown=$cor+500000;
		
		#check each gene see if it falls in the 500kb+/- region of this site
		foreach $refgene (@$chr) {
			#calculate if this gene is within $cor500kbup-$cor500kbdown;
			@annotations = split(/\t/, $refgene);
			#$refID = $annotations[1];
			#$chr = $annotations[2];
			$txstart = $annotations[4];
			$txend = $annotations[5];
			$RNAseqval = $annotations[17];
			
			if($cor500kbup<$txstart and $txend<$cor500kbdown) {
				$count_of_genes_in_1MB++;
				$expval=$expval+$RNAseqval;
			}
		}
		
		#print "$line\t$expval\t$count_of_genes_in_1MB\n";
		push (@vals, $count_of_genes_in_1MB);
		$average=$average+$count_of_genes_in_1MB;
	}
	
	#calculate median
	#sort data points
	@vals = sort {$a<=>$b} (@vals);
	#test to see if there are an even number of data points
	if( @vals % 2 == 0){
	#if even then:
	  $sum = $vals[(@vals/2)-1] + $vals[(@vals/2)];
	  $med = $sum/2;
	}
	else{
	#if odd then:
	  $med=$vals[@vals/2];
	}
	
	$average=$average/(@vals);
	
	print "$med\t$average\n";
	print OUT "median gene density in 1Mbp window\t$med\taverage gene density in 1Mbp window\t$average\n";
	close IN;
	
}

sub RAND1k {
	#get chr size information from chromInfo.txt file
	@chrhuman = (chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY);
	@chrhalias =@chrhuman;
	$totalgenome=0;

	#get chr size information from CHRINFO file
	@chrhalias =@chrhuman;
	while ($chrnum = shift(@chrhalias)) {
		open CHRINFO, "</Volumes/legacy/GenomeDatabase/hg19db/chromInfo.txt" or die;
		while ($line = <CHRINFO>) {
			if($line =~/$chrnum\t(\d+)/) {
				$chrname_size{$chrnum} = $1;
				push(@chrsize, $1);
			}
		}
		close CHRINFO;
	}

	#convert chr size to chrstart and chrstop info for each chr
	@chrhalias =@chrhuman;
	$chrstart =1;
	while ($chrname = shift(@chrhalias)){
		$chrsize = shift(@chrsize);
		$totalgenome = $totalgenome+$chrsize;
		$chrstop = $totalgenome;
		$chrboundary{$chrname} = ["$chrstart","$chrstop"]; #Hash array to store chr boundary for each chr
		#print "$chrname\t$chrstart\t$chrstop\n";
		$chrstart = $totalgenome+1;
	}
	print "Total = $totalgenome\n";

	#example show how to access hash arrays
	#foreach $chr (keys %chrboundary) {
	#	print "$chr: @{$chrboundary{$chr}}\n";
	#}

	#generate large set of random number

	srand();
	open OUTRAND, ">Rand1MHg19.txt" or die;

	for ($i=1; $i<=2000; $i++) {
		$rcor = int($totalgenome*rand()) + 1; #generate a randnum evenly distributed across genome
	
		foreach $chr (keys %chrboundary) {		#convert each randnum to chr specific coordinate
			if ($rcor>= $chrboundary{$chr}[0] and $rcor <= $chrboundary{$chr}[1]){
				$cor=$rcor - $chrboundary{$chr}[0];
				$test = "$chr"."_gap.txt";
			
				open GAP, "</Volumes/legacy/GenomeDatabase/hg19db/gap/$test" or die;  #check if the cor is in a gap.  If not set the flag and print it to output
				$flag = "yes";
				while ($line = <GAP> and $line =~/chr.*\t(\d+)\t(\d+)\t(\d+)\tN/) {
					if ($cor >= $1 and $cor<=$2) { $flag = "no"; last;}  #falls within a gap and set flag to no
				}
				close GAP;
				
				if ($flag eq "yes") {
					$j++; # set counter for how many random cors to be generated
					print OUTRAND "random\t$chr\t+\t$cor\n";
				}
			}
		}
	
		if ($j eq 1000) {last;}  #If the counter is up to the limit, then exit loop
	}

	close OUTRAND;	
	#calculate gene density in surrounding 1Mb region (500kb up and 500kb down)
}