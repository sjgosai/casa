#!/usr/bin/perl 

#use warnings;
use strict;
#usage perl CRISPR_screen_ENCODE.pl <SAM mapping file> <TSV guide list> (output file_name) geneID ensembleGeneID TSS_chrom TSS_start TSS_end gene_strand("neg, pos")
#usage example perl CRISPR_screen_ENCODE.pl FADS2-rep8_low-sort_mapped.sam FADS_guides_controls.tsv FADS2_rep8_LS FADS1 ENSG00000149485 chr11 61817002 61817003 neg

my $SAM=$ARGV[0];
my $TSV=$ARGV[1];
my $name=$ARGV[2];
my $ID=$ARGV[3];
my $ENSG=$ARGV[4];
my $TSS_chrom=$ARGV[5];
my $TSS_start=$ARGV[6];
my $TSS_end=$ARGV[7];
my $TSS_strand=$ARGV[8];

my $output= $name."_ENCODE_guideQuant.bed";
open (SCORES,">$output");


if ($TSS_strand eq "pos"){
	 $TSS_strand = "+";
	}
elsif ($TSS_strand eq "neg"){
	 $TSS_strand = "-";
	}
elsif ($TSS_strand eq "NA"){
	 $TSS_strand = "NA";
	}
else{
print STDERR "please format gene strand as pos or neg (or NA). output will be + or - accordingly (or, you know, NA)";
exit;
}   

#load each Guide into a hash and save important details
my %Guides=();
my %Guide_coords=();
my %Guide_seqs=();

print "Loading Guide/Control details .\n";

if (!open(INPUT1,$TSV)){
print STDERR "Can't open guide file: ".$TSV." Make sure its there.\n";
exit;
}   

my $guideCount=0;while(my $tmp1=<INPUT1>){
	chomp ($tmp1);
	my ($ID,$coords,$seq,$spec,$offtarget,$summary)=split(/\s+/,$tmp1);
	if ($ID eq "NAME"){
		next;
	}
	else {
		my ($first,$second)=split(">",$ID);
		$Guides{$second}=0;
		$Guide_coords{$second}=$coords;
		$Guide_seqs{$second}=$seq;
		$guideCount++;
	}
}
delete $Guides{"*"};
close INPUT1;
#########
print "Loaded ".$guideCount." guides.\n";


#now load high_exp
my $readCount=0;
print "Loading mapped sam.\n";

if (!open(INPUT2,$SAM)){
print STDERR "Can't open file: ".$SAM." Make sure its there.\n";
exit;}   

while(my $tmp1=<INPUT2>){
	chomp ($tmp1);
	my ($ID,$name,$align)=split(/\s+/,$tmp1);
	if ($ID eq "\@SQ"){	}
	else{
		if (exists $Guides{$align}){
			$Guides{$align}=($Guides{$align})+1;
			$readCount++;
		}
	}
}
delete $Guides{"*"};
close INPUT2;
#########	
print "Analzing ".$guideCount." guides, to which ".$readCount." reads have been aligned.\n";


print "Formatting outputting tracks\n";

foreach my $guide (sort keys %Guides){
	my $ctrl="NA";
	my $chrom="NA";
	my ($start,$end,$strand) =0;
	
	my $guideSpacerSeq = $Guide_seqs{$guide};
	my $guideSeq = $Guide_seqs{$guide};
	my $firstLetter=  substr( $guideSpacerSeq, 0, 1 );	
	if (($firstLetter eq "G") || ($firstLetter eq "g")){
		$guideSeq = $Guide_seqs{$guide};
		}
	else{
			$guideSeq = "G".$guideSeq;
		}




	if (index($guide, "_NT") != -1) { 
		$ctrl="NT";
	}	
	elsif (index($guide, "BASSIK") != -1) { 
		$ctrl="ST";
		($chrom,$start,$end,$strand)=split(/[:-]+/,$Guide_coords{$guide});
		if ($strand ne "+"){
			$strand = "-";
			$end = $start - 1;
			$start = $start - 4;
			
		}
		else {
			$start = $end ;
			$end = $end + 3;
		}
	}	
	elsif (index($guide, "_ST") != -1) { 
		$ctrl="STT";
		($chrom,$start,$end,$strand)=split(/[:-]+/,$Guide_coords{$guide});
		my $corrected_STT_seq = substr( $guideSpacerSeq, -20);
		$Guide_seqs{$guide}= $corrected_STT_seq;	
		if ($strand ne "+"){
			$strand = "-";
			$start=$start+3;
			$end=$end +1;
			my $newCord = $chrom.":".$start."-".$end.":".$strand;
			$Guide_coords{$guide}=$newCord;
			$end = $start - 1;
			$start = $start - 4;
			
		}
		else {
			$end=$end -2;
			my $newCord = $chrom.":".$start."-".$end.":".$strand;
			$Guide_coords{$guide}=$newCord;
			$start = $end ;
			$end = $end + 3;
		}
	}	
	else{	
		($chrom,$start,$end,$strand)=split(/[:-]+/,$Guide_coords{$guide});
		if ($strand ne "+"){
			$strand = "-";
			$end = $start - 1;
			$start = $start - 4;
			
		}
		else {
			$start = $end ;
			$end = $end + 3;
		}
	}
	
	#discern control status
	if ($ctrl eq "NA"){
		my $output_name = $ID."|".$Guide_coords{$guide};
		print SCORES $chrom."\t".$start."\t".$end."\t".$output_name."\t".$Guides{$guide}."\t".$strand."\t".$Guide_coords{$guide}."\t".$TSS_chrom."\t".$TSS_start."\t".$TSS_end."\t".$TSS_strand."\t".$ID."\t".$ENSG."\t".$Guide_seqs{$guide}."\t".$guideSeq."\t"."targeting"."\t"."NA"."\n";
	}
	
	elsif ($ctrl eq "NT"){
		my $NT_num= $guide;
		$NT_num =~ s/[^0-9]//g;
		my $output_name = "NA|nt_".$NT_num;
		print SCORES "NA"."\t"."NA"."\t"."NA"."\t".$output_name."\t".$Guides{$guide}."\t"."NA"."\t"."nt_".$NT_num."\t".$TSS_chrom."\t".$TSS_start."\t".$TSS_end."\t".$TSS_strand."\t".$ID."\t".$ENSG."\t".$Guide_seqs{$guide}."\t".$guideSeq."\t"."negative_control"."\t"."NT"."\n";
	}
	
	
	elsif ($ctrl eq "STT"){
		my $output_name = "NA|".$Guide_coords{$guide};
		print SCORES $chrom."\t".$start."\t".$end."\t".$output_name."\t".$Guides{$guide}."\t".$strand."\t".$Guide_coords{$guide}."\t".$TSS_chrom."\t".$TSS_start."\t".$TSS_end."\t".$TSS_strand."\t".$ID."\t".$ENSG."\t".$Guide_seqs{$guide}."\t".$guideSeq."\t"."negative_control"."\t"."STT"."\n";
	}
	
	
	elsif ($ctrl eq "ST"){
		my $output_name = "NA|".$Guide_coords{$guide};
		print SCORES $chrom."\t".$start."\t".$end."\t".$output_name."\t".$Guides{$guide}."\t".$strand."\t".$Guide_coords{$guide}."\t".$TSS_chrom."\t".$TSS_start."\t".$TSS_end."\t".$TSS_strand."\t".$ID."\t".$ENSG."\t".$Guide_seqs{$guide}."\t".$guideSeq."\t"."negative_control"."\t"."ST"."\n";
	}

}
close SCORES;
print "ENCOD-ing is done-zo \n";















