#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Path;

my ($inputdir,$genomedir,$bismarkdir,$outputdir,$threads,$help);

GetOptions(
	"inputdir=s" => \$inputdir,
    "genomedir=s" => \$genomedir,
    "bismarkdir=s" => \$bismarkdir,
    "outputdir=s" => \$outputdir,
    "threads=s" => \$threads,
	"help!" => \$help,
);

my @samples = `find $inputdir -name "*_paired_R1.fastq.gz"`;
print join("\n",@samples)."\n";
foreach my $sample_p1 (@samples){
	chomp $sample_p1;
	$sample_p1 =~ /.*\/(.*)\_paired_R1.fastq.gz/;
    my $sample_id = $1;
    
	if(!-e "$outputdir/$bismarkdir/$sample_id"){
		mkpath("$outputdir/$bismarkdir/$sample_id",0644);
		if($@){
			print "Make path $outputdir/$bismarkdir/$sample_id failed:\n$@";
			exit(1);
		}
	}
    
	open(SH,">$outputdir/$bismarkdir/$sample_id/${sample_id}_methyanaly.sh") or die "$!\n";
    if(!-e "$inputdir/$sample_id/$sample_id\_1.fastq" || -z "$inputdir/$sample_id/$sample_id\_1.fastq"){
        print SH "bismark --parallel $threads -p $threads --rg_id $sample_id --rg_sample $sample_id --genome $genomedir -1 $inputdir/$sample_id\_paired_R1.fastq.gz -2 $inputdir/$sample_id\_paired_R2.fastq.gz -o $outputdir/$bismarkdir/$sample_id\n";
        print SH "samtools sort -n -@ $threads -o $outputdir/$bismarkdir/$sample_id/$sample_id\_paired_R1_bismark_bt2_pe_sorted.bam $outputdir/$bismarkdir/$sample_id/$sample_id\_paired_R1_bismark_bt2_pe.bam\n";
        print SH "deduplicate_bismark -p --output_dir $outputdir/$bismarkdir/$sample_id --outfile $sample_id\_paired_R1_bismark_bt2_pe_sorted $outputdir/$bismarkdir/$sample_id/$sample_id\_paired_R1_bismark_bt2_pe_sorted.bam\n";
        print SH "bismark_methylation_extractor --bedGraph --merge_non_CpG --comprehensive $outputdir/$bismarkdir/$sample_id/$sample_id\_paired_R1_bismark_bt2_pe_sorted.deduplicated.bam -o $outputdir/$bismarkdir/$sample_id\n";
    }
    close SH;
    
	my $out = system("sh $outputdir/$bismarkdir/$sample_id/${sample_id}_methyanaly.sh 1>>$outputdir/$bismarkdir/$sample_id/std.log 2>>$outputdir/$bismarkdir/$sample_id/error.log &");
	if($out==0){
		print "The task of $sample_id is successfully submitted\n";
	}
}

# perl fig1f_step1_dnameth.pl --inputdir /public/ZhangJinLab/project_enhancer/ncbi-WGS-DNA-seq-gaoshaolong --genomedir /public/ZhangJinLab/project_enhancer/annodata/GRCm38 --outputdir /public/ZhangJinLab/project_enhancer/ncbi-WGS-DNA-seq-gaoshaolong --bismarkdir bismarkfile --threads 2
