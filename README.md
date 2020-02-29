# NGS2_project
Genomic evaluation through accounting for recognized variant locations is one of the aims of bioinformatics analysis. Calling of variants either during alignment or after it may represents as one of the challenging for improving speed and accuracy of sequence alignment. The previous known Burrows–Wheeler alignment (BWA) and variant identification using Genome Analysis ToolKit (GATK) and SAMtools are compared with Findmap and Findvar programs which are used during the alignment through the present study. Otherwise, due to troubleshooting, different programs may be used in the comparison such as Freebayes. Strategies using human colorectal tumor sample and a whole human reference chromosome 18 are tested. Aligning reads and counting reference and alternate alleles for DNA source, outputs potential new single nucleotide polymorphysim, insertion, and deletion alleles are carried out. Tests are assessed according to times for processing for BWA, GATK, and SAMtools on one side, and for the other program on the other in minutes. Alignment programs requiring total memory (in GB) for each program are also well-thought-out. Mapping of reads and accuracy of calling alleles for known variants are valued. Advantages of either both programs through processing rapidness, precise alignment, usefulness of data summaries, compact output, and number of steps are reported in which calling known variants more efficiently and accurately are reflected.
Key words: Sequence alignment, Variant calling, GATK, BWA, Findmap, Findvar, Freebayes
Troubleshooting 
I-
Download latest SRA-tools version (Ubuntu)
                                                Do
##
Loading Error 
sratoolkit.2.9.6-ubuntu64
II-
Downloading samples_Data

##
Loading error from SRA by fastq-dump due to low space of environment and large data from SRA site
Loading by wget by data access 
III-
Download Reference Genome 

##
Loading error when load whole reference 
download hg19 chromosome fasta files
Load hg19 and subset 5 chromosome only by command line  
IV-
Create Reference Index

##
Failed index by BWA men 
Index by BWA bwtsw for long and paired reads
V-
Align to Reference Genome

##
aligning paired end reads
bwa sampe 5_chrbwaidx DB_DS05-2-19_GTGAAA_L002_R1_009.sai DB_DS05-2-19_GTGAAA_L002_R2_001.sai DB_DS05-2-19_GTGAAA_L002_R1_009.fastq DB_DS05-2-19_GTGAAA_L002_R2_001.fastq > DB_DS05-2-19_GTGAAA_L002_pe.sam
Results: 
[bwa_sai2sam_pe_core] convert to sequence coordinate... 
[infer_isize] fail to infer insert size: too few good pairs
[bwa_sai2sam_pe_core] time elapses: 1.70 sec
[bwa_sai2sam_pe_core] changing coordinates of 0 alignments.
[bwa_sai2sam_pe_core] align unmapped mate...
[bwa_sai2sam_pe_core] time elapses: 0.05 sec
[bwa_sai2sam_pe_core] refine gapped alignments... 0.18 sec
[bwa_sai2sam_pe_core] print alignments... [bwa_sai2sam_pe_core] paired reads have different names: "HWI-ST1148:187:C48MFACXX:2:2214:1850:4532", "HWI-ST1148:187:C48MFACXX:2:1101:1465:1992"

real    0m2.611s
user    0m2.011s
sys     0m0.447s





I didn’t kwon the solve for this trouble shooting 







VI-
Generate and sorted BAM files


samtools view -H DS05-2-19_GTGAAA_L002.sorted.bam
Result: @HD     VN:1.6  SO:coordinate
@SQ     SN:chr18        LN:78077248
@SQ     SN:chr19        LN:59128983
@SQ     SN:chr20        LN:63025520
@SQ     SN:chr21        LN:48129895
@SQ     SN:chr22        LN:51304566
@PG     ID:bwa  PN:bwa  VN:0.7.17-r1188 CL:bwa sampe 4_chrbwaidx DB_DS05-2-19_GTGAAA_L002_R1_009.sai DB_DS05-2-19_GTGAAA_L002_R2_001.sai DB_DS05-2-19_GTGAAA_L002_R1_009.fastq DB_DS05-2-19_GTGAAA_L002_R2_001.fastq
Check the file
VII-
Mapping QC


or bamFile in *.sorted.bam;do
  output=${bamFile%.sorted.bam}
  samtools depth $bamFile | awk '{{sum+=$3}} END {{print "Average = ",sum/NR, "No of covered Nuc = ", NR}}' > $output.cov
  samtools flagstat $bamFile > $output.stat
done
Resuls: 
2 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
0 + 0 mapped (0.00% : N/A)
2 + 0 paired in sequencing
1 + 0 read1
1 + 0 read2
0 + 0 properly paired (0.00% : N/A)
0 + 0 with itself and mate mapped
0 + 0 singletons (0.00% : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
Try to kwon Where the error 
VIII-
Mark duplicate


for sample in *.sorted.bam;do
 name=${sample%.sorted.bam}
 java  -Xmx2g -jar $picard_path/picard.jar MarkDuplicates INPUT=$sample OUTPUT=$name.dedup.bam METRICS_FILE=$name.metrics.txt;
done
METRICS CLASS        picard.sam.DuplicationMetrics
LIBRARY UNPAIRED_READS_EXAMINED READ_PAIRS_EXAMINED     SECONDARY_OR_SUPPLEMENTARY_RDS  UNMAPPED_READS  UNPAIRED_READ_DUPLICATES        READ_PAIR_DUPLICATES    READ_PAIR_OPTICAL_DUPLICATES    PERCENT_DUPLICATION     ESTIMATED_LIBRARY_SIZE
Unknown Library 0       0       0       2       0       0       0       0       

IX-
Install GATK


conda install -c bioconda gatk4


Indexing


# samples
for sample in *.dedup.bam;do
 #name=${sample%.dedup.bam}
 java -Xmx2g -jar $picard_path/picard.jar BuildBamIndex VALIDATION_STRINGENCY=LENIENT INPUT=$sample
done


# Reference 
java -Xmx2g -jar $picard_path/picard.jar CreateSequenceDictionary R=5_chr.fa O=5_chr.dict


Download known variants

##
Download known polymorphic sites

1.
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/GATK/common_all.vcf.gz -O common_all.vcf.gz
Difficult to select specific chromosome
2.
Time gunzip common_all.vcf.gz

3.
gatk IndexFeatureFile -I common_all.vcf


Recalibrate Bases BQSR

4.
for sample in *.dedup.bam;do
name=${sample%.dedup.bam}
gatk --java-options "-Xmx2G" BaseRecalibrator \
-R 5_chr.fa -I $sample --known-sites common_all.vcf \
-O $ name.report
A USER ERROR has occurred: Number of read groups must be >= 1, but is 0
5.
gatk --java-options "-Xmx2G" ApplyBQSR \
-R 4_chr.fa -I $sample -bqsr $name.report \
-O $name.bqsr.bam --add-output-sam-program-record --emit-original-quals
done

Joint variant calling using HaplotypeCaller
## assess genotype likelihood per-sample
for sample in *.bqsr.bam;do
 name=${sample%.bqsr.bam}
gatk --java-options "-Xmx2G" HaplotypeCaller \
  -R 4_chr.fa -I $sample \
  --emit-ref-confidence GVCF \
  --pcr-indel-model NONE \
  -O $ DS05-2-19_GTGAAA_L002.gvcf
Done

## Joint Genotyping
gatk --java-options "-Xmx60G" GenotypeGVCFs \
-R 4_chr.fa \
-V DS05-2-19_GTGAAA_L002.gvcf \
--max-alternate-alleles 6 \
-O DS05-2-19_GTGAAA_L002.vcf
A USER ERROR has occurred: Argument emit-ref-confidence has a bad value: Can only be used in single sample mode currently. Use the --sample-name argument to run on a single sample out of a multi-sample BAM file
## annotated output
gatk --java-options "-Xmx60G" GenotypeGVCFs \
-R 5_chr.fa \
-V DS05-2-19_GTGAAA_L002.gvcf \
--max-alternate-alleles 6 \
--dbsnp common_all.vcf \
-O DS05-2-19_GTGAAA_L002_ann.vcf

# check how many variant got annotated
grep -v "^#" raw_vari DS05-2-19_GTGAAA_L002ants_ann.vcf | awk '{print $3}' | grep "^rs" | wc -l

VCF statistics

## index the VCF file
conda install -c bioconda tabix
bgzip -c DS05-2-19_GTGAAA_L002_ann.vcf > DS05-2-19_GTGAAA_L002_ann.vcf.gz
tabix -p vcf DS05-2-19_GTGAAA_L002_ann.vcf.gz 

conda install -c bioconda rtg-tools
rtg vcfstats DS05-2-19_GTGAAA_L002_ann.vcf.gz > stats.txt

Split SNPs and indels

gatk --java-options "-Xmx2G" SelectVariants \
-R 5_chr.fa \
-V DS05-2-19_GTGAAA_L002_ann.vcf \
--select-type-to-include SNP \
-O DS05-2-19_GTGAAA_L002_ann_SNP.vcf

gatk --java-options "-Xmx2G" SelectVariants \
-R 4_chr.fa \
-V DS05-2-19_GTGAAA_L002_ann.vcf \
--select-type-to-include INDEL \
-O DS05-2-19_GTGAAA_L002_ann_INDEL.vcf

XX- 
Results 


## Figures 
wget https://raw.githubusercontent.com/dib-lab/dogSeq/master/scripts/densityCurves.R
sudo Rscript -e "install.packages('ggplot2', contriburl=contrib.url('http://cran.r-project.org/'))"
for f in SNP.* INDEL.*;do
  Rscript densityCurves.R "$f"
done

