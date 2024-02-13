#check alignment quality of fastq files with FastQC; open the .html outputs in a web browser
parallel -j 4 ~/Apps/FastQCfastqc {} ::: *.fastq.gz 

#download latest reference genome at this link: https://ftp.ensembl.org/pub/current_fasta/drosophila_melanogaster/dna/
#the version I downloaded is the 28th update of the Drosophila release 6 genome; right now it'd be in its 45th update or so
wget https://ftp.ensembl.org/pub/current_fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.28.dna.toplevel.fa.gz

#index reference genome
cd ~/References/FlyBase_Ref/
bwa index -a bwtsw Drosophila_melanogaster.BDGP6.28.dna.toplevel.fa

#create reference sequence dictionary -- this is needed for the the later steps, when you actually call variants
gatk CreateSequenceDictionary -R Drosophila_melanogaster.BDGP6.28.dna.toplevel.fa -O Drosophila_melanogaster.BDGP6.28.dna.toplevel.dict

#align fastqs with reference
cd ~/FASTQs/Whole_Genome_Sequences
for i in $(ls *.fastq.gz | rev | cut -c 13- | rev | uniq)
do
bwa mem -t 8 ~/References/FlyBase_Ref/Drosophila_melanogaster.BDGP6.28.dna.toplevel.fa ${i}_R1.fastq.gz ${i}_R2.fastq.gz | samtools view --threads 8 -S -b > ~/BWA_alignments/${i}.bam
done

#add read groups to BAM files
cd ~/BWA_alignments
mkdir ../Read_groups
for i in $(ls *.bam | rev | cut -c 5- | rev)
do
gatk AddOrReplaceReadGroups -I ${i}.bam -O ../Read_groups/${i}_read_grps.bam -SO coordinate -RGID $i -RGLB $i -RGPL Illumina -RGSM $i -RGPU $i
done 

#mark duplicate reads
cd ../Read_groups
mkdir ../MarkDuplicates
parallel -j 4 gatk MarkDuplicates -I {} -O ../MarkDuplicates/{.}_dpl.bam -M ../MarkDuplicates/{.}_metrics.txt ::: *.bam
done

#sort by coordinates [NOT NEEDED IF ALREADY SORTED BY COORDINATE WHEN ADDING READ GROUPS]
#cd ../MarkDuplicates
#mkdir ../SortSam
#parallel -j 4 gatk SortSam -I {} -SO coordinate -O ../SortSam/{.}_sorted.bam ::: *.bam

#index sorted BAM file with samtools
parallel samtools index {} ::: *.bam

#call variants 
mkdir ../HaplotypeCaller
parallel -j 4 gatk HaplotypeCaller -R ~/References/FlyBase_Ref/Drosophila_melanogaster.BDGP6.28.dna.toplevel.fa -I {} -O ../HaplotypeCaller/{.}.g.vcf.gz -ERC GVCF ::: *.bam

cd ../HaplotypeCaller
mkdir ../GenotypeGVCFs
parallel gatk GenotypeGVCFs -R ~/References/FlyBase_Ref/Drosophila_melanogaster.BDGP6.28.dna.toplevel.fa -I {} -O ../GenotypeGVCFs/{.}.vcf.gz ::: *.g.vcf.gz 