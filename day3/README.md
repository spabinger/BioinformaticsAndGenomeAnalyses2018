# Day 3

In this practical you will get to know basic tools for SAM/BAM manipulation and call variants using different programs. Please check out the [Useful information](#useful-information) section.


## Data

* BAM file for Sniffles: [BAM_file_sniffles](https://hmd.ait.ac.at/bcga/na12878_pacbio_mts_ngmlr-0.2.3_mapped_small_sorted.bam)
* BAM file for SV calling: [BAM_file](https://hmd.ait.ac.at/bcga/HG00096.mapped.ILLUMINA.bwa.GBR.exome.20120522_chr1_50mio_with_chr.bam)
* Target file: [Target_file](https://hmd.ait.ac.at/bcga/target.bed)


## Required tools

* Bedtools - http://bedtools.readthedocs.io
* Freebayes - https://github.com/ekg/freebayes
* GATK - https://www.broadinstitute.org/gatk/download
* IGV - https://www.broadinstitute.org/igv/home
* Java (7) - http://www.oracle.com/technetwork/java/javase/downloads/index.html?ssSourceSiteId=otnjp
* Picard - http://broadinstitute.github.io/picard/
* SAMtools - http://samtools.sourceforge.net/
* Varscan 2 (2.3.6) - http://varscan.sourceforge.net/
* VarDict - https://github.com/AstraZeneca-NGS/VarDictJava
* VCFtools - https://vcftools.github.io/index.html
* VCFlib - https://github.com/vcflib/vcflib



## Information

* Original data from [1000 genomes](http://www.internationalgenome.org/data-portal/sample/NA10847)
* Stop GATK "calling home" - [information](http://gatkforums.broadinstitute.org/gatk/discussion/1250/what-is-phone-home-and-how-does-it-affect-me)  



## Exercise 1 - SV calling

Will will use Sniffles to call SV variants - check out the Readme page:
[https://github.com/fritzsedlazeck/Sniffles](https://github.com/fritzsedlazeck/Sniffles)

__(*)__ Download and install Sniffles

	wget https://github.com/fritzsedlazeck/Sniffles/archive/master.tar.gz -O Sniffles.tar.gz
	tar xzvf Sniffles.tar.gz
	cd Sniffles-master/
	mkdir -p build/
	cd build/
	cmake ..
	make


__(*)__ Test if it works

	cd ../bin/sniffles*
	./sniffles
	
	
__(*)__ Use Sniffles with the provided test-dataset<br/>
https://github.com/fritzsedlazeck/Sniffles/tree/master/test_set
	
	## Go into the bin directory
	./sniffles -m ../../test_set/reads_region.bam -v my_test.vcf
	
	## You should expect to find one deletion at 21:21492142-21492648


__(*)__ Use Sniffles to call SV variants on the provided file (see *Data*) section

	./sniffles -m na12878_pacbio_mts_ngmlr-0.2.3_mapped_small_sorted.bam -v my_sv_results.vcf

__(*)__ Check out the identified variants

	- How many variants were found?
	- On which chromosomes did you identify variants?
	- What types of variants did you find?
	- What would be your next steps?
	
	

## Exercise 2

We assume that we have a properly mapped BAM file from quality checked reads.
For some variant callers we use a [target file](target.bed) to shorten variant calling time.

#### Important

* After each step inspect the generated output (cat, less, head, grep, ...).
* Organize your data and results in folders.
* Check out the specified parameters. If you don't know the exact meaning, consult the help pages (-h, --help).


#### System information

    cat /proc/cpuinfo
    cat /proc/meminfo
    

#### SAMtools

__(*)__ Download the necessary files
    see Data section


__(*)__ Rename the BAM file

	mv HG00096.mapped.ILLUMINA.bwa.GBR.exome.20120522_chr1_50mio_with_chr.bam aln.bam
	

__(*)__ How big is the BAM file

    ls -lah <aln.bam>
    

__(*)__ Inspect the header of the BAM file
	
	module add samtools
    samtools ...
    samtools view -H aln.bam


__(*)__ View the BAM file

    module add samootls
    samtools view <bam.file> | less
    
    
__(*)__ How many reads are in the BAM file?<br/>
Is there another way to count the reads (check the samtools view parameters - look for -v)
   
    samtools view <file.bam> | grep -v "^#" | wc -l
    samtools flagstat <file.bam>
    
    
__(*)__ Answer the following questions by investigating the SAM file
* What version of the human assembly was used to perform the alignments?
* What version of bwa was used to align the reads?
* What is the name of the first read?
* At what position does the alignment of the read start?

    Use SAMtools for these questions
    samtools view -H aln.bam | less
    samtools view aln.bam | less

    
__(*)__ Sort the BAM file

    samtools sort -o sorted.bam -@ <THREADS> <file.bam>
    

__(*)__ Extract small file

    samtools view -h sorted.bam | head -n 50000 | samtools view -S -b -0 small_sorted.bam 
    
    
__(*)__ Index the bam file
    
    samtools index <sorted.bam>


#### Alignment stats
    samtools flagstat sorted.bam
    samtools idxstats sorted.bam


    
  
#### Prepare reference genome
__(*)__ Prepare dict index
    
    module add picard
    java -jar -Xmx2G picard-tools/picard.jar CreateSequenceDictionary R=hg19.fasta O=hg19.dict

__(*)__ Prepare fai index
    
    samtools faidx hg19.fasta 


#### BAM file preparations
__(*)__ Sort with Picard
    
    java -Xmx4g -jar /bcga2016/picard-tools/picard.jar SortSam I=aln.bam O=sorted_picard.bam SORT_ORDER=coordinate


__(*)__ Mark duplicates
     
    java -Xmx4g -jar /bcga2016/picard-tools/picard.jar MarkDuplicates I=sorted_picard.bam O=dedup.bam M=metrics.txt


__(*)__ Add ReadGroup
    
    java -Xmx4g -jar /bcga2016/picard-tools/picard.jar AddOrReplaceReadGroups I=dedup.bam O=deduprg.bam RGID=group1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=sample1


__(*)__ Index with Picard
    
    java -Xmx4g -jar /bcga2016/picard-tools/picard.jar BuildBamIndex I=deduprg.bam O=deduprg.bam.bai VALIDATION_STRINGENCY=SILENT


__(*)__ Collect insert size metrics
    
    module add R-3.2.4
    java -Xmx4g -jar /bcga2016/picard-tools/picard.jar CollectInsertSizeMetrics I=deduprg.bam O=insertSizeHistogram.txt H=insertSizeHistogram.pdf
    
    
__(*)__ View the PDF
    evince insertSizeHistogram.pdf


__(*)__ Questions
* How many reads were marked as duplicated? Look for flags.
* What are the other sorting possibilities for SortSam?
* Inspect the insert size metrics histogram.


#### SAMtools variant calling

__(*)__ Call

     module add bcftools-1.1
     samtools mpileup -uf hg19.fasta deduprg.bam | bcftools call -c -v -o samtools.vcf

__(*)__ Investigate result

    #How many variants were called
    grep -v "^#" samtools.vcf | wc -l
    #Print the variant that are between 1-300000 
    awk '!/^#/ && $2 < "300000"' samtools.vcf

#### FreeBayes variant calling

__(*)__ Call

     module add freebayes-xy
     freebayes -f hg19.fasta -t target.bed -v freebayes.vcf deduprg.bam

__(*)__ Investigate result
  
    #Perform the same procedures as done for samtools
    #Do you notice differences?
    grep -v "^#" freebayes.vcf | wc -l


#### VarDict variant calling

__(*)__ Call

     module add VarDict-1.4.5
     AF_THR="0.01" # minimum allele frequency
     /bcga2016/vardict/VarDictJava/build/install/VarDict/bin/VarDict -G hg19.fasta -f ${AF_THR} -N my_sample -b sorted.bam -z -c 1 -S 2 -E 3 -g 4 -R chr11:1-800000 | /bcga2016/vardict/VarDictJava/VarDict/teststrandbias.R | /bcga2016/vardict/VarDictJava/VarDict/var2vcf_valid.pl -N my_sample -E -f $AF_THR > vardict.vcf
     
     In case the software is not fully installed (just for the instructions)
     VarDict -G hg19.fasta -f 0.01 -N my_sample -b deduprg.bam -z -c 1 -S 2 -E 3 -g 4 -R chr11:1-800000  > vardict.var



#### Useful information


__(*)__ Make file executable

    chmod +x <file.name>
    
__(*)__ Extract xy from a file

    grep -i "xy" <file.name>
    
__(*)__ Extract information from a file excluding the pattern and display the first 100 lines

    grep -v "^#" <file.name> | head -n 100
	
	
	