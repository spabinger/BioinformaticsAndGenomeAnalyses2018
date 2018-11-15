# Day 3

In this practical you will get to know basic tools for SAM/BAM manipulation and call variants using different programs. Please check out the [Useful information](#useful-information) section.


## Data

* BAM file for Sniffles: [BAM_file_sniffles](https://hmd.ait.ac.at/bcga/na12878_pacbio_mts_ngmlr-0.2.3_mapped_small_sorted.bam)
* BAM file for variant calling: [BAM_file](https://hmd.ait.ac.at/bcga/HG00096.mapped.ILLUMINA.bwa.GBR.exome.20120522_chr1_50mio_with_chr_small.bam)
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



## Exercise 1

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
    
    Perform (expect weird output):
    - head HG00096.mapped.ILLUMINA.bwa.GBR.exome.20120522_chr1_50mio_with_chr_small.bam
    - head na12878_pacbio_mts_ngmlr-0.2.3_mapped_small_sorted.bam


__(*)__ Rename the BAM file

	mv HG00096.mapped.ILLUMINA.bwa.GBR.exome.20120522_chr1_50mio_with_chr_small.bam aln.bam
	

__(*)__ How big is the BAM file

    ls -lah aln.bam
    

__(*)__ Inspect the header of the BAM file
	
	samtools view -H aln.bam | less
    


__(*)__ View the BAM file

    samtools view -h aln.bam | less
    
    
__(*)__ How many reads are in the BAM file?<br/>
Is there another way to count the reads (check the samtools parameters)
   
    samtools view aln.bam | wc -l 
    samtools view aln.bam | grep -c "^S" 
    samtools flagstat aln.bam
    
    
__(*)__ Answer the following questions by investigating the SAM file

    * What version of the human reference assembly/genome was used to perform the alignments?
    * What version of bwa was used to align the reads?
    * What is the name of the first read?
    * At what position does the alignment of the read start?

    samtools view -H aln.bam | less
    samtools view aln.bam | less
    
    Answers:
    - NCBI37
    - bwa 0.5.9-r16      [samtools view -H aln.bam | grep "^@PG" | head -n1 ]
    - SRR077487.2342754     [samtools view  aln.bam | head -n1 | cut -f 1]
    - 10003   [samtools view  aln.bam | head -n1 | cut -f 4]

    
__(*)__ Sort the BAM file

    samtools sort -o sorted.bam -@ 2 aln.bam    [prefix with  time  to stop the duration of the command] 
    

__(*)__ Extract small file

    samtools view -h sorted.bam | head -n 50000 | samtools view -S -b - > small_sorted.bam
    or 
    samtools view -h sorted.bam | head -n 50000 | samtools view -S -b -o small_sorted.bam

    
    
__(*)__ Index the bam file
    
    samtools index sorted.bam
    samtools index small_sorted.bam
    
    
__(*)__ How many alignments contain a deletion (D)
    
    samtools view sorted.bam | cut -f 6 | grep -c "D"
    or
    samtools view sorted.bam | awk '$6 ~ /D/ {print}' | wc -l 
    
  
__(*)__ How many sequences are in the genome reference file

    samtools view -H sorted.bam | grep -c "SN:"


#### Alignment stats

    samtools flagstat sorted.bam
    samtools idxstats sorted.bam

  
#### Prepare reference genome
__(*)__ Prepare dict index *either with chr1 or with a full reference genome*
    
    ## If it is not working use
    java -jar /home/modules/picard/picard.jar instead of picard.jar
    
    java -jar -Xmx4G picard.jar CreateSequenceDictionary R=chr1.fa O=chr1.dict
    samtools dict chr1.fa -o chr1.fa.dict


__(*)__ Prepare fai index
    
    samtools faidx hg19.fasta 



## Exercise 2 - SV calling

Will will use Sniffles to call SV variants - check out the Readme page:
[https://github.com/fritzsedlazeck/Sniffles](https://github.com/fritzsedlazeck/Sniffles)

* http://samtools.github.io/hts-specs/VCFv4.3.pdf

__(*)__ Download and install Sniffles

    ## Open a new terminal
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
	
	ln -s sniffles_dir/Sniffles-master/bin/sniffles-core-1.0.10/sniffles sniffles_link
	./sniffles -m ../../test_set/reads_region.bam -v my_test.vcf
	
	## You should expect to find one deletion at 21:21492142-21492648


__(*)__ Use Sniffles to call SV variants on the provided file (see *Data*) section

	./sniffles -m na12878_pacbio_mts_ngmlr-0.2.3_mapped_small_sorted.bam -t 4 -v my_sv_results.vcf

__(*)__ Check out the identified variants

	- How many variants were found?
	  grep -v -c "^#" my_sv_results.vcf 
	  99	
	  
	- On which chromosomes did you identify variants?
	  grep -v "^#" my_sv_results.vcf | cut -f 1 | sort | uniq
	
	- What types of variants did you find?
      grep -v "^#" my_sv_results.vcf | cut -f 5 | grep "^<" | sort | uniq
      <DEL>
      <DUP>
      <INS>
      ...

	- Open the VCF and BAM file in IGV. Browse to a location of an identified variant.
	    samtools index na12878_pacbio_mts_ngmlr-0.2.3_mapped_small_sorted.bam
	    Change the VCF version from v4.3 to v4.2 in the VCF file
	    1:871662
	    
	- What would be your next steps?
	
	
## Exercise 3



#### BAM file preparations

__(*)__ Check out Picard

    Go to the Picard website and investigate the features.
    - How many methods are available?
    - What syntax/command needs to be used to specify parameters?
    Option=Value
    - In what programming language is Picard written?
    Java

__(*)__ Prepare alignment file for Picard
    
    samtools view aln.bam | awk '{if ($7 != "*") { print }}' OFS='\t' > aln_fixed_alignment.sam
    samtools view -H aln.bam > aln_fixed_header.sam
    cat aln_fixed_header.sam aln_fixed_alignment.sam > aln_fixed_all.sam
    samtools view -S -h -b aln_fixed_all.sam > aln_fixed.bam
    rm aln_fixed_header.sam
    rm aln_fixed_alignment.sam
    rm aln_fixed_all.sam


__(*)__ Sort with Picard
    
    java -Xmx4g -jar picard.jar SortSam I=aln_fixed.bam O=sorted_picard.bam SORT_ORDER=coordinate


__(*)__ Mark duplicates
     
    java -Xmx4g -jar picard.jar MarkDuplicates I=sorted_picard.bam O=dedup.bam M=metrics.txt


__(*)__ Replace ReadGroup as we already have an existing one
    
    java -Xmx4g -jar picard.jar AddOrReplaceReadGroups I=dedup.bam O=deduprg.bam RGID=group1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=sample1


__(*)__ Index with Picard
    
    java -Xmx4g -jar picard.jar BuildBamIndex I=deduprg.bam O=deduprg.bam.bai VALIDATION_STRINGENCY=SILENT


__(*)__ Collect insert size metrics
    
    module add R....
    java -Xmx4g -jar picard.jar CollectInsertSizeMetrics I=deduprg.bam O=insertSizeHistogram.txt H=insertSizeHistogram.pdf
    
    
__(*)__ View the PDF (choose any viewer you want)

    e.g. evince insertSizeHistogram.pdf


__(*)__ Questions

    * How many reads were marked as duplicated? Look for flags.
    samtools flagstat dedup.bam
    cat metrix.txt --> 30222*2 + 1577 
    
    * What are the other sorting possibilities for SortSam?
    queryname, coordinate, duplicate    
    
    * Inspect the insert size metrics histogram.


#### SAMtools variant calling

__(*)__ Call

     samtools mpileup -f chr1.fa deduprg.bam | bcftools call -c -v -o samtools.vcf

__(*)__ Investigate result

    #How many variants were called
    grep -v "^#" samtools.vcf | wc -l
    Answer: 7321
    
    #Print the variant that are between 1-300000 
    awk '!/^#/ && $2 < 20000' samtools.vcf > samtools_1_20000.vcf


#### FreeBayes variant calling

__(*)__ Check out the FreeBayes website

    - Does FreeBayes support using a target.bed file?
    - Can you use FreeBayes with multiple cores?
    - How can you ask questions to the developer/community?
    

__(*)__ Call

     module add freebayes....
     freebayes -f chr1.fa -v freebayes.vcf deduprg.bam

__(*)__ Investigate result
  
    #Perform the same procedures as done for samtools
    #Do you notice differences?
    grep -v "^#" freebayes.vcf | wc -l



#### Useful information


__(*)__ Make file executable

    chmod +x <file.name>
    
__(*)__ Extract xy from a file

    grep -i "xy" <file.name>
    
__(*)__ Extract information from a file excluding the pattern and display the first 100 lines

    grep -v "^#" <file.name> | head -n 100
	
	
	
