# Day 4

In this practical you will get to know more tools for variant calling. Furthermore, VCF files will be annotated, filtered, and displayed. Please check out the [Useful information](#useful-information) section.


## Data

* BAM file for SV calling (but we will use the BAM file from yesterday): [BAM_file](https://hmd.ait.ac.at/bcga/HG00096.mapped.ILLUMINA.bwa.GBR.exome.20120522_chr1_50mio_with_chr.bam)
* 1000g VCF file: https://hmd.ait.ac.at/bcga/1000G_phase1.indels.hg19.vcf
* Mills VCF file: https://hmd.ait.ac.at/bcga/Mills_and_1000G_gold_standard.indels.hg19.vcf *Not needed in this practical*


## Exercise 1


__(*)__ Preparation

    Download the dataset
    

#### GATK variant calling

__(*)__ Known indel sites are here specified as variables - either copy the whole path or use variables as well

    KNOWN_INDELS_1="1000G_phase1.indels.hg19.vcf"
    KNOWN_INDELS_2="Mills_and_1000G_gold_standard.indels.hg19.vcf" (We will skip this file today)


__(*)__ Prepare the VCF file

    grep "^#" 1000G_phase1.indels.hg19.vcf | head -n 33 > 1000G_first_head.vcf
    grep "contig=<ID=chr1,length" 1000G_phase1.indels.hg19.vcf > 1000G_second_head.vcf
    grep "^#" 1000G_phase1.indels.hg19.vcf | tail -n -3 > 1000G_third_head.vcf
    grep -v "^#" 1000G_phase1.indels.hg19.vcf | grep -w "chr1" > 1000G_chr1_body.vcf
    cat 1000G_first_head.vcf 1000G_second_head.vcf 1000G_third_head.vcf 1000G_chr1_body.vcf > 1000G_chr1.vcf
    rm 1000G_first_head.vcf
    rm 1000G_second_head.vcf
    rm 1000G_third_head.vcf
    rm 1000G_chr1_body.vcf
    
    If it does not work make sure that you have a correct dict for your reference genome (chr1.fa)
    Furthermore, check that the specified chr1 length in the 1000G_chr.vcf file matches your index.
    
    
__(*)__ Index the VCF file

    ./gatk IndexFeatureFile -F 1000G_chr1.vcf


__(*)__ Base quality recalibration
    
    ./gatk BaseRecalibrator -R chr1.fa -I deduprg.bam --known-sites 1000G_chr1.vcf -O recal_data_table.txt


__(*)__ Apply (former Print) to get recalibrated reads
    
    ./gatk ApplyBQSR -R chr1.fa -I deduprg.bam -bqsr recal_data_table.txt -O deduprg_recal.bam


__(*)__ Second pass of recalibration - to see changes

     ./gatk BaseRecalibrator -R chr1.fa -I deduprg_recal.bam --known-sites 1000G_chr1.vcf -O recal_data_table_after.txt


__(*)__ Generate before after plots (requires R and ggplot2)

    Fix missing R packages
    R
    Inside R call
    install.packages(c('reshape','gplots','gsalib'))
    
    ./gatk AnalyzeCovariates -before recal_data_table.txt -after recal_data_table_after.txt -plots recalibration_plots.pdf
    


__(*)__ Before variant calling - investigate the features

    * What does the option "--genotyping-mode" do?
    * Investigate the other options
    * Do you find a option to require a minimum base quality score?
    * Can you limit variant calling to a certain region?
    * Can you use multiple cores?
    
    
__(*)__ Start variant caller

    ./gatk HaplotypeCaller -R chr1.fa -I deduprg_recal.bam --genotyping-mode DISCOVERY -O gatk.vcf


__(*)__ Questions

    * Check out the before and after plots.
    * How many variants were called?



#### SeattleSeq Annotation

__(*)__ Access<br/>

    http://snp.gs.washington.edu/SeattleSeqAnnotation/

__(*)__ Annotate VCF file

    Upload one of the VCF files (e.g. freebayes.vcf)
    Specify VCF as return type
    submit
    You should receive an annotated VCF file to the specified email address
    

#### Display files in IGV

    (Download and open) IGV
    Load the BAM file and the VCF files into IGV
    Look at the mapping on Chr 1
    Check out the results of the different variant calling programs.



#### Merge VCFs and VCF stats

__(*)__ VCFlib - merge

    vcfcombine freebayes.vcf gatk.vcf samtools.vcf > vcf_lib_merged.vcf

__(*)__ VCFlib - stats - shown here for one VCF file - repeat for all 3

    vcfstats freeb_call.vcf > freeb_call_vcf_lib_stats.txt



#### Use of VCFtools

__(*)__ VCFtools

    export PERL5LIB=<full-path>/vcftools_0.1.12a/perl
    export PATH=<full-path>/tabix/tabix-0.2.6:$PATH

    ## Index (tabix) and pack files
    cp gatk.vcf gatk_tab.vcf
    bgzip gatk_tab.vcf
    tabix -p vcf gatk_tab.vcf.gz

    ## repeat for the other two VCF files

    vcf-merge freebayes_tab.vcf.gz gatk_tab.vcf.gz samtools_tab.vcf.gz > vcf_tools_merged.vcf
    vcf-stats freebayes_tab.vcf.gz > freebayes_tab_stats.txt
    
    
    All commands at once
    export PERL5LIB=/bcga2016/vcftools-0.1.14/share/perl5
    export PATH=/bcga2016/tabix-0.2.6:$PATH

    cp freebayes.vcf freebayes_tab.vcf
    bgzip freebayes_tab.vcf 
    tabix -p vcf freebayes_tab.vcf.gz 
    cp samtools.vcf samtools_tab.vcf
    bgzip samtools_tab.vcf 
    tabix -p vcf samtools_tab.vcf.gz 
    
     
    vcf-merge samtools_tab.vcf.gz freebayes_tab.vcf.gz > merged.vcf





#### Filter variants *If there is time*
__(*)__ Using vcfutils
     
     <full-path>/bcftools/bin/vcfutils.pl varFilter -Q 20 -d 5 -D 200 samtools.vcf > samtools_filtered.vcf

__(*)__ Questions
* What other parameters can you specify for filtering variants?
* How many variants were filtered?



  
    
#### Useful information

__(*)__ Determine number of cores

    cat /proc/cpuinfo  

__(*)__ Determine memory size
    
    cat /proc/meminfo

__(*)__ Make file executable

    chmod +x <file.name>
    
__(*)__ Extract information from a file

    grep -i "info" <file.name>
    
__(*)__ Extract information from a file excluding the pattern and display the first 100 lines

    grep -v "^#" <file.name> | head -n 100
    
    
    