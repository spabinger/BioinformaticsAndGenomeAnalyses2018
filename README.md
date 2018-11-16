# Bioinformatics and Genome Analyses 2018

Here you find instructions and material to the **Tools for variant analysis of next-generation genome sequencing data** section of the **Bioinformatics and genomes analyses 2018** course.

[Course information](https://webext.pasteur.fr/tekaia/BCGAIPT2018/BCGAIPT2018_Prog.html)

**NOTE: Please note that the information might change as the course progresses!**

## Lecture Notes
Here you will find the lecture notes of the course:<br/>
* [Lecture 1](lectures/lecture1.pdf) 

## Practicals
Please use the links below to get information about each day:

* [Day1](day1)
* [Day2](day2)
* [Day3](day3)
* [Day4](day4)

#### Resources for learning Python
* [learnpython.org](https://www.learnpython.org/)
* [Python - getting started](https://www.python.org/about/gettingstarted/)
* [Python - interactive shell](https://www.python.org/shell/)
* [RealPython](https://realpython.com/)
* [Python for Biologists](https://pythonforbiologists.com/)
* [The Python Guru](https://thepythonguru.com/)


#### Helpful links

###### Genotyping
* [GATK](https://software.broadinstitute.org/gatk/)
* [Comparison of pipelines](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4619817/)

###### Somatic mutations
* [In-depth comparison of somatic point mutation callers (NatSciRep)](https://www.nature.com/articles/srep36540)
* [Lancet Tool](https://github.com/nygenome/lancet)
* [Evaluation of low-cov callers](https://www.nature.com/articles/srep43169)

###### Structural variation
* [NatRev Paper](https://www.ncbi.nlm.nih.gov/pubmed/21358748)
* [Sniffles](https://github.com/fritzsedlazeck/Sniffles)
* [MAVIS - merging, annotation, validation](http://mavis.bcgsc.ca/)
* [GRIDSS - rearrangements](https://genome.cshlp.org/content/early/2017/11/02/gr.222109.117.abstract)
* [tardis](https://github.com/BilkentCompGen/tardis)
* [SV-Bay - structural variant detection in cancer](https://github.com/InstitutCurie/SV-Bay)
* [Mapping and phasing using Nanopore - Nat](https://www.nature.com/articles/s41467-017-01343-4)

###### Annotation
* [Annovar](http://annovar.openbioinformatics.org/en/latest/)
* [SeattleSeqAnnotation](http://snp.gs.washington.edu/SeattleSeqAnnotation)
* [SVScore](http://www.github.com/lganel/SVScore)

###### Pipeline
* [bcbio](https://github.com/bcbio/bcbio-nextgen)


#### Installation of packages
Please use the following commands to install packages (will be determined during the course):

```
module load anaconda3
conda install -c bioconda pysam 
conda install -c bioconda pybedtools
conda install -c bioconda pyvcf

```

