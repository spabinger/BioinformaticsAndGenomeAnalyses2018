import sys
import os
import csv
import pysam
import pybedtools


__author__ = 'XXX'

##
## sam_examiner
##
## The tasks are to
## - use the files with gene_information
## - examine the given SAM/BAM file
## - output the identified information
##


class SamExaminer:
    


    def __init__(self, bam_file_path, file_name):
        ## Your gene of interest
        self.gene = "PRMT3" # You can choose any gene you are interested in, as long as it is on chr11
        
        print("Using the BAM file: " + str(bam_file_path))
        self.bam_file_path = bam_file_path
                
        print("Fetching gene coordinates")
        
        ## Check if the gene coordinates file exists
        if os.path.exists(file_name):
            with open(file_name) as fh_in:
                for row in csv.reader(fh_in, delimiter="\t"):
                    if row[0] == self.gene:
                        print("found row")
                        info = row
                        
                        self.gene = info[0]
                        self.name = info[1]
                        self.chrom = info[2]
                        self.txStart = int(info[3])
                        self.txEnd = int(info[4])
                        self.strand = info[5]
                        self.exonCount = info[6]
                        self.exonStarts = info[7] # bytestring
                        self.exonEnds = info[8] # bytestring
                    
            print("Done loading fetched data")
            
        
        print("Done with init")
                
                
    def get_sam_header(self):
        print("todo")
        
        
    def get_properly_paired_reads_of_gene(self):
        print("todo")
        
                
    def get_gene_reads_with_indels(self):
        print("todo")
        
        
    def get_number_mapped_reads(self):
        print("todo")
                
        
    def get_gene_symbol(self):
        print("todo")
        
        
    def get_region_of_gene(self):
        print("todo")
        
        
    def get_number_of_exons(self):
        print("todo")
        
    
    def print_summary(self):
        print("Print all results here")
    
        
if __name__ == '__main__':
    print("SamExaminer")
    sam_examiner = SamExaminer("", "gene_coordinates.txt") ## Change the bam file here
    
    
    sam_examiner.print_summary()
    
    

