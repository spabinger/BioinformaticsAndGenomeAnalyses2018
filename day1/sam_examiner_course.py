import sys
import os
import csv
import pysam
import pybedtools


__author__ = 'Stephan'

##
## sam_examiner
##
## The tasks are to
## - use the files with gene_information
## - examine the given SAM/BAM file
## - output the identified information
##


class SamExaminer:
    
    def __init__(self, bam_file_path, gene_coordinates_file):
        ## Your gene of interest
        self.gene = "PRMT3" # You can choose any gene you are interested in, as long as it is on chr11
        
        print("Using the BAM file: " + str(bam_file_path))
        self.bam_file_path = bam_file_path

        ## Check that BAM file exist
        if os.path.exists(self.bam_file_path):
            print("BAM file exists")
        else:
            print("ERROR!: BAM file does not exist")
                
        print("Fetching gene coordinates")
        
        ## Check if the gene coordinates file exists
        if os.path.exists(gene_coordinates_file):
            with open(gene_coordinates_file) as fh_in:
                for row in csv.reader(fh_in, delimiter="\t"):
                    if row[0] == self.gene:
                        print("found gene: " + str(row[0]))
                        self.gene = row[0]
                        self.name = row[1]
                        self.chrom = row[2]
                        self.txStart = int(row[3])
                        self.txEnd = int(row[4])
                        self.strand = row[5]
                        self.exonCount = row[6]
                        self.exonStarts = row[7] # bytestrings
                        self.exonEnds = row[8] # bytestrings
                    
            print("Done loading fetched data")
            
        
        print("Done with init")
                
                
    def get_header_hd(self):
        ## Outline
        ## - Open the BAM file
        ## - Return HD element of header

        ## Opening the BAM file
        samfile = pysam.AlignmentFile(self.bam_file_path, "rb")

        ## Return HD header
        return(samfile.header['HD'])
        
        
    def get_num_gene_prop_paired_reads(self):
        ## Outline
        ## - Open the BAM file
        ## - Fetch reads of gene
        ## - Extract reads that are properly paired
        ## - Count these reads

        ## Opening the BAM file
        samfile = pysam.AlignmentFile(self.bam_file_path, "rb")

        ## Remove 'chr' prefix as we use a GRCh reference genome
        chrom_wo_chr = self.chrom[3:]
        
        ## Counter to number of reads that are properly paired
        pp_counter = 0

        ## Get all reads mapped to the gene region 
        for read in samfile.fetch(chrom_wo_chr, self.txStart, self.txEnd):
            if read.is_proper_pair:
                pp_counter = pp_counter + 1

        return(pp_counter)
        
                
    def get_num_gene_reads_w_indels(self):
        ## Outline
        ## - Open the BAM file
        ## - Fetch reads of gene
        ## - Extract the CIGAR string
        ## - Count insertions and deletions

        ## Opening the BAM file
        samfile = pysam.AlignmentFile(self.bam_file_path, "rb")
        
        ## Remove 'chr' prefix as we use a GRCh reference genome
        chrom_wo_chr = self.chrom[3:]
        print("chrom_wo_chr: " + chrom_wo_chr)

        ## Counter to count the number of reads with an insertion or deletion
        indel_counter = 0

        ## Get all reads mapped to the gene region 
        for read in samfile.fetch(chrom_wo_chr, self.txStart, self.txEnd):
            if not read.is_unmapped:
                ## Check all entries in cigar tuples
                for cigar_tuple in read.cigar: 
                    if cigar_tuple[0] == 1 or cigar_tuple[0] == 2:
                        indel_counter = indel_counter + 1
                        
        ## Return the number of reads with an INDEL
        return(indel_counter)
        
    def get_number_mapped_reads(self):
        ## Opening the BAM file
        samfile = pysam.AlignmentFile(self.bam_file_path, "rb")
        
        ## Count all reads and the mapped reads
        counter_all = 0
        counter_mapped_reads = 0        
        for read in samfile:
            counter_all = counter_all + 1
            if not read.is_unmapped:
                counter_mapped_reads = counter_mapped_reads + 1

        #print("TESTING: number of reads: %d of %d " % (counter_mapped_reads, counter_all))

        return(counter_mapped_reads)
                
        
    def get_gene_id(self):
        return(self.name)

        
    def get_gene_symbol(self):
        return(self.gene)

        
    def get_region_of_gene(self):
        region_of_gene = "%s:%d-%d" % (self.chrom, self.txStart, self.txEnd) 
        return(region_of_gene)
        
        
    def get_number_of_exons(self):
        return(self.exonCount)
        
    
    def print_summary(self):
        print("Results:")
        print("Region of gene: " + self.get_region_of_gene())
        print("Number of exons: " + self.get_number_of_exons())
        print("Number of mapped reads: " + str(self.get_number_mapped_reads()))
        print("Gene id: " + str(self.get_gene_id()))
        print("Gene symbol: " + str(self.get_gene_symbol()))
        print("Number of gene reads with INDEL: " + str(self.get_num_gene_reads_w_indels()))
        print("Number of proper pair reads: " + str(self.get_num_gene_prop_paired_reads()))
        print("Header HD element: " + str(self.get_header_hd()))


if __name__ == '__main__':
    print("SamExaminer")

    bam_file_base = ""
    sam_examiner = SamExaminer(bam_file_base + "HG00096_chrom11_small.bam", "gene_coordinates.txt")     
    
    sam_examiner.print_summary()
    
    


