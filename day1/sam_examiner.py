import mysql.connector

__author__ = 'XXX'

##
## sam_examiner
##
## The tasks are to
## - download gene information from UCSC
## - examine the given SAM/BAM file
## - output the identified information
##


class SamExaminer:
    


    def __init__(self):
        ## Your gene of interest
        self.gene = "PRMT3" # You can choose any gene you are interested in, as long as it is on chr11
		

    
    def fetch_gene_coordinates(self, genome_reference, file_name):
        
        print "Connecting to UCSC to fetch data"
        
        ## Open connection
        cnx = mysql.connector.connect(host='genome-mysql.cse.ucsc.edu', user='genomep', passwd='password', db=genome_reference)
        
        ## Get cursor
        cursor = cnx.cursor()
        
        ## Build query fields
        query_fields = ["refGene.name2",
                        "refGene.name",
                        "refGene.chrom",
                        "refGene.txStart",
                        "refGene.txEnd",
                        "refGene.strand",
                        "refGene.exonCount",
                        "refGene.exonStarts",
                        "refGene.exonEnds"]
        
        ## Build query
        query = "SELECT DISTINCT %s from refGene" % ",".join(query_fields)
        
        ## Execute query
        cursor.execute(query)
        
        ## Write to file
        ## TODO this may need some work 
        with open(file_name, "w") as fh:
            for row in cursor:
                fh.write(str(row) + "\n")
    
            
        ## Close cursor & connection
        cursor.close()
        cnx.close()
        
        print "Done fetching data"
                
    def get_sam_header(self):
        print "todo"
        
    def get_properly_paired_reads_of_gene(self):
        print "todo"
        
    def get_gene_reads_with_indels(self):
        print "todo"
        
    def get_number_mapped_reads(self):
        print "todo"
        
    def get_gene_symbol(self):
        print "todo"
        
    def get_region_of_gene(self):
        print "ads"
        
    def get_number_of_exons(self):
        print "ads"
    
    def print_summary(self):
        print "Print all results here"
    
        
if __name__ == '__main__':
    print "Practical 1"
    sam_examiner = SamExaminer()
    
    sam_examiner.fetch_gene_coordinates("hg19", "MYFILE.TXT") ## TODO change filename
    sam_examiner.print_summary()
    
    

