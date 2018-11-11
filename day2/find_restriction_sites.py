import sys
import os
import csv
import argparse
from Bio.Restriction import AllEnzymes
from Bio.Restriction import RestrictionBatch
from pyfaidx import Fasta
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA


__author__ = 'XY'


## Aim of this program is to find restriction enzyme cut sites (e.g., use in GBS sequencing)
## - Reference genome
## - Regions where to look for cut sites
## - List of enzymes to be considered
##
## We will use Biopython and pyfaidx
## Required as an exiting and indexed reference genome (samtools faidx)
## Check out
## - IUPACAmbiguousDNA
## - Seq
## - RestrictionBatch





parser = argparse.ArgumentParser(description="Find restriction enzyme cut sites.")
parser.add_argument('--reference', help="Path to the reference genome", type=str, required=True)
parser.add_argument('--input_file', help="BED file with regions and annotation column, file needs to have a header (chr, start, end) and be tab delimited", type=argparse.FileType('r'), required=True)
parser.add_argument('--res_enzymes', help="Restriction enzymes separated by ','", type=str, required=True)


##
## Calculates restriction cut sites for a given list of input
##


## TODO -> use this
def check_header(row):
    header_ok = True
    if row[0].lower() != "chr":
        header_ok = False
    if row[1].lower() != "start":
        header_ok = False
    if row[2].lower() != "end":
        header_ok = False

    return header_ok



def find_restriction_sites(sequence, enyzme_list):
    ## TODO Create a Seq object (use IUPACAmbiguousDNA and Seq)
    my_seq = ""

    ## TODO Create restriction batch
    rb = ""
    print("\nSequence: %s" % str(my_seq))
    
    ## Search the sequence using the restriction batch
    rb_result = None
    
    ## Parse the search result and store it in a list
    enzyme_details = []
    for enzyme in rb_result:
        if rb_result[enzyme]:
            ## TODO
            
    ##return the enzyme_details



def find_sequence(ref, chr, start, end):
    ## TODO use ref
    ## Check out https://pypi.org/project/pyfaidx/
    return "" ## TODO


def read_input(hg_file, my_input, enzyme_list):
    ## TODO Open output file - try to implement a feature that allows the user to choose the output file name
    with my_input as my_fh, open("restriction_enzymes_result.txt", "w") as my_ofh:
        
        csv_writer = csv.writer(my_ofh, delimiter="\t")

        ## Process the regions
        for row in csv.reader(my_fh, delimiter="\t"):
            
            ## TODO
            
            
            ## Write row
            csv_writer.writerow("")


def print_all_enyzmes():
    print("Printing all %d enzymes" % len(AllEnzymes))
    for x in sorted(AllEnzymes):
        print(x)

## TODO
def check_restriction_enzymes(enzyme_list):
    ## TODO
    ## Use AllEnzymes


## Check the format of the input file
def check_input_file_format(my_input):
    print("Checking the format of the input")
    with my_input as my_fh:
        for row in csv.reader(my_fh, delimiter="\t"):
            if len(row) != 3:
                return "Number of columns in the input file is wrong (%d instead of %d)" % (len(row), 3)

    return None


def main():

    ## TODO Parse the arguments
    args = parser.parse_args()

    ## TODO assign/use correct variable
    print("ref_file_path: " + str(ref_file_path))

    ## Read the ref file
    ref_file = Fasta(ref_file_path)
    
    ## TODO Input file
    my_input = ""
    
    print("MyInput: " + str(my_input))
    
    if not my_input:
        print("No input given")
        return 0
    
    ## Check the input file format
    check_result = check_input_file_format(my_input)
    if check_result is not None:
        print(check_result)
        sys.exit(1)

        
    ## Check the provided list of enzymes
    print("Got these enzymes: " + str(args.res_enzymes))
    if not args.res_enzymes:
        print("No enzymes provided")
        return 0
    

    ## Check that valid enzymes were given
    check_restriction_enzymes(enzyme_list)
    
    read_input(ref_file, my_input, enzyme_list)
        
    return 0


if __name__ == "__main__":
    sys.exit(main())



