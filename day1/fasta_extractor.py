import argparse
import sys
import csv
import subprocess
import os

__author__ = 'XY'

## Aim of this program is to extract a sequence given the following arguments:
## - Reference genome
## - chr
## - start
## - end
##
## We will use SAMTOOLS for extracting the sequence
## Required as an exiting and indexed reference genome (samtools faidx)
## If enough time left, consider extending the script to use a BED file as input (multiple regions)


## TODO Define here the argument (see above)
parser = argparse.ArgumentParser(description="Extract a FASTA file from a given region")
#parser.add_argument('long_name_of_argument', help='help_text', default="default_value" type=type_of_the_argument)


def process_input(parameters, ref_file):
        print("Extracting: " + str(parameters))

        ## TODO define the output filename
        out_file_name = "" + ".fasta"

        with open(out_file_name, "w") as outfile:
            ## TODO Subprocess call


def main():

        ## Parse the arguments
        
        ## Specify the reference version
        ref_file = "XX"
        
        print("Using reference file: " + str(ref_file))

        ## Process only one line if no input_file
        process_input(row, ref_file)
        
        print("All done")
        return 0


if __name__ == "__main__":
        sys.exit(main())


