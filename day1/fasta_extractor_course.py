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
parser.add_argument("--ref", help='Reference FASTA file', type=str)
parser.add_argument('--chr', help='Chromosome of the region.', type=str, default="chr1")
parser.add_argument('--start', help='Start of the region.', type=int, default=1000)
parser.add_argument('--end', help='End of the region.', type=int, default=2000)

def process_input(parameters, ref_file):
        print("Extracting: " + str(parameters))

        ## TODO define the output filename
        out_file_name = "" + ".fasta"

        with open(out_file_name, "w") as outfile:
                print("")
        ## TODO Subprocess call


def main():

        ## Parse the arguments
        args = parser.parse_args()

        ## Print the arguments
        print("Using reference file:\t" + str(args.ref))
        print("Using chr:\t\t" + str(args.chr))
        print("Using start:\t\t" + str(args.start))
        print("Using end:\t\t" + str(args.end))


        if args.ref is None or not os.path.exists(args.ref):
            print("ERROR!: Specified reference path does not exist (" + str(args.ref) + ")")
            sys.exit(1) 

        ## Process only one line if no input_file
        process_input("", args.ref)

        print("All done")
        return 0


if __name__ == "__main__":
        sys.exit(main())


