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


## Define here the argument (see above)
parser = argparse.ArgumentParser(description="Extract a FASTA file from a given region")
parser.add_argument("--ref", help='Reference FASTA file', type=str)
parser.add_argument('--chr', help='Chromosome of the region.', type=str, default="chr1")
parser.add_argument('--start', help='Start of the region.', type=int, default=1000)
parser.add_argument('--end', help='End of the region.', type=int, default=2000)

def process_input(ref, chrom, start, end):
        print("Extracting: %s %d-%d from %s " % (chrom, start, end, ref))

        ## Define the output filename
        out_file_name = chrom + "_" + str(start) + "_" + str(end) + ".fasta"
        print("Writing file called: %s" %out_file_name)

        ## Perform the subprocess call
        with open(out_file_name, "w") as outfile:
                print("Calling SAMTOOLS")
                region = chrom + ":" + str(start) + "-" + str(end)
                print("Extracting region: %s" % region)
                try:
                    subprocess.check_call(["samtools", "faidx", ref, region], stdout=outfile )
                except subprocess.CalledProcessError as e:
                    print(e)                     
                    print("ERROR!: SAMTOOLS could not extract the region")
                    print("        Please make sure that the chromosome exists")
                    sys.exit(1)


def main():

        ## Parse the arguments
        args = parser.parse_args()

        ## Print the arguments
        print("Using reference file:\t" + str(args.ref))
        print("Using chr:\t\t" + str(args.chr))
        print("Using start:\t\t" + str(args.start))
        print("Using end:\t\t" + str(args.end))

        found_error = False
        if args.ref is None or not os.path.exists(args.ref):
            print("ERROR!: Specified reference path does not exist (" + str(args.ref) + ")")
            found_error = True

        ## Check that the 'end' position is bigger or equal to the 'start' position
        if args.end < args.start:
            print("ERROR!: End cannot be smaller than the start. Please check your parameters")
            found_error = True

        ## Exit the program if an error was found
        if found_error is True:
            sys.exit(1)

        ## Process only one line if no input_file
        process_input(args.ref, args.chr, args.start, args.end)

        print("All done")
        return 0


if __name__ == "__main__":
        sys.exit(main())


