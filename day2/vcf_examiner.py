import argparse
import os
import sys
import vcf

__author__ = 'XXX'

## Aim of this program is to extract information from a given VCF file; arguments:
## - VCF file
##

## Define the argument parser
parser = argparse.ArgumentParser(description="Examine the content of a VCF file")
## Specify an input to pass the path of a VCF file to the program
#parser.add_argument('long_name_of_argument', help='help_text', default="default_value" type=type_of_the_argument)

class VCF_Examiner:
    
    def __init__(self):
        ## Check if pyvcf is installed
        print("PyVCF version: %s" % vcf.VERSION)
    
    
    def get_average_quality_of_file(self):
        '''
            Get the average PHRED quality of all variants
            :return:
            '''
        print("TODO")
    
    
    def get_total_number_of_variants_of_file(self):
        '''
            Get the total number of variants
            :return: total number of variants
            '''
        print("TODO")
    
    
    def get_variant_caller_of_vcf(self):
        '''
            Return the variant caller name
            :return:
            '''
        print("TODO")
    
    
    def get_human_reference_version(self):
        '''
            Return the genome reference version
            :return:
            '''
        print("TODO")
    
    
    def get_number_of_indels(self):
        '''
            Return the number of identified INDELs
            :return:
            '''
        print("TODO")
    
    
    def get_number_of_snvs(self):
        '''
            Return the number of SNVs
            :return:
            '''
        print("TODO")
    
    def get_number_of_heterozygous_variants(self):
        '''
            Return the number of heterozygous variants
            :return:
            '''
        print("TODO")
    
    
    def print_summary(self):
        print("Print all results here")


if __name__ == '__main__':
    print("VCF_Examiner")

    ## Read the arguments (argument parser)

        
    ## Pass the argument to the VCF_Examiner object

    vcf_examiner = VCF_Examiner()
    vcf_examiner.print_summary()
    
    print("All done.")



