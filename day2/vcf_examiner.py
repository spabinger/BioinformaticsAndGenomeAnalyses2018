#! /usr/bin/env python3

import vcf

__author__ = 'XXX'


class VCF_Examiner:
    
    def __init__(self):
        ## Check if pyvcf is installed
        print("PyVCF version: %s" % vcf.VERSION)
    
    
    def get_average_quality_of_son(self):
        '''
            Get the average PHRED quality of all variants
            :return:
            '''
        print("TODO")
    
    
    def get_total_number_of_variants_of_son(self):
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
    assignment1 = VCF_Examiner()
    assignment1.print_summary()
    
    print("All done.")



