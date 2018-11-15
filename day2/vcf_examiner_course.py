import argparse
import os
import sys
import vcf

__author__ = 'Stephan'

## Aim of this program is to extract information from a given VCF file; arguments:
## - VCF file
##

## Define the argument parser
parser = argparse.ArgumentParser(description="Examine the content of a VCF file")

## Specify the input to pass the path of a VCF file to the program
parser.add_argument('--vcf', help='The VCF file', type=str, required=True)

class VCF_Examiner:
    
    def __init__(self, vcf_path):
        ## Check if pyvcf is installed
        print("PyVCF version: %s" % vcf.VERSION)

        ## Create a variable in the object
        self.vcf_path = vcf_path
    
    
    def get_average_quality_of_file(self):
        '''
            Get the average PHRED quality of all variants
            :return:
            '''
        ## Outline
        ## - Open the file
        ## - Introduce counter (so we know how many variants are present)
        ## - Introduce variable for sum of quality
        ## - Go over each record and fill values (use record.QUAL)
        ## - Calculate the average

        total_counter = 0
        sum_quality = 0
        with open(self.vcf_path) as my_vcf_fh:
            vcf_reader = vcf.Reader(my_vcf_fh)
            for record in vcf_reader:
                total_counter = total_counter + 1
                sum_quality = sum_quality + record.QUAL

        ## Build the average
        average = sum_quality / total_counter

        return("%.2f" % average)


    
    
    def get_total_number_of_variants_of_file(self):
        '''
            Get the total number of variants
            :return: total number of variants
            '''
        total_counter = 0
        with open(self.vcf_path) as my_vcf_fh:
            vcf_reader = vcf.Reader(my_vcf_fh)
            for record in vcf_reader:
                total_counter = total_counter + 1

        return(total_counter)
    
    
    def get_variant_caller_of_vcf(self):
        '''
            Return the variant caller name
            :return:
            '''
        with open(self.vcf_path) as my_vcf_fh:
            vcf_reader = vcf.Reader(my_vcf_fh)

            ## Print the '_header_lines'    
            for header_line in vcf_reader._header_lines:
                if header_line.startswith("##source"):
                    return(header_line.strip())
    
    
    def get_human_reference_version(self):
        '''
            Return the genome reference version
            :return:
            '''
        ## Outline
        ## - Open the file
        ## - Look for "##reference" in the header_lines
        ## - Return the first found entry

        with open(self.vcf_path) as my_vcf_fh:
            vcf_reader = vcf.Reader(my_vcf_fh)

            ## Print the '_header_lines'    
            for header_line in vcf_reader._header_lines:
                if header_line.startswith("##reference"):
                    return(header_line.strip().split("/")[-1])


    
    
    def get_number_of_indels(self):
        '''
            Return the number of identified INDELs
            :return:
            '''
        indel_counter = 0
        with open(self.vcf_path) as my_vcf_fh:
            vcf_reader = vcf.Reader(my_vcf_fh)
            for record in vcf_reader:
                if record.is_indel:
                    indel_counter = indel_counter + 1

        return(indel_counter)
    
    
    def get_number_of_snvs(self):
        '''
            Return the number of SNVs
            :return:
            '''
        snv_counter = 0
        with open(self.vcf_path) as my_vcf_fh:
            vcf_reader = vcf.Reader(my_vcf_fh)
            for record in vcf_reader:
                if record.is_snp:
                    snv_counter = snv_counter + 1

        return(snv_counter)

    
    def get_number_of_heterozygous_variants(self):
        '''
            Return the number of heterozygous variants
            :return:
            '''
        ## Outline
        ## - open the file
        ## - introduce counter
        ## - use num_het
        het_counter = 0
        with open(self.vcf_path) as my_vcf_fh:
            vcf_reader = vcf.Reader(my_vcf_fh)
            for record in vcf_reader:
                het_counter = het_counter + record.num_het

        return(het_counter)
        
    
    def print_summary(self):
        print("Results:")


        print("Variant caller: " + str(self.get_variant_caller_of_vcf()))
        print("Reference: " + str(self.get_human_reference_version())) 
        print("Number of variants: " + str(self.get_total_number_of_variants_of_file()))
        print("Number of INDELs: " + str(self.get_number_of_indels()))
        print("Number of SNVs: " + str(self.get_number_of_snvs()))
        print("Number of het: " + str(self.get_number_of_heterozygous_variants()))
        print("Avg QUAL of variants: " + str(self.get_average_quality_of_file()))
        


if __name__ == '__main__':
    print("VCF_Examiner")

    ## Parse the arguments
    args = parser.parse_args()
    
    ## Check that a valid VCF file was given
    if args.vcf is not None and os.path.exists(args.vcf):
        print("Found VCF file: " + str(args.vcf))
    else:
        print("VCF file not present")
        sys.exit(1)
        
    ## Pass the argument to the VCF_Examiner object

    vcf_examiner = VCF_Examiner(args.vcf)
    vcf_examiner.print_summary()
    
    print("All done.")



