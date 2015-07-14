#! /usr/bin/env python

import sys, vcf

def vcf2bed(inputFilePath, outputFilePath):

    """Generate bed file for checking the nonmatched-control bam files

    Args:
        inputFilePath (str): the path to the input vcf file
    
        outputFilePath (str): the path to the output bed file
    """
 
    vcf_reader = vcf.Reader(open(inputFilePath, 'r'))
    hOUT = open(outputFilePath, 'w')

    for record in vcf_reader:
        print >> hOUT, record.CHROM + '\t' + str(int(record.POS) - 1) + '\t' + str(record.POS)

    # vcf_reader.close()
    hOUT.close()


if __name__ == "__main__":
    import sys
    vcf2bed(sys.argv[1], sys.argv[2])
 
