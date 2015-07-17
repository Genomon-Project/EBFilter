#! /usr/bin/env python

import vcf, os, subprocess

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


def vcf2pileup(inputFilePath, outputFilePath, controlBamPathList, mapping_qual_thres, base_qual_thres):

    vcf_reader = vcf.Reader(open(inputFilePath, 'r'))
    hOUT = open(outputFilePath, 'w')
    FNULL = open(os.devnull, 'w')

    for record in vcf_reader:

        mutReg = record.CHROM + ":" + str(record.POS) + "-" + str(record.POS)

        print ' '.join(["samtools", "mpileup", "-B", "-d", "10000000", "-q", str(mapping_qual_thres), "-Q", str(base_qual_thres), "-b", controlBamPathList, "-r", mutReg])
        # mpileup =  pysam.mpileup("-B", "-d", "10000000", "-q", str(mapping_qual_thres), "-Q", str(base_qual_thres), "-b", controlBamPathList, "-r", mutReg)
        # print >> hOUT, mpileup

        subprocess.call(["samtools", "mpileup", "-B", "-d", "10000000", "-q", str(mapping_qual_thres), "-Q", str(base_qual_thres), "-b", controlBamPathList, "-r", mutReg], stdout = hOUT, stderr = FNULL)

    FNULL.close()
    hOUT.close()


if __name__ == "__main__":
    import sys
    vcf2bed(sys.argv[1], sys.argv[2])
 
