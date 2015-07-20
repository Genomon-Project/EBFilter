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


def partition_vcf(inputFilePath, outputFilePrefix, partitionNum):

    vcf_reader1 = vcf.Reader(filename = inputFilePath)
    recordNum = 0
    for record in vcf_reader1:
        recordNum += 1

    partitionNum_mod = min(recordNum, partitionNum)
    eachPartitionNum = recordNum / partitionNum_mod

    currentPartition = 0
    currentRecordNum = 0
    
    vcf_reader2 = vcf.Reader(filename = inputFilePath)
    vcf_writer = vcf.Writer(open(outputFilePrefix + "0", 'w'), vcf_reader2)
    for record in vcf_reader2:
        vcf_writer.write_record(record)
        currentRecordNum += 1
        if currentRecordNum >= eachPartitionNum and currentPartition < partitionNum_mod - 1:
            currentPartition += 1
            currentRecordNum = 0
            vcf_writer.close()
            vcf_writer = vcf.Writer(open(outputFilePrefix + str(currentPartition), 'w'), vcf_reader2) 

    vcf_writer.close()

    return partitionNum_mod


def vcf2pileup(inputFilePath, outputFilePath, bamPath, mapping_qual_thres, base_qual_thres, is_multi):

    vcf_reader = vcf.Reader(open(inputFilePath, 'r'))
    hOUT = open(outputFilePath, 'w')
    FNULL = open(os.devnull, 'w')

    for record in vcf_reader:

        mutReg = record.CHROM + ":" + str(record.POS) + "-" + str(record.POS)

        print ' '.join(["samtools", "mpileup", "-B", "-d", "10000000", "-q", str(mapping_qual_thres), "-Q", str(base_qual_thres), "-b", bamPath, "-r", mutReg])
        # mpileup =  pysam.mpileup("-B", "-d", "10000000", "-q", str(mapping_qual_thres), "-Q", str(base_qual_thres), "-b", controlBamPathList, "-r", mutReg)
        # print >> hOUT, mpileup

        if is_multi == True:
            subprocess.call(["samtools", "mpileup", "-B", "-d", "10000000", "-q", str(mapping_qual_thres), "-Q", str(base_qual_thres), "-b", bamPath, "-r", mutReg], stdout = hOUT, stderr = FNULL)
        else:
            subprocess.call(["samtools", "mpileup", "-B", "-d", "10000000", "-q", str(mapping_qual_thres), "-Q", str(base_qual_thres), bamPath, "-r", mutReg], stdout = hOUT, stderr = FNULL)

    FNULL.close()
    hOUT.close()


if __name__ == "__main__":
    import sys
    vcf2bed(sys.argv[1], sys.argv[2])
 
