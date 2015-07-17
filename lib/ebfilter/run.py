#! /usr/bin/env python

import process_vcf
import control_count
import os, subprocess
import vcf, pysam

def main(args):

    # should add validity check for arguments
    targetMutationFile = args.targetMutationFile
    targetBamPath = args.targetBamPath
    controlBamPathList = args.controlBamPathList
    outputPath = args.outputPath

    controlFileNum = sum(1 for line in open(controlBamPathList, 'r'))
    mapping_qual_thres = args.q
    base_qual_thres = args.Q

    process_vcf.vcf2pileup(targetMutationFile, outputPath + '.tmp.pileup', controlBamPathList, mapping_qual_thres, base_qual_thres)


    vcf_reader2 = vcf.Reader(open(targetMutationFile, 'r'))
    hIN_PILEUP = open(outputPath + '.tmp.pileup', 'r')

    for vcf_record, pileup in zip(vcf_reader2, hIN_PILEUP):
        F_pileup = pileup.rstrip('\n').split('\t')

        controlCounts = ["0,0,0,0"] * controlFileNum
        for i in range(controlFileNum):
            controlCounts[i] = control_count.varCountCheck(str(vcf_record.ALT[0]), F_pileup[3 + 3 * i], F_pileup[4 + 3 * i], F_pileup[5 + 3 * i], base_qual_thres)
        print vcf_record.CHROM + '\t' + str(vcf_record.POS) + '\t' + '\t'.join(controlCounts)


    hIN_PILEUP.close()



