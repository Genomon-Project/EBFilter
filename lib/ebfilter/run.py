#! /usr/bin/env python

import process_vcf
import control_count
import beta_binomial
import os, subprocess
import vcf, pysam, numpy

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

        varCounts_p = [0] * controlFileNum
        varCounts_n = [0] * controlFileNum
        depthCounts_p = [0] * controlFileNum
        depthCounts_n = [0] * controlFileNum

        for i in range(controlFileNum):
            varCounts_p[i], depthCounts_p[i], varCounts_n[i], depthCounts_n[i] = control_count.varCountCheck(str(vcf_record.ALT[0]), F_pileup[3 + 3 * i], F_pileup[4 + 3 * i], F_pileup[5 + 3 * i], base_qual_thres)
            # print vcf_record.CHROM + '\t' + str(vcf_record.POS) + '\t' + '\t'.join(controlCounts)
        # print ','.join(str(varCounts_p)) + '\t' + ','.join(str(depthCounts_p)) + '\t' + ','.join(str(varCounts_n)) + ','.join(str(depthCounts_n))
 
        alpha_p, beta_p = beta_binomial.fit_beta_binomial(numpy.array(depthCounts_p), numpy.array(varCounts_p))
        alpha_n, beta_n = beta_binomial.fit_beta_binomial(numpy.array(depthCounts_n), numpy.array(varCounts_n))

        pvalue_p = beta_binomial.beta_binom_pvalue([alpha_p, beta_p], 20, 5)
        pvalue_n = beta_binomial.beta_binom_pvalue([alpha_n, beta_n], 20, 5)

        print str(alpha_p) + '\t' + str(beta_p)
        print str(alpha_n) + '\t' + str(beta_n)
        print str(pvalue_p) + '\t' + str(pvalue_n)

    hIN_PILEUP.close()



