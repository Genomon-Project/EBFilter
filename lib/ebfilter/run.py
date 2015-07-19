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

    ##########
    # generate pileup files
    process_vcf.vcf2pileup(targetMutationFile, outputPath + '.tmp.target.pileup', targetBamPath, mapping_qual_thres, base_qual_thres, False)
    process_vcf.vcf2pileup(targetMutationFile, outputPath + '.tmp.control.pileup', controlBamPathList, mapping_qual_thres, base_qual_thres, True)
    ##########
    
    ##########
    # load pileup files
    pos2pileup_target = {}
    pos2pileup_control = {}

    hIN = open(outputPath + '.tmp.target.pileup')
    for line in hIN:
        F = line.rstrip('\n').split('\t')
        pos2pileup_target[F[0] + '\t' + F[1]] = '\t'.join(F[3:])
    hIN.close()

    hIN = open(outputPath + '.tmp.control.pileup')
    for line in hIN:
        F = line.rstrip('\n').split('\t')
        pos2pileup_control[F[0] + '\t' + F[1]] = '\t'.join(F[3:])
    hIN.close()
    ##########


    vcf_reader2 = vcf.Reader(open(targetMutationFile, 'r'))

    for vcf_record in vcf_reader2:
        F_target = pos2pileup_target[str(vcf_record.CHROM) + '\t' + str(vcf_record.POS)].split('\t')
        F_pileup = pos2pileup_control[str(vcf_record.CHROM) + '\t' + str(vcf_record.POS)].split('\t')


        varCounts_target_p, depthCounts_target_p, varCounts_target_n, depthCounts_target_n = control_count.varCountCheck(str(vcf_record.ALT[0]), F_target[0], F_target[1], F_target[2], base_qual_thres)

        varCounts_control_p = [0] * controlFileNum
        varCounts_control_n = [0] * controlFileNum
        depthCounts_control_p = [0] * controlFileNum
        depthCounts_control_n = [0] * controlFileNum

        for i in range(controlFileNum):
            varCounts_control_p[i], depthCounts_control_p[i], varCounts_control_n[i], depthCounts_control_n[i] = control_count.varCountCheck(str(vcf_record.ALT[0]), F_pileup[3 * i], F_pileup[1 + 3 * i], F_pileup[2 + 3 * i], base_qual_thres)
            # print vcf_record.CHROM + '\t' + str(vcf_record.POS) + '\t' + '\t'.join(controlCounts)
        # print ','.join(str(varCounts_p)) + '\t' + ','.join(str(depthCounts_p)) + '\t' + ','.join(str(varCounts_n)) + ','.join(str(depthCounts_n))
 
        alpha_p, beta_p = beta_binomial.fit_beta_binomial(numpy.array(depthCounts_control_p), numpy.array(varCounts_control_p))
        alpha_n, beta_n = beta_binomial.fit_beta_binomial(numpy.array(depthCounts_control_n), numpy.array(varCounts_control_n))

        pvalue_p = beta_binomial.beta_binom_pvalue([alpha_p, beta_p], depthCounts_target_p, varCounts_target_p)
        pvalue_n = beta_binomial.beta_binom_pvalue([alpha_n, beta_n], depthCounts_target_n, varCounts_target_n)

        print str(alpha_p) + '\t' + str(beta_p)
        print str(alpha_n) + '\t' + str(beta_n)
        print str(pvalue_p) + '\t' + str(pvalue_n)




