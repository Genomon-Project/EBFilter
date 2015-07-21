#! /usr/bin/env python

import process_vcf
import control_count
import beta_binomial
import utils
import sys, os, subprocess, math, multiprocessing 
import vcf, pysam, numpy


def EBFilter_worker(targetMutationFile, targetBamPath, controlBamPathList, outputPath, mapping_qual_thres, base_qual_thres):

    controlFileNum = sum(1 for line in open(controlBamPathList, 'r'))

    ##########
    # generate pileup files
    process_vcf.vcf2pileup(targetMutationFile, outputPath + '.target.pileup', targetBamPath, mapping_qual_thres, base_qual_thres, False)
    process_vcf.vcf2pileup(targetMutationFile, outputPath + '.control.pileup', controlBamPathList, mapping_qual_thres, base_qual_thres, True)
    ##########

    ##########
    # load pileup files
    pos2pileup_target = {}
    pos2pileup_control = {}

    hIN = open(outputPath + '.target.pileup')
    for line in hIN:
        F = line.rstrip('\n').split('\t')
        pos2pileup_target[F[0] + '\t' + F[1]] = '\t'.join(F[3:])
    hIN.close()

    hIN = open(outputPath + '.control.pileup')
    for line in hIN:
        F = line.rstrip('\n').split('\t')
        pos2pileup_control[F[0] + '\t' + F[1]] = '\t'.join(F[3:])
    hIN.close()
    ##########

    vcf_reader = vcf.Reader(open(targetMutationFile, 'r'))
    vcf_reader.infos['EB'] = vcf.parser._Info('EB', 1, 'Float', "EBCall Score", "EBCall", "ver0.1.0")
    vcf_writer =vcf.Writer(open(outputPath, 'w'), vcf_reader)


    for vcf_record in vcf_reader:
        F_target = pos2pileup_target[str(vcf_record.CHROM) + '\t' + str(vcf_record.POS)].split('\t')
        F_pileup = pos2pileup_control[str(vcf_record.CHROM) + '\t' + str(vcf_record.POS)].split('\t')

        # obtain the mismatch numbers and depths of target sequence data for positive and negative strands
        varCounts_target_p, depthCounts_target_p, varCounts_target_n, depthCounts_target_n = control_count.varCountCheck(str(vcf_record.ALT[0]), F_target[0], F_target[1], F_target[2], base_qual_thres)

        varCounts_control_p = [0] * controlFileNum
        varCounts_control_n = [0] * controlFileNum
        depthCounts_control_p = [0] * controlFileNum
        depthCounts_control_n = [0] * controlFileNum

        # obtain the mismatch numbers and depths (for positive and negative strands) of control sequence data
        for i in range(controlFileNum):
            varCounts_control_p[i], depthCounts_control_p[i], varCounts_control_n[i], depthCounts_control_n[i] = control_count.varCountCheck(str(vcf_record.ALT[0]), F_pileup[3 * i], F_pileup[1 + 3 * i], F_pileup[2 + 3 * i], base_qual_thres)

        # estimate the beta-binomial parameters for positive and negative strands
        alpha_p, beta_p = beta_binomial.fit_beta_binomial(numpy.array(depthCounts_control_p), numpy.array(varCounts_control_p))
        alpha_n, beta_n = beta_binomial.fit_beta_binomial(numpy.array(depthCounts_control_n), numpy.array(varCounts_control_n))

        # evaluate the p-values of target mismatch numbers for positive and negative strands
        pvalue_p = beta_binomial.beta_binom_pvalue([alpha_p, beta_p], depthCounts_target_p, varCounts_target_p)
        pvalue_n = beta_binomial.beta_binom_pvalue([alpha_n, beta_n], depthCounts_target_n, varCounts_target_n)

        # perform Fisher's combination methods for integrating two p-values of positive and negative strands
        EB_pvalue = utils.fisher_combination([pvalue_p, pvalue_n])
        EB_score = 0
        if EB_pvalue < 1e-60:
            EB_score = 60
        else:
            EB_score = - round(math.log10(EB_pvalue), 3)


        # add the score and write the vcf record
        vcf_record.INFO['EB'] = EB_score
        vcf_writer.write_record(vcf_record)

    vcf_writer.close()


    # delete intermediate files
    subprocess.call(["rm", outputPath + '.target.pileup'])
    subprocess.call(["rm", outputPath + '.control.pileup'])



def main(args):

    # should add validity check for arguments
    targetMutationFile = args.targetMutationFile
    targetBamPath = args.targetBamPath
    controlBamPathList = args.controlBamPathList
    outputPath = args.outputPath

    mapping_qual_thres = args.q
    base_qual_thres = args.Q
    thread_num = args.t


    if thread_num == 1:
        # non multi-threading mode
        EBFilter_worker(targetMutationFile, targetBamPath, controlBamPathList, outputPath, mapping_qual_thres, base_qual_thres)
    else:
        # multi-threading mode
        ##########
        # partition vcf files
        process_vcf.partition_vcf(targetMutationFile, outputPath + ".tmp.input.vcf.", thread_num)

        jobs = []
        for i in range(thread_num):
            process = multiprocessing.Process(target = EBFilter_worker, args = \
                (outputPath + ".tmp.input.vcf." + str(i), targetBamPath, controlBamPathList, outputPath + "." + str(i), mapping_qual_thres, base_qual_thres))
            jobs.append(process)
            process.start()

        # wait all the jobs to be done
        for i in range(thread_num):
            jobs[i].join()

        # merge the individual results
        process_vcf.merge_vcf(outputPath + ".tmp.input.vcf.", outputPath, thread_num)

        # delete intermediate files
        for i in range(thread_num):
            subprocess.call(["rm", outputPath + ".tmp.input.vcf." + str(i)])
            subprocess.call(["rm", outputPath + "." + str(i)])



