#! /usr/bin/env python

import generate_bed

def main(args):

    # should add validity check for arguments
    targetMutationFile = args.targetMutationFile
    targetBamPath = args.targetBamPath
    controlBamPathList = args.controlBamPathList
    outputPath = args.outputPath

    mapping_qual_thres = args.q
    base_qual_thres = args.Q


    generate_bed.vcf2bed(targetMutationFile, outputPath + '.mutpos.bed')


