#! /usr/bin/env python

import pysam


def checkCount(targetBam, controlBams, mutPosList, Params):

    # mapping quality threshould
    # TH_MAPPING_QUAL = Params["TH_MAPPING_QUAL"]
    TH_MAPPING_QUAL = 20

    # base quality threshould
    TH_BASE_QUAL = Params["TH_BASE_QUAL"]

    """
    for line in pysam.mpileup("-B", "-d", "10000000", "-q", str(TH_MAPPING_QUAL), "-l", mutPosList, "-Q", "0", controlBams)

        print line

 
    """


