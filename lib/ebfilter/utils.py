#! /usr/bin/env python

import scipy.stats
import math

def fisher_combination(pvalues):

    return 1 - scipy.stats.chi2.cdf(sum([-2 * math.log(x) for x in pvalues]), len(pvalues))

