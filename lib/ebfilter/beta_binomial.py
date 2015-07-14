#! /usr/bin/env python

import math
import numpy

def beta_binomial_loglikelihood(params, data):

    alpha = params[0]    
    beta = params[1]

    As = numpy.array(data[0::2])
    Bs = numpy.array(data[1::2])

    ML = 0
    ML = ML + reduce(lambda a, b: a + math.lgamma(b), numpy.r_[0, As + Bs + 1])
    ML = ML - reduce(lambda a, b: a + math.lgamma(b), numpy.r_[0, As + 1])
    ML = ML - reduce(lambda a, b: a + math.lgamma(b), numpy.r_[0, Bs + 1])
    
    ML = ML - reduce(lambda a, b: a + math.lgamma(b), numpy.r_[0, As + Bs + alpha + beta])
    ML = ML + reduce(lambda a, b: a + math.lgamma(b), numpy.r_[0, As + alpha])
    ML = ML + reduce(lambda a, b: a + math.lgamma(b), numpy.r_[0, Bs + beta])

    ML = ML + len(As) * (math.lgamma(alpha + beta) - math.lgamma(alpha) - math.lgamma(beta))

    return(ML)

 

