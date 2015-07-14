#! /usr/bin/env python

import math
import numpy
import scipy.optimize

def beta_binomial_loglikelihood(params, As, Bs):

    """Calculating log-likelihood of beta-binomial distribution

    Args:
        params (List[float]): the parameter of beta distribution ([alpha, beta])
    
        As (numpy.array([int])): the counts for success
        
        Bs (numpy.array([int])): the counts of trials

    """

    alpha = params[0]    
    beta = params[1]

    ML = 0
    ML = ML + reduce(lambda a, b: a + math.lgamma(b), numpy.r_[0, As + Bs + 1])
    ML = ML - reduce(lambda a, b: a + math.lgamma(b), numpy.r_[0, As + 1])
    ML = ML - reduce(lambda a, b: a + math.lgamma(b), numpy.r_[0, Bs + 1])
    
    ML = ML - reduce(lambda a, b: a + math.lgamma(b), numpy.r_[0, As + Bs + alpha + beta])
    ML = ML + reduce(lambda a, b: a + math.lgamma(b), numpy.r_[0, As + alpha])
    ML = ML + reduce(lambda a, b: a + math.lgamma(b), numpy.r_[0, Bs + beta])

    ML = ML + len(As) * (math.lgamma(alpha + beta) - math.lgamma(alpha) - math.lgamma(beta))


    # Here, we set the penalty term of alpha and beta (0.5 is slightly arbitray...)
    ML = ML - 0.5 * math.log(alpha + beta)
    return(-ML)

 

def fit_beta_binomial(As, Bs):

    """Obtaining maximum likelihood estimator of beta-binomial distribution

    Args:
        As (numpy.array([int])): the counts for success
        
        Bs (numpy.array([int])): the counts of trials

    """

    result = scipy.optimize.fmin_l_bfgs_b(beta_binomial_loglikelihood,
                                          [20, 20],
                                          args = (As, Bs),
                                          approx_grad = True,
                                          bounds = [(0.1, 10000000), (1, 10000000)])

    return result[0]




