#-------------------------------------------------------------------------------
# Name:        simple
# Purpose:
#
# Author:      jackchen
#
# Created:     29/08/2011
# Copyright:   (c) jackchen 2011
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#!/usr/bin/env python

import itertools as it
import numpy as np
import time
from algorithm import *

def evaluate(alg):
    a = np.array(alg.a)
    Algorithm.evaluate(a, alg.expected.marginals, alg.expected.interactions, out=1)

def run(choices, n, m, seed):
    alg = Algorithm(n, m, seed)
    # each generation of the population
    for i in range(n):
        # for each iteraton, pre-compute the error for each attribute
        alg.generate()
        # next, run the evaluation for each of the 2^m choices
        # this should be done as an iterator, since we don't keep the results here
        for c in choices:
            # the called function should have some method of choosing the smallest error
            alg.analyze(c)
        # explicitly call "update" to update the algorithm
        alg.update()
        if (i+1)%(n/10) == 0: print('%s completed' %((i+1)/float(n)))

    evaluate(alg)
    return alg

def replace1(alg, choices, p, seed):
    # now replace the first p "people", 1 at a time
    for i in range(p):
        alg.removeAt(i)
        alg.generate()
        for c in choices:
            alg.analyze(c)
        alg.update()

    evaluate(alg)

def test(choices, n, m, seed):
    alg = Algorithm(n, m, seed)
    alg.setArray(alg.original)
    evaluate(alg)
    return alg

def test_random(choices, n, m, seed):
    alg = Algorithm(n, m, seed)
    rnd = random.Random()
    pop = [choices[rnd.randint(0, 1023)] for i in range(n)]
    alg.setArray(pop)
    evaluate(alg)
    return alg

def main():
    n = 1000
    m = 10
##    seed = random.Random().randint(1,n)
    seed = 1000
    choices = list(it.product([0,1], repeat=m))
##    alg = test(choices, n, m, seed)
##    alg = test_random(choices, n, m, seed)
    alg = timed(lambda:run(choices, n, m, seed))
    replace1(alg, choices, 100, seed)

if __name__ == '__main__':
    os.chdir('d:/CMU/Summer.2011/dev/attr2')
    main()
