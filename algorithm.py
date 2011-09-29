#-------------------------------------------------------------------------------
# Name:        algorithm
# Purpose:
#
# Author:      jackchen
#
# Created:     05/08/2011
# Copyright:   (c) jackchen 2011
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#!/usr/bin/env python

import itertools as it
import functools as ft
import math
import numpy as np
import os
import random
import time
from collections import defaultdict

class Algorithm:
    """creates an instance of the algorithm"""
    def __init__(self, n, m, seed):
        self.n = n
        self.m = m
        self.seed = seed
        self.original, marg, intr, self.mlist, self.ilist = Algorithm.build(self.n, self.m, self.seed)
        self.expected = Bunch(marginals=marg, interactions=intr)
        self.a = []
        self.counter = 0
        self.actual = Bunch(marginals=np.array([0]*len(marg)), \
                            interactions=dict((k,np.array([0]*4)) for k in self.ilist))
        self.mlookup = [0]*self.m
        self.ilookup = np.repeat([np.array([0]*4)], len(self.ilist), axis=0)

    def generate(self):
        """generate the calculated marginals and interactions based on the possible
        incremental updates for easy lookup during the analysis phase"""
        # marginal generation functions
        def fe(i):
            a = self.expected.marginals[i]
            return lambda x:Algorithm.__wm_weight__*error(x, a)
        f = lambda i: (self.actual.marginals[i]/float(self.counter+1),
                       (self.actual.marginals[i]+1)/float(self.counter+1),
                       fe(i))
        # interaction generation functions
        def ge(i,j):
            a = self.expected.interactions[i,j]
            return lambda x:Algorithm.__wi_weight__*error(Algorithm.__interaction__(x), a)
        g = lambda i,j: (self.actual.interactions[i,j] + Algorithm.__c0__,
                         self.actual.interactions[i,j] + Algorithm.__c1__,
                         self.actual.interactions[i,j] + Algorithm.__c2__,
                         self.actual.interactions[i,j] + Algorithm.__c3__,
                         ge(i,j))

        self.mlookup = [(fn(a),fn(b)) for a,b,fn in (f(i) for i in self.mlist)]
        self.ilookup = [(gn(a),gn(b),gn(c),gn(d)) for a,b,c,d,gn in
                        (g(i,j) for i,j in self.ilist)]
        self.minimum = (-1,0)

    def analyze(self, choice):
        """runs an analysis (evaluaton) for a specific population choice using the
        generated lookup values. Once a minimum is identified, either replace the previous
        minimum, or keep a collection of min values to select from later"""
        # lookup on marginals
        add = lambda x,y: x+y
        mresult = ft.reduce(add, (
            m[c] for c,m in zip(choice, self.mlookup)))
        # lookup on interactions
        iresult = ft.reduce(add, (
            m[Algorithm.__match__[x,y]] for (x,y),m in
                zip(((choice[i], choice[j]) for i,j in self.ilist), self.ilookup)))

        result = mresult + iresult

        # update the minimum
        if self.minimum[0] == -1: self.minimum = (result, choice)
        elif result < self.minimum[0]: self.minimum = (result, choice)
        # TODO: track "multiple" most fit, for now, first one wins

    def update(self):
        """completes the iteration and saves the value of the most 'fit' person"""
        choice = self.minimum[1]
        self.a += [choice]
        self.counter += 1
        self.actual.marginals += choice
        for i,j in self.actual.interactions:
            self.actual.interactions[i,j] += Algorithm.__convert__[choice[i], choice[j]]
        return choice

    def removeAt(self, i):
        """removes the choice from the array and updates the actual marginals"""
        choice = self.a[i]
        del self.a[i]
        self.counter -= 1
        self.actual.marginals -= choice
        for i,j in self.actual.interactions:
            self.actual.interactions[i,j] -= Algorithm.__convert__[choice[i], choice[j]]

    def setArray(self, a):
        """primes the array to some preset value"""
        self.a = a
        self.counter = len(a)
        temp = np.array(a)
        self.actual.marginals = np.sum(temp, axis=0)
        for i,j in self.ilist:
            _,_,self.actual.interactions[(i,j)] = \
                Algorithm.__extrapolate__(temp, i, j, self.counter, raw=1)

    # Static properties
    __c0__ = [1,0,0,0]
    __c1__ = [0,1,0,0]
    __c2__ = [0,0,1,0]
    __c3__ = [0,0,0,1]
    __convert__ = {(0,0):__c0__,
                   (0,1):__c1__,
                   (1,0):__c2__,
                   (1,1):__c3__}
    __match__ = {(0,0):0,
                 (0,1):1,
                 (1,0):2,
                 (1,1):3}
    __wm_weight__ = 4
    __wi_weight__ = 1

    # Static functions
    @staticmethod
    def __marginal__(name, rnd, pr=-1):
        """returns a marginal with a pre-defined probability"""
        pr = pr if pr > 0 else rnd.uniform(0,1)
        return (name, pr)

    @staticmethod
    def __interaction__(cells):
        """defines an pairwise interaction between attributes i and j"""
        log = math.log
        f = lambda x: x if x>0 else 0.001
        value = log(f(cells[0]))+log(f(cells[3]))-log(f(cells[1]))-log(f(cells[2]))
        value = logit(value)
        return value

    @staticmethod
    def __extrapolate__(a, i, j, n=-1, raw=0):
        """extrapolates the interaction value of attributes i and j"""
        n = n if n>0 else len(a)
        # first pull i,j out of a
        ai = a[:,i]
        aj = a[:,j]
        # zip to create tuple values of (ci, cj) for c in (1..n)
        c = zip(ai,aj)
        # translate c into a contingency table
        vals = [Algorithm.__convert__[(ci,cj)] for ci,cj in c]
        result = np.sum(vals, axis=0)
        # TODO: use p_hat [ni/float(n) for ni in result] or just ni?
        if raw: return i, j, result
        return Algorithm.__interaction__(cells=result)

    @staticmethod
    def build(n, m, seed):
        """builds the initial population marginals and interactions"""
        rnd = random.Random(seed)
        mlist = range(m)
        ilist = list(it.combinations(range(m),2))
        marginals = [Algorithm.__marginal__(i, rnd) for i in mlist]
        f = lambda pr: 1 if pr >= rnd.uniform(0,1) else 0
        pop = [[f(pr) for _,pr in marginals] for x in range(n)]
        a = np.array(pop)
        marginals = dict([(name,sum(a[:,i])/float(n))
                          for i,(name,_) in enumerate(marginals)])
        interactions = dict([((i,j), Algorithm.__extrapolate__(a,i,j,n)) for i,j in ilist])
        return pop, marginals, interactions, mlist, ilist

    @staticmethod
    def evaluate(a, marginals, interactions, n=-1, out=0):
        """evaluates the fitness of the array"""
        n = n if n>0 else len(a)
        # first evaluate the marginals
        f = lambda i,pr: sum(a[:,i])/float(n)
        fw = lambda x: Algorithm.__wm_weight__*x
        # next evaluate the interactions
        g = lambda idx,i,j: Algorithm.__extrapolate__(a,i,j,n)
        gw = lambda x: Algorithm.__wi_weight__*x

        if out:
            print ("name\twt\ttarget\t\tactual\t\tweighted diff^2\t")
            print ('-'*50)
            results = [(name,Algorithm.__wm_weight__,pr,f(i,pr),fw(error(pr,f(i,pr)))) \
                       for i,(name,pr) in enumerate(marginals.items())]
            for e in results[:len(marginals)]: print ('attr%s\t%s\t%f\t%f\t%f\t' %(e))
            print ('-'*50)
            results += [('i_%s%s' %(i,j),Algorithm.__wi_weight__,val,g(idx,i,j),gw(error(val,g(idx,i,j)))) \
                        for idx,((i,j),val) in enumerate(interactions.items())]
            print ("name\twt\ttarget\t\tactual\t\tweighted diff^2\t")
            print ('-'*50)
            for e in results[len(marginals):]: print ('%s\t%s\t%f\t%f\t%f\t' %(e))
            print ("sum of square errors:", sum([e[-1] for e in results]))
        else:
            results = [fw(error(pr,f(i,pr))) for i,(name,pr) in enumerate(marginals)]
            results += [gw(error(val,g(idx,i,j))) for idx,((i,j),val) in enumerate(interactions)]
##        results = [fw(error(pr,f(i,pr))) for i,(name,pr) in enumerate(marginals.items())]
##        results += [gw(error(val,g(idx,i,j))) for idx,((i,j),val) in enumerate(interactions.items())]
##        err = sum(results)
##        print ("sum of square errors:", err)

        err = sum([e[-1] for e in results]) if out else sum(results)
        return err,a

class Bunch:
    def __init__(self, **kwds):
        self.__dict__.update(kwds)

def logit(x):
    """Logit function"""
    return 1.0 / (1.0 + np.exp(-x))

def error(x, y):
    """Error function"""
    return (x-y)**2

def timed(fun):
    start = time.time()
    try:
        return fun()
    finally:
        elapsed = time.time() - start
        print ('-'*50)
        print ("time elapsed:",elapsed)

def main():
    alg = Algorithm(1000, 10, 1000)
##    for m in alg.expected.marginals: print(m)
##    for e in alg.ilist: print(e)
    for e in alg.actual.interactions.items(): print(e)

if __name__ == '__main__':
    main()
