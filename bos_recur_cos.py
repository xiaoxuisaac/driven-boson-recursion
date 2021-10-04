#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 16 21:27:50 2021

@author: xiaoxu
"""


import sympy as smp
from sympy import Number, Add, Mul, Pow, Basic, latex, Symbol, I
from sympy.physics.secondquant import  Dagger
from copy import deepcopy
import pickle
from terms import Term, Terms, Annilator, Creator, Exp, t

class Kamiltonian(Terms):

    Ks = {}

    @classmethod
    def get(cls, n, k):
        """
        get Kamiltonian K(n)_[k] according to the recursive formula.
        """
        key = str([n,k])
        # print(key)

        if key in cls.Ks.keys():
            return cls.Ks[key]

        if n == 1 and k == 0:#base case
            raise Exception("Base case of recursion not specified")


        elif k == 1:
            term = S(n).dot()
            for m in range(1, n):
                term += S(n-m).lie_der(K(m,0))

            cls.Ks[key] = term.simplify()
            #this can be alternatively expressed as static of term2 - rotating
            #parts of other K(n,k'!=k)

        elif k > 1 and k <= n+1:
            # print("general")
            terms = Terms()
            for m in range(1,n):
                # print(m,k-1)
                terms += S(n-m).lie_der(K(m,k-1))/smp.S(k)
            cls.Ks[key] = terms.simplify()
        else:
            # print("zero case")
            return Terms()
        # cls.Ks[key].simplify()
        return cls.Ks[key]

    @classmethod
    def set_H(cls, H):
        for n, Hn in enumerate(H):
            cls.Ks[str([n+1,0])] = Hn.simplify()

    @classmethod
    def dump(cls, fdir):
        with open(fdir,'wb') as f:
            pickle.dump(cls.Ks,f)

    @classmethod
    def load(cls, fdir):
        with open(fdir,'rb') as f:
            cls.Ks = pickle.load(f)


class Generator(Terms):
    Ss = {}
    gauge = 0

    @classmethod
    def get(cls, n):
        """
        get generator S(n) according to the recursive formula.
        """
        if n in cls.Ss.keys(): return cls.Ss[n]

        # print("S"+str(n))
        osc = Terms()
        if n == 0: return Terms()

        osc += smp.S(-1)*K(n,0).rot()
        for m in range(1, n):
            osc += smp.S(-1)*S(n-m).lie_der(K(m,0)).rot()
        for k in range(2, n+1):
            # print(k, n-1)
            osc += smp.S(-1)*K(n,k).rot()
        osc = osc.integrate()
        cls.Ss[n] = osc.simplify()

        if cls.gauge != 0: # separate A into slow and fast
            sta = Terms()
            for k in range(2, n+1):
                sta += Eta.get(n,k).sta()
            sta = sta.lie_integrate(Term.modes[0])
            # sta = sta.hermiticize()


            cls.Ss[n] += sta.simplify()
        # print("S"+str(n)+"done")
        return cls.Ss[n]

    @classmethod
    def dump(cls, fdir):
        with open(fdir,'wb') as f:
            pickle.dump(cls.Ss,f)

    @classmethod
    def load(cls, fdir):
        with open(fdir,'rb') as f:
            cls.Ss = pickle.load(f)

class Eta(Terms):
    Etas = {}

    @classmethod
    def get(cls, n, k):
        key = str([n,k])
        if key in cls.Etas.keys():
            return cls.Etas[key]
        if k == 0:
            return Terms()
        if k == 1:
            #TODO: multimode
            cls.Etas[key] = S(n).lie_der(Term.modes[0])
            return cls.Etas[key]
        if k > 1 and k < n+1:
            terms = Terms()
            for m in range(1, n):
                terms += S(n-m).lie_der(cls.get(m, k-1))
            cls.Etas[key] = terms
            return cls.Etas[key]
        return Terms()

def K(n,k = -1):
    if k!= -1: return Kamiltonian.get(n,k)
    Kn = Terms()
    for ki in range(n+2):
        Kn += Kamiltonian.get(n,ki)
    return Kn.simplify()

def S(n):
    return Generator.get(n)

