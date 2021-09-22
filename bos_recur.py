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
from printable import Basic as Printable,  B, Bd
from terms import Term, Terms, Annilator, Creator, Exp, t
import pickle

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

        if n == 0 and k == 0:#base case
            raise Exception("Base case of recursion not specified")

        elif n == 0 and k == 1:
            cls.Ks[key] = S(1).dot().simplify()

        elif n != 0 and k == 1:
            term1 = S(n+1).dot()
            term2 = S(n).lie_der(K(0,0))

            cls.Ks[key] = (term1 + term2).simplify()
            #this can be alternatively expressed as static of term2 - rotating
            #parts of other K(n,k'!=k)

        elif k > 1 and k <= n+1:
            # print("general")
            terms = Terms()
            for m in range(n):
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
        cls.Ks[str([0,0])] = H.simplify()

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

    @classmethod
    def get(cls, n):
        """
        get generator S(n) according to the recursive formula.
        """
        if n in cls.Ss.keys(): return cls.Ss[n]

        # print("S"+str(n))
        terms = Terms()
        if n == 0: return Terms()
        if n == 1:
            terms += smp.S(-1)*K(0,0)
        if n > 1:
            terms += smp.S(-1)*S(n-1).lie_der(K(0,0))
        for k in range(2, n+1):
            # print(k, n-1)
            terms += smp.S(-1)*K(n-1,k)
        terms = terms.rot().integrate()

        cls.Ss[n] = terms.simplify()
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


def K(n,k = -1):
    if k!= -1: return Kamiltonian.get(n,k)
    Kn = Terms()
    for ki in range(n+2):
        Kn += Kamiltonian.get(n,ki)
    return Kn

def S(n):
    return Generator.get(n)

# Term.set_lambda(smp.S(1))


omegad = Symbol('\omega_d')
omega0 = Symbol('\omega_a')
g3 = Symbol('g_3')
g4 = Symbol('g_4')
g5 = Symbol('g_5')
g6 = Symbol('g_6')
delta = Symbol('\delta')


a = Annilator('a')
ar = a*Exp(-1*I*omega0*t)
ad = Creator('a')
adr = ad*Exp(1*I*omega0*t)
pi = Symbol('\Pi ', complex=True)
pir = pi*Exp(-1*I*omegad*t)
pis = pi.conjugate()
pisr = pis*Exp(1*I*omegad*t)
H = delta*ad*a + g4*(ar+adr+pir+pisr)**4
H = g3*(ar+adr)**3+ g4*(ar+adr)**4+ g5*(ar+adr)**5#+ g6*(ar+adr)**6

Kamiltonian.set_H(H)

# import pickle
# with open('k5.txt','rb') as f:
#     k5 = pickle.load(f)

# with open('k4.txt','rb') as f:
#     k4 = pickle.load(f)

# with open('k3.txt','rb') as f:
#     k3 = pickle.load(f)

# with open('k2.txt','rb') as f:
#     k2 = pickle.load(f)

# with open('k1.txt','rb') as f:
#     k1 = pickle.load(f)

# with open('k0.txt','rb') as f:
#     k0 = pickle.load(f)

