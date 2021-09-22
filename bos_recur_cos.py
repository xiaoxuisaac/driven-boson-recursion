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
            sta = sta.hermiticize()


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

omegad = Symbol('\omega_d', real = True)
# omega0 = Symbol('\omega_0', real = True)
g4 = Symbol('g_4', real = True)/4
# g3 = Symbol('g_3', real = True)
# g5 = Symbol('g_5', real = True)
g6 = Symbol('g_6', real = True)/6
# g7 = Symbol('g_7', real = True)
g8 = Symbol('g_8', real = True)/8
g10 = Symbol('g_{10}', real = True)/10
# delta = Symbol('\delta', real = True)
# Generator.gauge= 'A'
Term.set_phase_space()

a = Annilator('a')
ar = a*Exp(-I*3*omegad*t)
ad = Creator('a')
adr = ad*Exp(I*3*omegad*t)
pi = Symbol('\Pi ', complex=True)
pir = pi*Exp(-I*7*omegad*t)
pis = pi.conjugate()
pisr = pis*Exp(I*7*omegad*t)

# H1 = delta*ad*a + g3*(ar+adr+pir+pisr)**3
H2 = g4*(ar+adr+pir+pisr)**4
# H3 = g5*(ar+adr+pir+pisr)**5
H4 = g6*(ar+adr+pir+pisr)**6
# H5 = g7*(ar+adr+pir+pisr)**7
H6 = g8*(ar+adr+pir+pisr)**8
H8 = g10*(ar+adr+pir+pisr)**10

Kamiltonian.set_H([H2, H4, H6, H8])



# H1 = g3*(ar+adr)**3
# H2 = g4*(ar+adr)**4
# H3 = g5*(ar+adr)**5
# H4 = g6*(ar+adr)**6
# H5 = g7*(ar+adr)**7
# H6 = g8*(ar+adr)**8

# Kamiltonian.set_H([H1, H2, H3, H4, H5, H6])



#two mode
# b = Annilator('b')
# br = b*Exp(-I*omegad*t)
# bd = Creator('b')
# bdr = bd*Exp(I*omegad*t)
# Term.set_phase_space(['a','b'])


# H1 = g3*(ar+adr+br+bdr)**3
# H2 = g4*(ar+adr+br+bdr)**4
# Kamiltonian.set_H([H1, H2])

# # K(6)


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

