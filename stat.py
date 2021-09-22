#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 16 22:02:24 2021

@author: xiaoxu
"""
import sympy as smp
from sympy import Number, Add, Mul, Pow, Basic, latex, Symbol, I
from sympy.physics.secondquant import  Dagger
from copy import deepcopy
import pickle
from terms import Term, Terms, Annilator, Creator, Exp, t
from bos_recur_cos import K, Kamiltonian


def test():
    omegad = Symbol('\omega_d', real = True)
    omega0 = Symbol('\omega_0', real = True)
    g4 = Symbol('g_4', real = True)
    delta = Symbol('\delta', real = True)
    # Generator.gauge= 'A'
    Term.set_phase_space()

    a = Annilator('a')
    ar = a*Exp(-I*omega0*t)
    ad = Creator('a')
    adr = ad*Exp(I*omega0*t)
    pi = Symbol('\Pi ', complex=True)
    pir = pi*Exp(-I*omegad*t)
    pis = pi.conjugate()
    pisr = pis*Exp(I*omegad*t)
    H = delta*ad*a + g4*(ar+adr+pir+pisr)**4

    Kamiltonian.set_H([H])
    K(3)


import cProfile, pstats
profiler = cProfile.Profile()
profiler.enable()
test()
profiler.disable()
stats = pstats.Stats(profiler)
stats.strip_dirs()
stats.sort_stats('cumtime')
stats.print_stats(50)