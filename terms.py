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
from no import normal
from printable import Basic as Printable,  B, Bd, C, Cd
import pickle

hbar = Symbol("hbar")
t = Symbol('t')



class Term(Printable): #in Hilbert space: normal-ordered. or in Phase space, Moyal product
    is_Term = True
    is_Terms = False
    lam = Symbol("hbar")
    space = "Hilbert"
    star_map = {}
    lie_der_map = {}

    def __mul__(self, other):
        if self.space == "Hilbert":
            return self.__mul__Hilbert(other)
        return self.__mul__phase(other)

    def __rmul__(self, other):
        if self.space == "Hilbert":
            return self.__rmul__Hilbert(other)
        return self.__rmul__phase(other)

    def lie_der(self, term):
        """
        compute L_{self} term

        Parameters
        ----------
        term : Term
            DESCRIPTION.
        Returns
        -------
        Terms
        """
        if self.space == "Hilbert":
            return self.lie_der_Hilbert(term)
        return self.lie_der_phase(term)

    def dot(self):
        freq = self.frequency
        if freq == 0: return Terms()

        if freq.is_Mul:
            factor = smp.S(1)
            omegas = smp.S(1)
            for m in freq.args:
                if m.is_Number:
                    factor *= m
                else:
                    omegas *= m
            frac = Prefactor((I*factor, smp.S(1)/omegas))
        else:
            frac = Prefactor((I, smp.S(1)/freq))

        return Term(self.args[0]*frac, *self.args[1:])
        return self*(I*freq)



    def integrate(self, var = 't'):
        if var == 't':
            return self._integrate_t()

    def _integrate_t(self):
        freq = self.frequency
        if freq == 0:
            raise Exception("secular term generated from time integration")

        if freq.is_Mul:
            factor = smp.S(1)
            omegas = smp.S(1)
            for m in freq.args:
                if m.is_Number:
                    factor *= m
                else:
                    omegas *= m
            frac = Prefactor((smp.S(1)/factor/I, omegas))
        else:
            frac = Prefactor((smp.S(1)/I, freq))

        return Term(self.args[0]*frac, *self.args[1:])

    def lie_integrate(self, mode):
        if self.space == "Hilbert":
            return self.lie_integrate_Hilbert(mode)
        return self.lie_integrate_phase(mode)

    @property
    def frequency(self):
        if self.args[2] == 1: return smp.S(0)
        if isinstance(self.args[2], smp.exp):
            return self.args[2].args[0]/I/t
        if self.args[2].is_Mul:
            freq = smp.S(0)
            for freq_compo in self.args[2].args:
                freq += freq_compo.args[0]
            return (freq/I/t).expand()

    def _is_same_type(self, other):
        # check if two terms contain same type of operator and oscillating part
        return self.args[1] == other.args[1] and self.args[2] == other.args[2]


    def __add__(self, other):
        if other.is_Term:
            if self._is_same_type(other):
                return Term(self.args[0]+other.args[0], *self.args[1:])
            return Terms(self, other)

        elif isinstance(other, Terms):
            return other + self
        else:
            raise NotImplementedError

    def __sub__(self, other):
        return self + -1*other

    def __pow__(self, power):
        result = self
        for i in range(power-1):
            result *= self
        return result

    def _latex(self, printer):
        e = ''
        if self.args[0].expression != 1:
                e += latex(self.args[0])
        e += latex(self.args[1]) if self.args[1] != 1 else ''
        e += latex(self.args[2]) if self.args[2] != 1 else ''
        if e == '': return '1'
        return e

    def __str__(self):
        e = self.args[0]*self.args[1]
        return str(e) if self.args[-1] == 1 else str(e) + str(self.args[-1])

    def as_coeff_Mul(self, *deps, **kwargs):
        if deps:
            if not self.has(*deps):
                return self, tuple()
        return smp.S.One, (self,)

    def __truediv__(self, scalar):
        '''
        TODO: what if in phase space?
        '''
        if isinstance(scalar, (float, int)) or scalar.is_commutative:
            return Term(self.args[0]/scalar, *self.args[1:])
        raise Exception("illegal division, " + str(scalar))


    def select(self, term):
        if term.is_Terms:
            if len(term) == 1:
                term = term.args[0]
            else:
                raise Exception("too many terms")
        if self.args[1] != term.args[1]: return Terms()
        #TODO: prefactor
        if term.args[0].expression.is_Number:
            return self

        variables = []
        if term.args[0].expression.is_Mul:
            for ti in term.args[0].expression.args:
                if ti.is_Number:
                    pass
                elif ti.is_Pow:
                    variables.append(ti.args[0])
                else:
                    variables.append(ti)
        else:
            ti = term.args[0].expression
            if ti.is_Pow:
                    variables.append(ti.args[0])
            else:
                variables.append(ti)

        if self.args[0].is_Add:
            result = []
            for factor in self.args[0].args:
                test = factor / term.args[0]
                select = True
                for v in variables:
                    if test.has(v): select = False
                if select: result.append(factor)
            if len(result) == 0: return Terms()
            return Term(Add(*result), *self.args[1:])
        else:
            test = self.args[0] / term.args[0]
            select = True
            for v in variables:
                if test.has(v): return Terms()
            return self

    @classmethod
    def set_lambda(cls, lam):
        cls.lam = lam

    def evalf(self, *args):
        return Term(self.args[0].expression.evalf(*args), *self.args[1:])

    def expand(self):
        return Term(self.args[0].expand(), *self.args[1:])

    #### Hilbert-space-specific methods
    @property
    def is_commutative(self):
        return self.args[1].is_commutative

    def __mul__Hilbert(self, other): #multiply and then normal order the product
        if not isinstance(other, Basic):
            # if the multiplicant is just a number
            return Term(self.args[0]*other, *self.args[1:])


        if not isinstance(other, (Terms,Term)) and isinstance(other, Basic):
            # if the multiplicand is a sympy object.
            # assumed to be normal ordered
            if other.is_commutative:
                return Term(self.args[0]*other, *self.args[1:])
            if not other.is_Mul:
                # if other is a simple non-commutative operator, self is assumed
                # to have 1 in args[1]
                return self*Term(Prefactor((smp.S(1),smp.S(1))), other, smp.S(1))
            else:
                c_part = []
                string1 = []
                for factor in other.args:
                    if factor.is_commutative:
                        c_part.append(factor)
                    else:
                        string1.append(factor)

                if self.is_commutative:
                   return Term(self.args[0]*Mul(*c_part), Mul(*string1),
                               *self.args[2:])

                return self*Term(Mul(*c_part), Mul(*string1), smp.S(1))

        if other.is_Terms: return other.__rmul__(self)

        if self.is_commutative or other.is_commutative:
            # if one of the multiplicant is commutative
            return Term(self.args[0]*other.args[0], self.args[1]*other.args[1],
                               self.args[2]*other.args[2])

        #nontrivial case, normal ordering
        factor_term = Term(self.args[0]*other.args[0], smp.S(1),
                           self.args[2]*other.args[2])

        no_terms = normal(self.args[1]*other.args[1], lam = self.lam)
        no_terms = no_terms.args if no_terms.is_Add else [no_terms]

        new_terms = []
        # construct list directly because all no_terms are unique
        for no_term in no_terms:
            new_terms.append(factor_term*no_term)
        return Terms(*new_terms)

    def __rmul__Hilbert(self, other):
        if isinstance(other, (int, float)):
            # if the multiplicanc is just a number
            return Term(self.args[0]*other, *self.args[1:])

        if isinstance(other, Basic) and other.is_commutative:
            return Term(self.args[0]*other, *self.args[1:])

        if isinstance(other, Basic) and not isinstance(other, Term):
            raise Exception("right multiplication with sympy obj \
                            shouldn't happend")

        return other.__mul__(self)


    def lie_der_Hilbert(self, term):
        """
        compute [self, term]/i\hbar

        Parameters
        ----------
        term : Term
            DESCRIPTION.
        Returns
        -------
        Terms
        """
        if self.is_commutative or term.is_commutative:
            return Term(Prefactor(), smp.S(1), smp.S(1))
        if self.args[1] == term.args[1]:
            return Term(Prefactor(), smp.S(1), smp.S(1))
        return (self*term - term*self)/(I*hbar)

    def lie_integrate_Hilbert(self, mode):
        raise NotImplementedError


    #### phase-space-specific methods

    @classmethod
    def set_phase_space(cls, modes = ['a'], lam = Symbol("\hbar")):
        cls.space = "phase"
        cls.modes = [Annilator(mode) for mode in modes]


    @classmethod
    def _star_prod(cls, expr1, expr2):
        '''exp1 and exp2 are assumed to be 1 or product of mode coordinates
        '''
        # pow1 = cls._get_power(expr1)
        # pow2 = cls._get_power(expr2)

        # c = cls.modes[0].args[1]
        # cd = cls.modes[0].args[1]._dagger_()

        # if (pow1, pow2) in cls.star_map.keys():
        if (expr1, expr2) in cls.star_map.keys():
            return cls.star_map[(expr1, expr2)]

        pb_flag = True
        counter = 0
        pbs = smp.S(0)

        modes = [mode.args[1] for mode in cls.modes]
        modes += [mode.args[1]._dagger_() for mode in cls.modes]


        while(pb_flag):
                pb_flag = False
            # for i in range(counter+1):

                pb, pb_flag = cls.diff_multinomial(expr1, expr2, counter, modes)

                # pb = expr1.diff(c, i).diff(c._dagger_(), counter -i)\
                    # *expr2.diff(c._dagger_(), i).diff(c, counter -i)
                # pb*=(-1)**(counter-i)/smp.factorial(i)/smp.factorial(counter-i)
                pb*= (cls.lam/2)**counter
                pbs += pb
                # pb_flag = pb_flag or (pb!=0)
                counter += 1


        # while(pb_flag):
        #     pb_flag = False
        #     for i in range(counter+1):
        #         pb = expr1.diff(c, i).diff(c._dagger_(), counter -i)\
        #             *expr2.diff(c._dagger_(), i).diff(c, counter -i)
        #         pb*=(-1)**(counter-i)/smp.factorial(i)/smp.factorial(counter-i)
        #         pb*= (cls.lam/2)**counter
        #         pbs += pb
        #         pb_flag = pb_flag or (pb!=0)
        #     counter += 1

        pbs = pbs.expand()

        if isinstance(pbs, Add):
            pbs = list(pbs.args)
        else:
            pbs = [pbs]

        terms = Terms()
        f = Term(Prefactor((smp.S(1),smp.S(1))),  smp.S(1), smp.S(1))
        for pb in pbs:
            terms += f._mul_pb(pb)
        # pbs = [f._mul_pb(pb) for pb in pbs]


        cls.star_map[(expr1, expr2)] = terms
        # cls.star_map[(pow2, pow1)] = [-1*pb for pb in pbs]
        return terms

    @classmethod
    def diff_multinomial(cls, expr1, expr2, n, modes):
        ''' return expr1 ({<-, ->})^n expr2, flag
        flag is a boolean variable labeling whether n is still with in the range
        such that NOT all generated summands are zero
        '''
        if n == 0: return expr1 * expr2, True
        if len(modes) == 1:
            expr1 = expr1.diff(modes[0], n)
            expr2 = expr2.diff(modes[0]._dagger_(), n)
            r = expr1*expr2
            if isinstance(modes[0], Cd):
                r = r*smp.S(-1)**n
            r = r/smp.factorial(n)
            return r, r!=0

        rs = smp.S(0)
        flag = False
        for i in range(n+1):
            mode = modes[0]
            new_expr1 = expr1.diff(mode, i)
            new_expr2 = expr2.diff(mode._dagger_(), i)
            if new_expr1 == 0 or new_expr2 == 0:
                #if at this order already zeor, higher order in differentiation has
                # to be zero.
                break
            r = smp.S(0)
            if expr1 !=0 and expr2 != 0:
                r, flgi = cls.diff_multinomial(new_expr1, new_expr2, n-i, modes[1:])
            if isinstance(mode, Cd):
                r *= (-1)**i
            r = r/smp.factorial(i)
            rs += r
            flag = flag or (r != 0)
        return rs, flag


    def __mul__phase(self, other):
        #TODO: multi-mode
        #TODO: optimize
        # if is Term, implement Greonewold product
        if isinstance(other, Term):
            if self.args[1] == 1 or other.args[1] == 1:
                return Term(self.args[0]*other.args[0], self.args[1]*other.args[1], self.args[2]*other.args[2])

            f = Term(self.args[0]*other.args[0],  smp.S(1), self.args[2]*other.args[2])
            pbs = self._star_prod(self.args[1], other.args[1])
            return pbs*f

        if isinstance(other, Terms): return other.__rmul__(self)

        # if the multiplicanc is just a number
        return Term(self.args[0]*other, *self.args[1:])


    def __rmul__phase(self, other):
        if isinstance(other, (int, float)):
            # if the multiplicanc is just a number
            return Term(self.args[0]*other, *self.args[1:])

        if isinstance(other, Basic) and not isinstance(other, Term):
            return Term(self.args[0]*other, *self.args[1:])

        return other.__mul__(self)

    def _mul_pb(self, pb):
        # multiple by a sympy expression from Poisson Bracket.
        # pb should only contain numbers and annilator/creator.
        # Lemma: any term from PB is not an Add instance.
        if pb == 0:
            return Term(Prefactor(), smp.S(1), smp.S(1))

        factor = smp.S(1)
        sym = smp.S(1)
        if pb.is_Mul:
            for arg in pb.args:
                if arg.is_Number or arg.has(self.lam):
                    factor *= arg
                else:
                    sym *= arg
            return Term(self.args[0]*factor, self.args[1]*sym, *self.args[2:])

        if pb.is_Number or pb.has(self.lam):
            return Term(self.args[0]*pb, smp.S(1), *self.args[2:])

        return Term(self.args[0], self.args[1]*pb, *self.args[2:])


    def lie_der_phase(self, term):
        """
        compute {{self, term}}/factor

        Parameters
        ----------
        term : Term
            DESCRIPTION.
        factor : TYPE, optional
            DESCRIPTION. The default is smp.S(1).

        Returns
        -------
        Terms
        """
        if self.args[1] == term.args[1]:
            return Term(Prefactor(), smp.S(1), smp.S(1))
        f =  Term(self.args[0]*term.args[0]/(I*hbar),  smp.S(1),
                  self.args[2]*term.args[2])
        if (self.args[1], term.args[1]) in self.lie_der_map.keys():
            mb = self.lie_der_map[(self.args[1], term.args[1])]
        else:
            mb = self._star_prod(self.args[1],term.args[1])\
                            - self._star_prod(term.args[1], self.args[1])
            self.lie_der_map[(self.args[1], term.args[1])] = mb
            self.lie_der_map[(term.args[1], self.args[1])] = -1*mb

        return f*mb

    def lie_integrate_phase(self, mode):
        # the factor I takes care of that in the moyal bracket
        f = Term(self.args[0]*I, smp.S(1), self.args[2])
        integ = self.args[1].integrate(mode.args[1]._dagger_())
        if isinstance(mode.args[1], Cd):
            integ *= smp.S(-1)
        return f._mul_pb(integ)

    def conjugate(self):
        return Term(self.args[0].conjugate(), Dagger(self.args[1]), self.args[2].conjugate())



    def weyl_transform(self):
        if self.args[1] == 1:
            return self

        terms = Terms(self)
        for mode in self.modes:
            c = mode.args[1]
            new_terms = Terms()
            for term in terms.args:
                new_terms += term.weyl_transform_mode(c)
            terms = new_terms
        return terms



    def weyl_transform_mode(self, c):
        #Weyl quantization to normal order
        if self.args[1] == 1:
            return self

        terms = Terms()
        f = Term(self.args[0],  smp.S(1),  *self.args[2:])
        term = self.args[1]
        j = smp.S(1)

        while(True):
            if term == 0:
                break
            terms += f._mul_pb(term)

            term *= self.lam/smp.S(2)/j
            term = smp.diff(term, c)
            term = smp.diff(term, c._dagger_())
            j += 1

        return terms




class Prefactor(Printable):
    '''prefactor of each term. for optimization.
    '''
    @property
    def factors(self):
        if '_factors' in self.__dict__.keys():
            return self._factors
        self._factors = {}
        for num, denom in self.args:
            self.factors[denom] = num
        return self._factors

    @property
    def expression(self):
        if '_expression' in self.__dict__.keys():
            return self._expression
        expr = smp.S(0)
        for  num, denom in self.args:
            expr += num/denom
        self._expression = expr
        return self._expression

    def _latex(self, printer):
        ltx = latex(self.expression)
        if len(self.args)>1 or (self.args != () and self.args[0][0].is_Add and self.args[0][1]==1):
            ltx = '('+ltx+')'
        return ltx

    def __str__(self):
        s = str(self.expression)
        if len(self.factors)>1 or (self.args != () and self.args[0][0].is_Add and self.args[0][1]==1):
            s = '('+s+')'
        return s

    def __mul__frac(self, frac=(smp.S(1), smp.S(1))):
        args = []
        for num, denom in self.args:
            if frac[0] != 1: num *= frac[0]
            if frac[1] != 1: denom *= frac[1]
            args.append((num, denom))
        return Prefactor(*args)

    @classmethod
    def from_dict(cls, factors):
        #TODO: combine common key that needs to be expanded.
        args = []
        for denom, num in factors.items():
            if num != 0:
                args.append((num, denom))
        p = cls.__new__(cls, *args)
        p._factors = factors
        return p

    def __mul__(self, other):
        if isinstance(other, Prefactor):
            #ontrivial case
            factors = {}
            for num1, denom1 in self.args:
                for num2, denom2 in other.args:
                    denom = denom1*denom2
                    num = num1*num2
                    if denom in factors.keys():
                        factors[denom]+=num
                    elif -1*denom in factors.keys():
                        factors[denom]-=num
                    else:
                        factors[denom] = num
            return self.from_dict(factors)

        else:
            return self.__mul__frac((other, smp.S(1)))

    def __rmul__(self, other):
        return self.__mul__(other)


    def __add__(self, other):
        if not isinstance(other, self.__class__):
            raise Exception("Cannot add a Prefactor with anoth non-Prefactor")
        factors = self.factors.copy()
        for num, denom in other.args:
            if denom in factors.keys():
                factors[denom]+=num
            elif -1*denom in factors.keys():
                factors[-1*denom]-=num
            else:
                factors[denom] = num
        return self.from_dict(factors)


    def __truediv__(self, other):
        return self.__mul__frac((smp.S(1)/other, smp.S(1)))

    def conjugate(self):
        args = [(arg[0].conjugate(), arg[1]) for arg in self.args]
        return Prefactor(*args)
    @property
    def is_zero(self):
        return len(self.args) == 0

    def expand(self):
        args = []
        for arg in self.args:
            args.append((arg[0].expand(), arg[1]))
        return Prefactor(*args)

class Annilator(Term):
    def __new__(cls, *args):
        # return super().__new__()
        if cls.space == 'Hilbert':
            return super().__new__(cls, Prefactor((smp.S(1), smp.S(1))), B(args[0]), smp.S(1))
        else:
            return super().__new__(cls, Prefactor((smp.S(1), smp.S(1))), C(args[0]), smp.S(1))


class Creator(Term):
    def __new__(cls, *args):
        if cls.space == 'Hilbert':
            return super().__new__(cls, Prefactor((smp.S(1), smp.S(1))), Bd(args[0]), smp.S(1))
        else:
            return super().__new__(cls, Prefactor((smp.S(1), smp.S(1))), Cd(args[0]), smp.S(1))


class Exp(Term):
    def __new__(cls, *args):
        # return super().__new__()
        return super().__new__(cls, Prefactor((smp.S(1), smp.S(1))), smp.S(1), smp.exp(args[0]))



class Terms(Printable):
    is_Terms = True
    is_Term = False

    def lie_der(self, terms):
        new_terms = Terms()
        if isinstance(terms, Term):
            terms = Terms(terms)
        for ti in self.args:
            for tj in terms.args:
                new_terms += ti.lie_der(tj)
        return new_terms

    def lie_integrate(self, mode):
        new_terms = Terms()
        for ti in self.args:
            new_terms += ti.lie_integrate(mode)
        return new_terms


    def dot(self):
        new_terms = []
        for ti in self.args:
            tidot = ti.dot()
            if tidot.args != (): new_terms.append(ti.dot())
        return Terms(*new_terms)

    def rot(self):
        new_terms = []
        for ti in self.args:
            if ti.args[2] != 1: new_terms.append(ti)
        return Terms(*new_terms)

    def sta(self):
        new_terms = []
        for ti in self.args:
            if ti.args[2] == 1: new_terms.append(ti)
        return Terms(*new_terms)

    def integrate(self, var = 't'):
        new_terms = []
        for ti in self.args:
            new_terms.append(ti.integrate(var))
        return Terms(*new_terms)



    def __add__(self,  terms):
        #TODO: if terms is just an ordinary sympy object
        if terms.is_Term:
            # base case. add through each term
            for i, term in enumerate(self.args):
                if term._is_same_type(terms):
                    term_i = term + terms
                    term_i = [term_i] if term_i.args[0] != 0 else []
                    if i == len(self.args) - 1:
                        return Terms(*self.args[:i], *term_i)
                    return Terms(*self.args[:i], *term_i, *self.args[i+1:])
            return Terms(*self.args, terms)

        if self.args == (): return terms

        result = self.copy()
        for term in terms.args:
            result += term
        return result

    def _latex(self, printer):
        e = ''
        for i, ei in enumerate(self.args):
            ei = latex(ei)
            e += '+' if ei[0] not in '+-' and i != 0 else ''
            e += ei
        return e

    def __str__(self):
        e = smp.S(0)
        for ei in self.args:
            e += Mul(ei.args[0].expression, *ei.args[1:])
        return str(e)

    def __len__(self):
        return len(self.args)


    def __mul__(self, other):
        #TODO if other is an regular sympy object. seems working
        if not isinstance(other, Terms):
            args = [other]
        else:
            args = other.args
        new_terms = Terms()
        for term1 in self.args:
            for term2 in args:
                new_terms += term1*term2
        return new_terms

    def __rmul__(self, other):
        #TODO if other is an regular sympy object. seems working
        if not isinstance(other, Terms):
            other = Terms(other)
        return other.__mul__(self)

    def __pow__(self, power):
        result = self
        for i in range(power-1):
            result *= self
        return result

    def __truediv__(self, scalar):
        '''
        TODO: what if Terms is in phase space?
        '''
        if isinstance(scalar, (float, int)) or scalar.is_commutative:
            terms = []
            for term in self.args:
                terms.append(term/scalar)
            return Terms(*terms)
        raise Exception("illegal division, " + str(scalar))

    def __sub__(self, other):
        return self + -1*other

    def simplify(self):
        # return self
        #TODO: prefactor
        terms = []
        for term in self.args:
            if (not term.args[0].is_zero) and term.args[1] != 1:
                terms.append(term.expand())
        return Terms(*terms)

    def expand(self):
        return self
        #TODO
        terms = []
        for term in self.args:
            term = Term(term.args[0].expand(), *term.args[1:])
            terms.append(term)
        return Terms(*terms)

    def evalf(self, *args):
        terms = self.simplify()
        result = []
        for term in terms.args:
            result.append(term.evalf(*args))
        return Terms(*result)

    def select(self, term):
        result = Terms()
        for ti in self.args:
            result += ti.select(term)
        return result

    #### phase-space-specific methods
    def weyl_transform(self):
        results = Terms()
        for term in self.args:
            results += term.weyl_transform()
        return results.simplify()

    def hermiticize(self):
        extr = Terms()
        for term in self.args:
            term_c = term.conjugate()
            if term_c != term:
                contain = [term_c == termi for termi in self.args]
                if not any(contain):
                    extr += term_c
        # print(self)
        # print(extr)
        return self+extr








