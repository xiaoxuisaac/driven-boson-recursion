from sympy import Basic, Symbol
from sympy.physics.secondquant import AnnihilateBoson, CreateBoson


class Printable(Basic):
    pass

class B(AnnihilateBoson):
    '''
    Hilbert space operator
    '''
    def _latex(self, printer):
        return "\\hat{%s}" % self.state.name

    def _dagger_(self):
        return Bd(*self.args)

    def conjugate(self):
        return self._dagger()


class Bd(CreateBoson):
    def _latex(self, printer):
        return "\\hat{%s}\\dagger " % self.state.name

    def _dagger_(self):
        return B(*self.args)
    def conjugate(self):
        return self._dagger()



class C(Symbol):
    '''
    phase space variable
    '''

    def __new__(cls, *args, **kwargs):
        # return super().__new__()
        return super().__new__(cls, *args, **kwargs)

    def _latex(self, printer):
        return "\\tilde{%s}" % self.name

    def _dagger_(self):
        return Cd(self.name)
    
    def conjugate(self):
        return self._dagger()



class Cd(Symbol):

    def __new__(cls, *args, **kwargs):
        # return super().__new__()
        if args[0][-1] == 'd':
            return super().__new__(cls, *args, **kwargs)
        return super().__new__(cls, args[0]+'d',  *args[1:], **kwargs)

    def _latex(self, printer):
        return "\\tilde{%s}^d" % self.name[:-1]

    def _dagger_(self):
        return C(self.name[:-1])

    def conjugate(self):
        return self._dagger()

