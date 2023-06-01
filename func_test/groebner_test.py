import sympy as sy
import numpy as np
from funcs.func_utils import *
from funcs.hilbert import *
from funcs.groebner2 import *
from sympy.polys.monomials import Monomial
sy.init_printing(use_unicode=False, wrap_line=False)

# Test the function using example from reference (CARLO TRAVERSO)#
x = sy.symbols('x_0:4')
poly1 = sy.polys.Poly.from_dict({(2,1,0,0):1., (0,0,3,0): -1.}, x)
poly2 = sy.polys.Poly.from_dict({(1,3,0,0):1., (0,0,0,4): -1.}, x)
F = [poly1, poly2]
dom=poly1.domain
#W = sy.Matrix([[1,1,1,1], [0,0,0,1], [0,0,1,0], [0,1,0,0]])
#hilb, l = hilbert([(0,0,3,0), (0,0,0,4)], W, x)
W = sy.Matrix([[1,1,1,1]])
l = sy.symbols('l_0:1')
hilb = ((1-l[0]**3)*(1-l[0]**4))/((1-l[0])**4)
s = 1
U = W[:s, 0:]
d = [10]
groebner, mes = dtrunc_groebner(W, s, U, d, F, hilb, x, l, dom)
print(groebner)
print(mes)
