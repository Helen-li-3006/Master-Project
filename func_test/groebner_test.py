import sympy as sy
from sympy.polys.monomials import Monomial
import numpy as np
from funcs.func_utils import *
from funcs.hilbert import *
sy.init_printing(use_unicode=False, wrap_line=False)

x = sy.symbols('x_0:4')
poly1 = sy.poly('x_0**2 + x_1**3 + x_3**2', x)
y = sy.symbols('y')
xs = sy.symbols('x_0:2')
pol2 = poly1.subs({x[0]:xs[0], x[1]:xs[1], x[2]:1, x[3]:1})


poly3 = sy.poly('2*x_0**2 + 2*x_1**3', xs)

print(poly3 - pol2)