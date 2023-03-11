import sympy as sy
from sympy.polys.monomials import Monomial
import numpy as np
from funcs.func_utils import *
from funcs.hilbert import *
sy.init_printing(use_unicode=False, wrap_line=False)

# Replicate example 1.2.17 from Gatermann for hibert series
x = sy.symbols('x_1:4')
Ims = [(2,0,2), (0,2,2), (1,0,3), (2,1,0), (1,2,0), (1,1,1)]
# Example uses natural grading
W = sy.Matrix([[1,1,1]])

hilb = hilbert(Ims, W, x)
print(sy.simplify(hilb))

