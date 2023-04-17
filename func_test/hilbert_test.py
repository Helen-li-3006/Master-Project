import sympy as sy
from sympy.polys.monomials import Monomial
import numpy as np
from funcs.func_utils import *
from funcs.hilbert import *
sy.init_printing(use_unicode=False, wrap_line=False)

# Replicate of example from a reference
x = sy.symbols('x_0:4')
Ims = [(0,0,3,0), (0,0,0,4)]
# Example uses Lex ordering
W = sy.Matrix([[1,1,1,1]])

hilb, l = hilbert(Ims, W, x)
hilbnum = hilbinumerator(Ims, W, x, sy.symbols('l'))
print('Numerator of the Hilbert series is:', hilbnum)
print('Full Hilbert series is ', hilb)
print('Expanded series up to order 5 is', hilb_expand(hilb, l, [5]))

