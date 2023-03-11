import sympy as sy
import numpy as np
from funcs.func_utils import *
sy.init_printing(use_unicode=False, wrap_line=False)

# Replicate the examples from fig 1.2 from Gatermann
W1 = sy.Matrix([[1,1],[1,2]])
W2 = sy.Matrix([[3,1],[1,-1]])
x = sy.symbols('x_0:2')
poly = sy.polys.Poly.from_dict({(2,1):3, (3,0):-1, (0,3):5, (1,2):-2}, x)
W1_ht, gens = leading_term(W1, poly, x)
W2_ht, gens = leading_term(W2, poly, x)


print(poly)
print(sy.polys.monomials.Monomial(tuple(W1_ht), gens))
print(sy.polys.monomials.Monomial(tuple(W2_ht), gens))


