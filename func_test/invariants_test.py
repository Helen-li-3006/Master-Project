import sympy as sy
from sympy.polys.monomials import Monomial
import numpy as np
from funcs.func_utils import *
from funcs.hilbert import *
from funcs.groebner2 import *
from funcs.invariants import *
sy.init_printing(use_unicode=False, wrap_line=False)

# Continue by using the same example C4, get the Molien series
G = [sy.Matrix([[1,0],[0,1]]), sy.Matrix([[0,1],[-1,0]]), sy.Matrix([[-1,0],[0,-1]]), sy.Matrix([[0,-1],[1,0]])]
id = sy.Matrix([[1,0],[0,1]])
l = sy.symbols('l_:1')
mol = molien(G, id, l)
print('Molien series for C4 is ', hilb_expand(mol, l, [5]))
x = sy.symbols('x_:2')
W = sy.Matrix([[1,1]]) # Elimination ordering
d = [5]
invs = inv_ring(G, mol, l, x, W, d)
print(invs)