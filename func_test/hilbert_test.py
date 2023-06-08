import sympy as sy
from sympy.polys.monomials import Monomial
import numpy as np
from funcs.func_utils import *
from funcs.hilbert import *
sy.init_printing(use_unicode=False, wrap_line=False)

# Replicate of example with single weight from 'Hilbert Functions and the Buchberger Algorithm' (CARLO TRAVERSO)
x = sy.symbols('x_0:4')
Ims = [(0,0,3,0), (0,0,0,4)]
# Example uses Lex ordering
W = sy.Matrix([[1,1,1,1]])
hilb, l = hilbert(Ims, W, x)
hilbnum = hilbinumerator(Ims, W, x, l)
print('Numerator of the Hilbert series is:', hilbnum)
print('Full Hilbert series is ', hilb)

# Example of multi-degree hilbert function from Multigraded Hilbert Functions and Buchberger Algorithm
x = sy.symbols('x_0:4')
W = sy.Matrix([[1,1,2,1], [3,1,1,1]])
# leading terms of the ideal 
Ims = [(1,0,0,1), (0,1,0,1)]
hilb, l = hilbert(Ims, W, x)
hilbnum = hilbinumerator(Ims, W, x, l)
print('Numerator of the Hilbert series is:', hilbnum)
print('Full Hilbert series is ', hilb)

# Test for Molien series using C4 under 2x2 matrix representation
G = [sy.Matrix([[1,0],[0,1]]), sy.Matrix([[0,1],[-1,0]]), sy.Matrix([[-1,0],[0,-1]]), sy.Matrix([[0,-1],[1,0]])]
id = sy.Matrix([[1,0],[0,1]])
l = sy.symbols('l_:1')
mol = molien(G, id, l)
print('Molien series for C4 is ', mol)
print('Molien series up to degree 4 is', hilb_expand(mol, l, [5]))
