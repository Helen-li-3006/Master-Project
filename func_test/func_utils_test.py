import sympy as sy
import numpy as np
from funcs.func_utils import *
sy.init_printing(use_unicode=False, wrap_line=False)

# Replicate the examples from fig 1.2 from Gatermann
W1 = sy.Matrix([[1,1],[1,2]])
W2 = sy.Matrix([[3,1],[1,-1]])
x = sy.symbols('x_0:2')
poly = sy.polys.Poly.from_dict({(2,1):3, (3,0):-1, (0,3):5, (1,2):-2}, x)
W1_ht = leading_term(W1, poly, x)
W2_ht = leading_term(W2, poly, x)


print('Original polynomial', poly)
print('Leading term wrt W1', sy.polys.monomials.Monomial(W1_ht, x))
print('Leading term wrt W2', sy.polys.monomials.Monomial(W2_ht, x))

# Test for S-polynomials function
x = sy.symbols('x_0:2')
f1 = sy.polys.Poly.from_dict({(3,1):3., (1,1):2., (0,1):-1.},x)
f2 = sy.polys.Poly.from_dict({(1,2):2., (0,3):-5.},x)
# Graded lex order x_0 > x_1
W = sy.Matrix([[1,1],[1,0]]) 
S = S_poly(f1,f2,W, x)
print('The S-polynomial is', S)

# Test for checking divisibility 
x = sy.symbols('x_0:2')
W = sy.Matrix([[1,1], [1,0]]) # Graded lex order
# Example 1: non-divisible case
# f1 = 3x_0^2*x_1 - x_0*x_1^2
f1 = sy.polys.Poly.from_dict({(2,1):3., (1,2):-1.}, x)
# f2 = x^2*x_1^2 - 2x_0^3*x_1^1
f2 = sy.polys.Poly.from_dict({(2,2):1., (3,1):-3.}, x)
print(f1, 'if divisible by ', f2, ': ', is_divisible(f1, f2, W, x))

# Example 2: divisible case
f2 = sy.polys.Poly.from_dict({(2,0):1., (1,1):-3.}, x)
print(f1, 'if divisible by ', f2, ': ', is_divisible(f1, f2, W, x))

# Test for division algorithm
# The 'Groebner basis'
x = sy.symbols('x_0:3')
F = [sy.polys.Poly.from_dict({(1,0,0):1., (0,0,1):1.},x), sy.polys.Poly.from_dict({(0,1,0):1., (0,0,1):-1.},x)]
# Weight system for lex order x_0>x_1>x_2
W = sy.Matrix([[1,0,0], [0,1,0], [0,0,1]])
# Polynomial to be reduced 
f = sy.polys.Poly.from_dict({(1,1,0):1.}, x)
g = normalf(W, F, f, x)
print('Normal form of ', f, ' and ', F, ' is ', g)

# Test for edge case of division algorithm (g=0)
x = sy.symbols('x')
F = [sy.Poly(x, x)]
W = sy.Matrix([1])
f = sy.Poly('x**2', x)
g = normalf(W, F, f, x)
print('Normal form of ', f, ' and ', F, ' is ', g)
