import sympy as sy
import numpy as np
from funcs.hilbert import *
from funcs.func_utils import *
sy.init_printing(use_unicode=False, wrap_line=False)

x = sy.symbols('x_1:4')
Ims = [(2,0,2), (0,2,2), (1,0,3), (2,1,0), (1,2,0), (1,1,1)]
W = sy.Matrix([[1,1,1],[2,2,2]])
r = W.shape[0]
l = sy.symbols('l_0:%d'%r)
print(np.eye(1,3,1).reshape(-1))
for i in range(3):
    print(sy.polys.Poly.from_dict({tuple(np.eye(1,3,i).astype(int).reshape(-1)):-1.},x))