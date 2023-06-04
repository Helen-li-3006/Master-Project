"""
This file will contain the algorithm implementation
for computing Hilbert series and Molien series
"""

import sympy as sy
from sympy.polys.monomials import Monomial
import numpy as np

from funcs.func_utils import *
sy.init_printing(use_unicode=False, wrap_line=False)

def hilbinumerator(Ims, W, gens, l):
    """
    Implementation of algorithm 2 in the report
    (1.2.16 in Gatermann, subroutine)

    Inputs:
    Ims: Monomial ideal to compute hilbert series (list of power tuples)
    W: Sypmy matrix Weight system of multi-grading
    l: Variable as l_0, ... ,l_r-1 for hilbert series

    Output:
    numi: Numerator of the hilbert series in terms of l_1, ... ,l_r
    """
    r = W.shape[0]
    ind = 0
    if len(Ims) == 0:
        numi = 1.
        return numi
    else:
        xalp = Ims[ind]
        ind += 1
        alp_deg = weighted_deg(W, xalp)
        numi = 1. - sy.polys.Poly.from_dict({tuple(alp_deg):1.}, l)
        for j in range(1,len(Ims)):
            Jms = [tuple(np.subtract(tup_lcm(Ims[i], Ims[j]), Ims[j])) for i in range(j)]
            # Find the weighted degree for all tuples in Jms
            J_deg = [tuple(weighted_deg(W, mon)) for mon in Jms]
            # List all non-linear monomials in Jms
            J1ms = [Jms[ind] for ind in range(len(J_deg)) if np.mean(J_deg[ind]) != 1]
            # List all linear monomials in Jms
            J_lin = [Jms[ind]for ind in range(len(J_deg)) if np.mean(J_deg[ind]) == 1]
            numJ1 = hilbinumerator(J1ms, W, gens, l)
            numJ = numJ1
            for lin_term in J_lin:
                # Using result 3, calculate numerator of ideal J
                pow = weighted_deg(W, lin_term)
                numJ = (1. - sy.polys.Poly.from_dict({tuple(pow):1.}, l)) * numJ
            alp_j = weighted_deg(W, Ims[j])
            numi = numi - sy.polys.Poly.from_dict({tuple(alp_j):1.}, l) * numJ
        return numi


def hilbert(Ims, W, gens):
    """
    Implementation of algorithm 3 in the report 
    (1.2.16 in Gatermann)

    Inputs:
    Ims: Monomial ideal to compute hilbert series 
    W: Weight system of multi-grading

    Output:
    hilb: hilbert series of the quotient ring K[x]/Ims as a symp expression
    """
    # Define the symbols for lambdas as l_0, ...,l_r-1
    r = W.shape[0]
    l = sy.symbols('l_0:%d'%r)
    hilb = hilbinumerator(Ims, W, gens, l)
    n = len(gens)
    for i in range(n):
        xi= tuple(np.eye(1,n,i).astype(int).reshape(-1))
        Wij = weighted_deg(W, xi)
        hilb = sy.factor(hilb/(1.- sy.polys.Poly.from_dict({Wij:1.},l)))
    return hilb.as_expr(l), l


# Change this function to be compatible for both multivariate and univariate case
def hilb_expand(hilb, ls, orders):
    """
    Expand the hilbert series up to an order. This can also be used for expanding molien series. 

    Inputs:
    hilb: Hilbert series
    ls: generators of the hilbert series (l_1, ..., l_r)
    orders: list of orders for each variables to expand up to (d_1, ..., d_r)
    
    Output: 
    exp_hilb: Expanded hilbert series with each variable up to d_i
    """
    exp_hilb = hilb
    for ind in range(len(orders)):
        exp_hilb = sy.series(exp_hilb, ls[ind], 0, orders[ind]).removeO()
    return exp_hilb
        


def molien(G, id, l):
    """
    Compute the Molien series for finite group

    Inputs:
    G: List of group elements (finite group) under some respresentation (matrix)
    id: identity of G under the same representation (matrix)
    l: variable for molien series (lambda)
    """
    mol = 0.0
    for rep in G:
        mol += 1.0/((id - l*rep).det())
    return (1/len(G))*mol
