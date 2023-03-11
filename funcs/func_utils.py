import sympy as sy
import numpy as np
sy.init_printing(use_unicode=False, wrap_line=False)

def weighted_deg(W, pow):
    """
    Given a weight system (multigrading when W has multiple rows)
    and a monomials, calculate the weighted degree
    under the weight system. 
    
    Input:
    W: mxn matrix, a valid weight system
    pow: n-length tuple, powers of a monomial

    Output:
    m-length array, weighted degree of monomial under each sub-grading
    """
    return W * sy.Matrix(pow)

def mod_merge(left, right):
    """
    Helper function for merge sort to be performed iteratively
    """
    if len(left) == 0 :
        return right
    
    if len(right) == 0:
        return left
    
    result = []
    index_left = index_right = 0

    while len(result) < len(left) + len(right):
        if max(left[index_left] - right[index_right]) <= 0:
            result.append(left[index_left])
            index_left += 1
        else:
            result.append(right[index_right])
            index_right += 1
        if index_left == len(left):
            result += right[index_right:]
            break
        if index_right == len(right):
            result += left[index_left:]
            break
    return result

def mod_merge_sort(deg):
    """
    Take a list of sympy arrays (weighted degree) and sort by the term 
    order induced by the weigh system via merge sort (efficient sorting 
    algorithm on python)

    Input:
    deg: list of weighted degrees

    Output:
    Sorted deg into ascending order, ending with the monomial of largest 
    weighted degree
    """
    if len(deg) < 2:
        return deg
    
    midpoint = len(deg)//2
    return mod_merge(mod_merge_sort(deg[:midpoint]), mod_merge_sort(deg[midpoint:]))

def leading_term(W, poly, gens):
    """
    For a polynomial, compute the leading term based on a weight system

    Input: 
    W : mxn matrix of valid weight system
    poly: Sympy polynomial expression 
    gens: n-len list of variables used by the poly

    Output: 
    ht: tuple of powers of the leading term monomial
    gens: n-len list of variables 
    """
    # Decompose the polynomial into a list of monomials
    monos = [mon for mon in poly.monoms()]
    degs = [weighted_deg(W, pow) for pow in monos]

    # Select and return monomial of the largest weighted deg
    if sy.shape(W)[0] > 1:
        sorted_deg = mod_merge_sort(degs)
        return monos[degs.index(sorted_deg[-1])], gens
    else:
        max_id = np.argmax(degs)
        ht = monos[max_id]
    return ht, gens

def tup_lcm(mon1, mon2):
    """
    Function that returns LCM of 2 monomials as a tple of powers

    Inputs:
    mon1, mon2: monomials powers in form of tuples

    Output:
    lcm: the LCM monomial in form of tuple of powers
    """
    lcm = []
    for i in range(len(mon1)):
        lcm.append(np.max([mon1[i], mon2[i]]))
    return tuple(lcm)
