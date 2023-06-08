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
    return tuple(W * sy.Matrix(pow))

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
        if max(np.subtract(left[index_left],right[index_right])) <= 0:
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
    order induced by the weight system via merge sort (efficient sorting 
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

def leading_term(W, poly, gens, ret_all = False, inc_coef = False):
    """
    For a polynomial, compute the leading term based on a weight system

    Input: 
    W : mxn matrix of valid weight system
    poly: Sympy polynomial expression 
    gens: n-len list of variables used by the poly

    Output: 
    ht: tuple of powers of the leading term monomial
    """
    # Decompose the polynomial into a list of monomials
    monos = [mon for mon in poly.monoms()]
    degs = [weighted_deg(W, pow) for pow in monos] # Weighted degrees

    # Select and return monomial of the largest weighted deg
    if sy.shape(W)[0] > 1:
        sorted_deg = mod_merge_sort(degs)[::-1]
        sorted_monom = [monos[degs.index(i)] for i in sorted_deg]
    else:
        degs = [deg[0] for deg in degs]
        sorted_deg = sorted(degs)[::-1]
        sorted_monom = [monos[degs.index(i)] for i in sorted_deg]
    if ret_all == True:
        return sorted_monom
    if inc_coef:
        # Return a dictionary of {ht(f):hc(f)} (tuple:coef)
        return {sorted_monom[0]:poly.as_dict()[sorted_monom[0]]}
    else:
        return sorted_monom[0]

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

def S_poly(f, g, W, gens, dom = None):
    """
    Compute the S-polynomial of polynomials f and g wrt weight system 
    W (which represents a term ordering)

    Inputs:
    f,g: Sympy polynomials with n variables
    W: nxr Sympy matrix that is a weight system
    gens: list of n variables defined in sympy class

    Output:
    S: S-polynomial S(f,g), sympy polynomial class
    """
    # Change into dictionary class for ease of computation
    f_dict = f.as_dict()
    g_dict = g.as_dict()
    # Leading terms as power tuple
    f_lt = leading_term(W, f, gens) # power tuple
    g_lt = leading_term(W, g, gens) # power tuple
    lt_lcm = tup_lcm(f_lt, g_lt) #power tuple
    num = sy.polys.Poly.from_dict({lt_lcm: 1.}, gens=gens, domain=dom)
    S = (num/sy.polys.Poly.from_dict({f_lt:f_dict.get(f_lt)},gens=gens, domain=dom))*f - (num/sy.polys.Poly.from_dict({g_lt:g_dict.get(g_lt)},gens=gens, domain=dom))*g
    if dom:
        return S.set_domain(dom)
    else:
        return S

def is_divisible(f1, f2, W, gens):
    """
    Return boolean variable for if ht(f2)|ht(f1)
    Inputs:
    f1, f2: Sympy polynomials with respect to the variables gens
    W: Sympy matrix that is a weight system 
    gens: variables x_i

    Output:
    is_div: Boolean if statement is satisfied
    """
    ht1 = leading_term(W, f1, gens)
    ht2 = leading_term(W, f2, gens)
    diff_bool = [0 if ht1[i]>=ht2[i] else 1 for i in range(len(ht1))]
    if sum(diff_bool) > 0:
        return False
    else:
        return True

def normalf(W, F, f, gens, dom=None):
    """
    Impelemnt division algorithm to check if polynomial f is in the span of F

    Input:
    f: The polynomial in sympy polynomial form - polynomial to apply top reduction
    F: List of sympy polynomials - 'Groebner basis'
    gens: List of symbols 
    dom: Domain of polynomial reduction

    Output:
    g: normal form of a polynomial w.r.t F
    """
    # Initialise variables
    g = f
    ind = 0
    while g != 0 and ind < len(F):
        if is_divisible(f,F[ind], W, gens):
            f1_dict = g.as_dict()
            f2_dict = F[ind].as_dict()
            f1_ht = leading_term(W, g, gens)
            f2_ht = leading_term(W,F[ind], gens)
            coef = f1_dict.get(f1_ht)/f2_dict.get(f2_ht)
            lcm_tup = tup_lcm(f1_ht, f2_ht)
            mult = tuple(map(lambda i, j: i - j, lcm_tup, f2_ht))
            g -= sy.polys.Poly.from_dict({mult:coef}, gens=gens, domain=dom) * F[ind]
            ind += 1
        else:
            ind += 1
            pass
    return g

def min_deg(poly, gens):
    """
    Return the minimum degree (vector for multigrading) for a polynomial (hilbert function)
    with respect to grading
    
    Inputs:
    poly: Hilbert function as a polynomial
    Outputs: 
    Power tuple of the minimal degree term

    """
    pows = [np.array(key) for key in poly.as_dict().keys()]
    pows_inc = mod_merge_sort(pows)
    return tuple(pows_inc[0]) #return minimal power tuple

def compare_degs(deg1, deg2):
    """
    Compare 2 weighted degrees with respect to the same weight system
    Inputs:
    deg1, deg2: tuples of the same length to compare degrees 
    
    Outputs:
    return True is deg1 <= deg2 and False o/w
    """
    diff = np.subtract(deg2, deg1)
    return max(diff) >= 0


