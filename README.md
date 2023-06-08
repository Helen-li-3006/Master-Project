# Master-Project
Repository that contains a package built with Sympy that implements algorithms to find ring of invariant polynomials for a polynomial ring for a group representation. The modules are separated by functionality (in the order of implementation):
func_utils ---------- contains all utility functions that is used in matrix form of graded algebra
hilbert ------------- contains the routine (hilbert series) and subruotine (hilbert series numerator) 
                      for computing hilbert series for a given weight system
groebner2 ----------- contains the function and helper functions used to implement algorithm 1.2.23 for 
                      computing a multi-truncated groebner basis up to multi-degree d given any weight 
                      system
invariants ---------- contains the function and helper fnctions used to implement algorithm 2.1.10 for 
                      computing the set of invariants (for a finite group) up to degree d using Molien 
                      series 