""" Helper functions related to maths

This is not yet a proper library, so a proper description will be created later.

Author: Lorenzo Rota
"""

from sage.all import *
import itertools


def constrained_partitions(n, k, min_elem, max_elem):
    """Iterable object over all non-negative integer decomposition k-tuples
    whose elements sum up to n.
    
    Args:
        n: the integer to be decomposed
        k: the number of integers to decompose it into
        min_elem: the smallest integer allowed in the decomposition
        max_elem: the largest integer allowed in the decomposition
    """
    
    allowed = range(max_elem, min_elem-1, -1)

    def helper(n, k, t):
        if k == 0:
            if n == 0:
                yield t
        elif k == 1:
            if n in allowed:
                yield t + (n,)
        elif min_elem * k <= n <= max_elem * k:
            for v in allowed:
                yield from helper(n - v, k - 1, t + (v,))

    return helper(n, k, ())

def constrained_bounded_partitions(n, k):
    """Iterable object over all non-negative integer decomposition k-tuples
    whose elements sum up to be less than n.
    
    Args:
        n: the integer to be decomposed
        k: the number of integers to decompose it into

    """    

    yield from itertools.chain.from_iterable(constrained_partitions(i, k, 0, i) for i in range(0, n))

def enum_tuple_ring(m, R):
    """Iterable object over all m-tuples of elements in the ring R.
    
    Arg:
    	m: The Number of indeterminates.
    	R: A ring.
    """
    
    if m==0:
        yield ()
    else:
        for i in R:
            for L in enum_tuple_ring(m-1, R):
                yield (i,) + L

def all_multi_index_geq(indices, tup, m):
    """Returns all m-tuple multi-indices that are pointwise greater than or equal to a particular m-tuple
    
    Args:
        indices: The vector of multi-indices.
        tup: The tuple that each multi-index is compared to.
        
    Returns:
        out (tuple): The tuple of multi-indices that are greater than or equal to tup

    """
    
    out = []
    for multi_index in indices:
        if all(multi_index[i] >= tup[i] for i in range(m)):
            out.append(multi_index)
    return tuple(out)

def all_multi_index_leq(indices, tup, m):
    """Returns all m-tuple multi-indices that are pointwise less than or equal to a particular m-tuple
    
    Args:
        indices: The vector of multi-indices.
        tup: The tuple that each multi-index is compared to.
        
    Returns:
        out (tuple): The tuple of multi-indices that are less than or equal to tup

    """
    
    out = []
    for multi_index in indices:
        if all(multi_index[i] <= tup[i] for i in range(m)):
            out.append(multi_index)
    return tuple(out)

def all_multi_index_bdd(indices, low, upp, m):
    """Returns all m-tuple multi-indices that are pointwise less than or equal to a particular m-tuple
    
    Args:
        indices: The vector of multi-indices.
        upp: The upper tuple that each multi-index is compared to.
        low: The lower tuple that each multi-index is compared to.
        
    Returns:
        out (tuple): The tuple of multi-indices that are bounded by upp and low

    """
    
    out = []
    for multi_index in indices:
        if all(low[i] <= multi_index[i] <= upp[i] for i in range(m)):
            out.append(multi_index)
    return tuple(out)

def to_GF(q, elem):
    """Convert integer to element from GF(q)
    
    Args:
        q: Field size.
        elem: Integer.
    
    Returns: Element in GF(q)

    """

    if type(GF(q)).__name__ == 'FiniteField_givaro_with_category':
        return GF(q).fetch_int(elem % q)
    else:
        return elem

def from_GF(q, elem):
    """Convert element from GF(q) to integer modulo q

    Args:
        q: Field size.
        elem: Field element.

    Returns:
        new_num: Integer modulo q

    """

    if type(GF(q)).__name__ == 'FiniteField_givaro_with_category':
        return elem.integer_representation()
    else:
        return elem

def GF_tuple_to_decimal(q, m, tup):
    """Convert tuple of radix GF(q) to decimal

    Args:
        q: Field size.
        m: Size of tuple.
        tup: m-tuple of elements in GF(q).

    Returns:
        new_num: Integer representation of tuple.

    """
    new_num = 0
    for j in range(m):
        pos_val = from_GF(q, tup[j])
        new_num += pos_val * q**(m - (j+1))
    return new_num

def GF_tuple_from_decimal(q, dim, num):
    """Convert decimal to tuple of radix GF(q)

    Args:
        q: Field size.
        dim: Size of tuple.
        num: Integer.

    Returns:
        tuple(base_q): m-tuple representation of num with elements in GF(q)

    """

    pos_val = num
    base_q = []
    base_q.append(to_GF(q, pos_val % q))
    for _ in range(dim - 1):
        pos_val = pos_val // q
        base_q.append(to_GF(q, pos_val % q))
    base_q.reverse()
    return tuple(base_q)

def multi_index_binomial(multi_index_a, multi_index_b, m):
    prod = 1
    for i in range(m):
        prod *= binomial(multi_index_a[i], multi_index_b[i])
    return prod

def hasse_derivative(monomial_basis, monomial_indices_dict, coefs, derivative_order, point, m):
    """ Computes the Hasse derivative w.r.t a specified monomial, given by its multi-index
    
    Args:
        monomial_basis (vector): A monomial basis.
        monomial_indices_dict (dict): A dictionary that associates a multi-index (vector) to a position (int) in the monomial_basis.
        coefs (vector): The coefficients in GF(q) of the polynomial.
        derivative_order (vector): A multi-index that indicates the partial derivative.
        point (tuple): The point at which the derivative is evaluated.
        m (int): The number of indeterminates.
    
    Returns:
        eval_deriv: The computed derivative in GF(q)
    
    """
        
    monomial_indices = list(monomial_indices_dict.keys())
    larger_indices = all_multi_index_geq(monomial_indices, derivative_order, m)
    
    # compute the differences between all larger_indices vectors and the derivative_order vector
    diffs = []
    for multi_index in larger_indices:
        diff_vec = multi_index - derivative_order
        diffs.append(diff_vec)
        
    # suppose the coefficients are called f_j in GF(q) and the multi-index binomial is multi_index_binom(j,i) and call it B
    # we now we sum the coefficients B many times (where B is an int) inside GF(q) and express it as f_j * B.
    # moreover, we consider f_j = coefs[index] where index is obtained from the multi-index
    summed_coefs = []
    for multi_index in larger_indices:
        coef = coefs[monomial_indices_dict[multi_index]]
        binomial = multi_index_binomial(multi_index, derivative_order, m)
        summed_coefs.append(coef*binomial)
    summed_coefs = vector(summed_coefs)
    
    select_monomials = []
    for multi_index in diffs:
        monomial = 1
        for i in range(m):
            #monomial *= interpol.variables[i]**multi_index[i]
            monomial *= point[i]**multi_index[i]
        select_monomials.append(monomial)
    select_monomials = vector(select_monomials)
    
    return(summed_coefs * select_monomials)

def comp_indices_dict(poly, m):
    """Construct a lookup table that associates a multi-index to
        the corresponding monomials from a given polynomial.

    Args:
        poly: The polynomial that we want to create this lookup table for.
        m (int): The number of indetermines of the polynomial.

    Returns:
        indices_dict (dict): The lookup table.
    
    """
    if m == 1:
        exponents = [vector([exp]) for exp in poly.exponents()]
    else:
        exponents = [vector(exp) for exp in poly.exponents()]
    for exp in exponents:
        exp.set_immutable()
    return dict(zip(exponents,range(len(exponents))))

def divided_difference(i, k, Xn, Yn, derivs):
    """Computes the k-th divided difference at position i
    
    Examples:
        [y_i y_i+1 ... y_i+k] = ([y_i+1 y_i+2 ... y_i+k] - [y_i y_i+1 ... y_i+k-1]) / x_i+k - x_i
    
    """
    
    # if k+i > |Xn|
    if i + k >= len(Xn):
        return 0
    # 0th divided difference
    elif k == 0:
        return Yn[i]
    elif Xn[i+k] == Xn[i]:
        return (derivs[i//3][k])/factorial(k)
    # recurrence
    else:
        return (divided_difference(i+1, k-1, Xn, Yn, derivs) - divided_difference(i, k-1, Xn, Yn, derivs)) / (Xn[i+k]-Xn[i])

def hermite_interpolation(data, var):
    """Compute Hermite interpolation given a table of data points and their derivatives
    
    Args:
        data (matrix): Columns represent [x f(x) f'(x) ...]
        var (sage symbol): Symbolic variable that is used to construct the polynomial
    """
    
    xn = []
    yn = []
    for i in range(data.nrows()):
        for j in range(data.ncols() - 1):
            xn.append(data[i][0])
            yn.append(data[i][1])
    
    term = 0
    for k in range(3*n):
        diff = divided_difference(0,k,xn,yn,data[:,1:])
        for i in range(k):
            diff *= (var - xn[i])
        term += diff
    return term