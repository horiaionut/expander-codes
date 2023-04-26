""" Custom multivariate polynomial interpolation library

This is not yet a proper library, so a proper description will be created later

Author: Lorenzo Rota.
"""

from sage.all import *
from reconstructible_code import ReconstructibleCode
import random

class ReedMullerCode(ReconstructibleCode):
    """A container that handles local decoding for Reed-Muller codes
    through polynomial interpolation.
    
    The multivariate polynomial that encodes a message is determined by the 
    parameters m, d, q, and allows for k = (m+d) choose d elements in q to be 
    encoded. Typically a value of q = 256 would be chosen such that each
    ascii character can be mapped to a field element in the polynomial ring.
    
    Example: Pick m=3,d=17,q=256 to encode a message of 1140 bytes.

    A local decoder

    Attributes:
        m (int): The number of indeterminates.
        d (int): The maximum degree of each monomial.
        q (int): The finite field order.
        k (int): The total number of monomials.
        query_complexity (int): The number of queries used in the local reconstruction.
        length (int): Length of the codeword.
        alphabet (list): The list of symbols from GF(q).
        interpolating_set (list): The collection of sets of k points (m-tuples)
            that map to the message.
        reconstruction_set (list): The collection of reconstructed symbols from the interpolating set.
        msg (str): The message to be encoded or which has been retrieved
        polynomial (sage polynomial): The polynomial that stores the message.

    """

    def __init__(self, m, d, q):
        self.m = m
        self.d = d
        self.q = q
        self.k = binomial(m+d,d)
        self.query_complexity = q-1
        self.length = q**m
        self.alphabet = list(GF(q))
        _, self.variables = GF(q)[','.join('x%s'%i for i in range(m))].objgens()


    def build_query_set(self, pos):
        if pos < 0 or pos >= self.q**self.m:
            raise ValueError("position needs to be in between 0 and " + str(self.q**self.m) + " but it is " + str(pos))

        point        = self.pos_to_point(pos)
        rand_dir     = vector(GF(self.q), [random.choice(self.alphabet[1:]) for _ in range(self.m)])
        query_points = [point + (t * rand_dir) for t in GF(self.q)][1:]
        query_set    = [self.point_to_pos(point) for point in query_points]

        return query_set


    def reconstruct(self, pos, symbols):
        """Locally reconstructs a symbol given the evaluations from a symbol set
        """
        m, d, q = self.m, self.d, self.q

        if pos < 0 or pos >= q**m:
            raise ValueError("position needs to be in between 0 and " + str(q**m) + " but it is " + str(pos))

        point = self.pos_to_point(pos)
        evals = vector(GF(q), symbols)
        C = codes.ReedSolomonCode(GF(q), Integer(q-1), Integer(d+1))
        decoder = C.decoder("BerlekampWelch")
        recon_poly = decoder.decode_to_message(evals)
        recon_symbol = recon_poly(0)

        return recon_symbol
    

    def monomial_basis(self):
        """ Generates a monomial basis.
        
        This function will be useful for constructing polynomials with known coefficients
        or to construct generalized Vandermonde matrices to solve for the coefficients.
        
        Returns:
            monomials (vector): A monomial basis.
            indices_dict (dict): A dictionary that associates multi-indices with the corresponding monomial basis position
            
        """
        
        m, d, q = self.m, self.d, self.q
        monomials = []
        counter = 0
        for n in range(d+1):
            for powers in constrained_partitions(n, m, 0, n):
                prod = 1
                for i, j in enumerate(powers):
                    prod = prod * self.variables[i]**j
                monomials.append(prod)

        return vector(monomials)


    def random_codeword(self):
        monomials = self.monomial_basis()
        coefs = vector(GF(self.q), [random.choice(self.alphabet) for _ in range(self.k)])
        polynomial = monomials * coefs

        codeword = []
        for i in range(self.q**self.m):
            codeword.append(polynomial(*self.pos_to_point(i)))

        return codeword


    def pos_to_point(self, pos):
        q, m  = self.q, self.m
        point = vector(GF(q), m)

        for idx in range(m):
            point[idx] = GF(q).fetch_int(pos % q)
            pos = pos // q

        return point


    def point_to_pos(self, point):
        pos = 0
        q   = 1

        for idx in range(self.m):
            pos += point[idx].integer_representation() * q
            q   *= self.q
        
        return pos


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