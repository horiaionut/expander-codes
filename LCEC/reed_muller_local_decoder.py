""" Custom multivariate polynomial interpolation library

This is not yet a proper library, so a proper description will be created later

Author: Lorenzo Rota.
"""

from sage.all import *
from helpers.math import *
from helpers.misc import *
from multivar_poly.abstract_reconstructible_code import AbstractReconstructibleCode
import random
import concurrent.futures


class ReedMullerLocalDecoder(AbstractReconstructibleCode):
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
        self.alphabet = list()
        self.interpolating_set = None
        self.reconstruction_set = None
        self.msg = ''
        # construct alphabet
        for elem in GF(q):
            self.alphabet.append(elem)


    def GF_to_int_tuple(self, lst, m=0):
        """ Converts the field elements of the interpolating set from GF(q) to an integer modulo q """
        
        q, k = self.q, self.k
        new_lst = []
        
        if m == 0:
            for elem in lst:
                new_lst.append(from_GF(q, elem))
        else:
            for i in range(k):
                tup = []
                for j in range(m):
                    tup.append(from_GF(q, lst[i][j]))
                new_lst.append(tuple(tup))
        
        return tuple(new_lst)


    def int_to_GF_tuple(self, lst, m=0):
        """ Converts the integers of the interpolating set to elements from GF(q) """
            
        q, k = self.q, self.k
        new_lst = []
        
        if m == 0:
            for elem in lst:
                new_lst.append(to_GF(q, elem))
        else:
            for i in range(k):
                tup = []
                for j in range(m):
                    tup.append(to_GF(q, lst[i][j]))
                new_lst.append(tuple(tup))

        return tuple(new_lst)


    def resolve_query(self, point):
        """Resolve the server and index associated with some query point

        Args:
            point (m-tuple): The point of evaluation.
       
        Returns:
            server_no: The server that stores the query evaluation.
            idx: The index at which the evaluation is stored.

        """

        m, d, q = self.m, self.d, self.q

        # determine server and local index
        server_no = from_GF(q, point[-1])
        idx = GF_tuple_to_decimal(q, m, point) // q

        return server_no, idx


    def carry_out_query(self, server_no, idx, block_number):
        """Carries out the query by reading from the respective database and returning the evaluation

        Args:
            server_no (int): The server that stores part of the database.
            idx (int): The index at which the evaluation is stored.
            block_number (int): The location that stores a particular encoding of the db.
        
        Returns:
            eval (GF(q)): The result of the query. 

        """
        # SERVER SHOULD MAYBE DO THIS

        q = self.q
        div_count = int(log(256) / log(q)) # represents the number of segments a byte can be divided into given the size q

        with open('servers/server{}/{}.bin'.format(server_no, block_number), 'rb') as file:
            file.seek(idx // div_count, 0)
            read_byte = int.from_bytes(file.read(1), "big")
            read_byte = GF_tuple_from_decimal(q, div_count, read_byte)
        # if index corresponds to leading symbol
        if idx % 2 == 0:
            symbol = read_byte[0]
        else:
            symbol = read_byte[1]
        
        return symbol


    def test_queries(self):
        """Test all queries and see if they correspond to the polynomial evaluation.
        """

        m, d, q = self.m, self.d, self.q

        error = False
        fast_poly = fast_callable(self.polynomial, vars=self.variables)
        for i in range(q**(m-1)):
            # resolving data
            point = self.interpolating_set[i]
            
            server_no, idx= self.resolve_query(point)
            symbol = self.carry_out_query(server_no, idx, 0)
            
            if fast_poly(*point) != symbol:
                error = True
        if not error:
            print('Success')


    def build_query_set(self, position):
        """Builds a query set, which in the case of Reed-Muller is a line to which
        the polynomial is restricted.

        """
        m, d, q = self.m, self.d, self.q

        position_vec = vector(GF(q), position)

        for attempt in range(100):
            rand_dir = vector(GF(q), GF_tuple_from_decimal(q, m, random.randint(1, self.length)))
            # check that rand_dir is not embedded in any hyperplane
            if rand_dir[-1] != to_GF(q, 0):
                break
        if attempt == 99:
            raise Exception('Unable to find appropriate query set')

        query_set = [position_vec + (t * rand_dir) for t in GF(q)][1:]
        return query_set


    def eval_query_set(self, query_set):
        """Evaluates the query sets
        """

        m, d, q = self.m, self.d, self.q

        evaluations = []
        for point in query_set:
            server_no, idx= self.resolve_query(point)
            symbol = self.carry_out_query(server_no, idx, 0)
            evaluations.append(symbol)

        return evaluations


    def reconstruct(self, position, symbols):
        """Locally reconstructs a symbol given the evaluations from a query set
        """

        m, d, q = self.m, self.d, self.q

        evals = vector(GF(q), symbols)
        C = codes.ReedSolomonCode(GF(q), Integer(q-1), Integer(d+1))
        decoder = C.decoder("BerlekampWelch")
        recon_poly = decoder.decode_to_message(evals)
        recon_symbol = recon_poly(0)

        return recon_symbol
    

    def reconstruct_multiple(self, start=None, end=None):
        """Reconstructs a contiguous subsequence of the encoded message

        Args:
            start (int): The starting point of the subsequence.
            end (int): The ending point of the subsequence.

        Returns:
            reconstruction (list): List of reconstructed points of the subsequence.

        """
        q = self.q
        div_count = int(log(256) / log(q)) # represents the number of times 8 bits can be divided given q

        subsequence_points = []
        if start is None or end is None:
            subsequence_points = self.interpolating_set
        else:
            start_is, end_is = start * div_count, end * div_count
            subsequence_points = self.interpolating_set[start_is:end_is]

        reconstruction = []
        for point in subsequence_points:
            query_set = self.build_query_set(point)
            evaluations = self.eval_query_set(query_set)
            symbol = self.reconstruct(point, evaluations)
            reconstruction.append(symbol)

        return reconstruction
    

    def retrieve_message(self, msg):
        """ Retrieves the message by evaluating the polynomial over the interpolating set"""

        q, k = self.q, self.k
        div_count = int(log(256) / log(q)) # represents the number of times 8 bits can be divided given q

        if k % div_count:
            raise ValueError('Message length should be {} characters long'.format(k // div_count))
        elif 8 % div_count:
            raise ValueError('Discrete logarithm of field size must divide 8')

        msg_in_GF = msg
        msg = ''
        for i in range(len(msg_in_GF) // div_count):
            # convert each subsequence of div_count many elements to a byte (in decimal representation)
            # then convert it to a char for the 'message' representation
            sub_seq = msg_in_GF[i*div_count:(i+1)*div_count]
            byte = GF_tuple_to_decimal(q, div_count, sub_seq)
            msg += chr(byte)
        
        self.msg = msg


    def build_key_from_interpolating_set(self):
        """Saves the interpolating set as a key
        
        Returns:
            key (str): A hex representation of the interpolating set

        """

        key = self.GF_to_int_tuple(self.interpolating_set, self.m)
        
        # create hex dump
        payload = ''
        for tup in key:
            for coord in tup:
                payload += '{:0>2X}'.format(coord)
        key = payload

        return key


    def build_interpolating_set_from_key(self, key):
        """Load the key and generate the corresponding interpolating set."""
        
        m, d, q, k = self.m, self.d, self.q, self.k
        self.interpolating_set = None

        # create interpolating set
        lst_of_tuples = []
        for i in range(k):
            tup = []
            for j in range(m*i, m*(i+1)):
                tup.append(int(key[2*j:2*j+2],16))
            lst_of_tuples.append(tuple(tup))
        
        self.interpolating_set = self.int_to_GF_tuple(lst_of_tuples, self.m)