from lcec import *
from codeword_reader import CodewordReader
from reed_muller import ReedMullerCode


def tryit():
    with open("inner_codeword.obj", 'rb') as file:
        loaded_inner_codeword = pickle.load(file)

    inner_code = ReedMullerCode(2, 2, 4)

    query_set = inner_code.build_query_set(3)
    print(query_set)

    symbols = [inner_codeword[pos] for pos in query_set]
    print(loaded_inner_code.reconstruct(3, symbols))