class CodewordReader:
    def __init__(self, codeword):
        self.codeword = codeword

    def fill(self, pos_to_symbol):
        for pos in pos_to_symbol:
            print('Reading pos: ', pos, 'from codeword. Symbol: ', self.codeword[pos])
            pos_to_symbol[pos] = self.codeword[pos]