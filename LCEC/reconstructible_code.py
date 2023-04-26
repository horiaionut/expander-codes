import io

class ReconstructibleCode:
    """Abstract class defining the main functionality of a reconstructible code"""

    def __init__(self, query_complexity, length, alphabet):
        self.query_complexity = query_complexity
        self.length           = length
        self.alphabet         = alphabet

    def build_query_set(self, position):
        raise NotImplementedError

    def reconstruct(self, position, symbols):
        raise NotImplementedError


class DummyReconstructibleCode(ReconstructibleCode):
    def __init__(self, query_complexity, length, alphabet):
        super().__init__(query_complexity, length, alphabet)

    def build_query_set(self, position):
        return range(self.length)

    def reconstruct(self, position, symbols):
        return self.alphabet[0]