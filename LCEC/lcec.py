import math 

from sage.graphs.bipartite_graph import BipartiteGraph

from reconstructible_code import ReconstructibleCode

class LCECode():
    """Class implementing the localy correctable expander code."""

    class BGraph:
        class Edge:
            def __init__(self, pos):
                self.global_pos = pos[2]
                self.local_pos  = (pos[0], pos[1])
                self.vertices   = [None, None]

        class Vertex:
            def __init__(self, degree):
                self.edges = [None] * degree

        def __init__(self, double_cover, lcec):
            self.edges    = [None] * double_cover.size()
            self.vertices = [None] * (2  * len(self.edges) // lcec.dc_degree)

            idx = 0
            for vertex in double_cover.left:
                self.vertices[idx] = LCECode.BGraph.Vertex(lcec.dc_degree)

                for e in double_cover.edges_incident(vertex):
                    if self.edges[e[2][2]] != None:
                        raise Exception('Two different edges can not have the same global position.')
                    
                    if self.vertices[idx].edges[e[2][0]] != None:
                        raise Exception('Two different edges can not have the same local position.')

                    self.edges[e[2][2]] = LCECode.BGraph.Edge(e[2])
                    self.vertices[idx].edges[e[2][0]] = self.edges[e[2][2]]
                    self.edges[e[2][2]].vertices[0] = self.vertices[idx]

                idx += 1

            for vertex in double_cover.right:
                self.vertices[idx] = LCECode.BGraph.Vertex(lcec.dc_degree)

                for e in double_cover.edges_incident(vertex):
                    if self.vertices[idx].edges[e[2][1]] != None:
                        raise Exception('Two different edges can not have the same local position.')

                    self.vertices[idx].edges[e[2][1]] = self.edges[e[2][2]]
                    self.edges[e[2][2]].vertices[1] = self.vertices[idx]

                idx += 1


    class Tree:
        def __init__(self, edge, lcec):
            self.children   = []
            self.gedge      = edge
            self.lcec       = lcec
            self.symbol     = None

            if edge == None:
                raise ValueError('None edge detected.')
        

        def make_tree(self, side, depth):
            assert(depth >= 1)

            if depth == 1:
                return [self]

            leaves = []

            vertex = self.gedge.vertices[side]

            pos_to_be_reconstructed = self.gedge.local_pos[side]

            query_set = self.lcec.inner_code.build_query_set(pos_to_be_reconstructed)

            for query_pos in query_set:
                
                child = LCECode.Tree(vertex.edges[query_pos], self.lcec)
                self.add_child(child)

                leaves += child.make_tree(1 - side, depth - 1)
            
            return leaves


        def correct_lower_tree(self, side):
            for child in self.children:
                child.correct_lower_tree(1 - side)

            lcec = self.lcec
            self.path_dist  = dict()

            if not self.children:
                for symbol in lcec.alphabet:
                    self.path_dist[symbol] = int(symbol != self.symbol)
            else:
                if self.gedge.local_pos[side] not in lcec.inner_code_rec_results:
                    lcec.inner_code_rec_results[self.gedge.local_pos[side]] = dict()

                    build_symbol_to_sets(lcec.inner_code,
                                self.gedge.local_pos[side], 
                                lcec.inner_code_rec_results[self.gedge.local_pos[side]])

                symbol_to_sets = lcec.inner_code_rec_results[self.gedge.local_pos[side]]

                for reconstructed_symbol in lcec.alphabet:
                    min_path_dist = float('inf')

                    if reconstructed_symbol in symbol_to_sets:
                        for symbols in symbol_to_sets[reconstructed_symbol]:
                            max_path_dist = 0

                            for idx, symbol in enumerate(symbols):
                                max_path_dist = max(max_path_dist, 
                                                    self.children[idx].path_dist[symbol])
                            
                            min_path_dist = min(max_path_dist, min_path_dist)
                    
                    self.path_dist[reconstructed_symbol] = min_path_dist + 1

                min_path_dist = self.path_dist[lcec.alphabet[0]] + 1

                for reconstructed_symbol in lcec.alphabet:
                    if self.path_dist[reconstructed_symbol] < min_path_dist:
                        min_path_dist = self.path_dist[reconstructed_symbol]
                        self.symbol   = reconstructed_symbol


        def reconstruct_top_tree(self, side):
            if self.symbol != None:
                return
        
            symbols = []    
            for child in self.children:
                child.reconstruct_top_tree(1 - side)
                symbols.append(child.symbol)

            self.symbol = self.lcec.inner_code.reconstruct(self.gedge.local_pos[side], symbols)


        def add_child(self, node):
            assert isinstance(node, LCECode.Tree)
            self.children.append(node)


        def add_positions_to_dict(self, pos_to_symbol):
            pos_to_symbol[self.gedge.global_pos] = None
            
            for child in self.children:
                child.add_positions_to_dict(pos_to_symbol)


        def fill_symbols(self, pos_to_symbol):
            self.symbol = pos_to_symbol[self.gedge.global_pos]

            for child in self.children:
                child.fill_symbols(pos_to_symbol)


    def __init__(self, inner_code : ReconstructibleCode, double_cover : BipartiteGraph):
        if not double_cover.is_regular():
            raise ValueError("the double cover needs to be regular.")
        
        self.dc_degree = double_cover.degree()[0]

        if self.dc_degree != inner_code.length:
            raise ValueError("the double cover's regularity degree is not equal to the length of the inner code (" + 
                                str(self.cd_degree) + " != " + str(inner_code.length) + ").")

        self.inner_code = inner_code
        self.alphabet   = inner_code.alphabet
        self.length     = 2 * double_cover.size()
        self.dc         = LCECode.BGraph(double_cover, self)

        self.inner_code_rec_results = dict()


    def correct(self, pos, codeword_reader):
        root, top_leaves, bottom_leaves = self.build_reconstruction_tree(pos)
        
        self.fill_bottom_trees_with_symbols(top_leaves, codeword_reader)

        for vertex in top_leaves:
            vertex.correct_lower_tree(1)

        root.reconstruct_top_tree(1)
        return root.symbol


    def build_reconstruction_tree(self, pos):
        depth1 = max(2, int(math.log(len(self.dc.vertices))) // self.dc_degree)
        #TODO: find the right value
        depth2 = depth1

        edge        = self.dc.edges[pos]
        root        = LCECode.Tree(edge, self)
        top_leaves  = root.make_tree(1, depth1)

        bottom_leaves = []
        for leaf in top_leaves:
            bottom_leaves += leaf.make_tree(1, depth2)

        return root, top_leaves, bottom_leaves


    def fill_bottom_trees_with_symbols(self, roots, codeword_reader):
        pos_to_symbol = dict()
        for vertex in roots:
            vertex.add_positions_to_dict(pos_to_symbol)
        
        codeword_reader.fill(pos_to_symbol)

        for vertex in roots:
            vertex.fill_symbols(pos_to_symbol)


def build_symbol_to_sets(reconstructible_code, pos, symbol_to_sets, symbols=[]):
    if len(symbols) == reconstructible_code.query_complexity:
        symbol = reconstructible_code.reconstruct(pos, symbols)
        
        if symbol not in symbol_to_sets:
            symbol_to_sets[symbol] = []
        
        symbol_to_sets[symbol].append(copy(symbols))
    else:
        for symbol in reconstructible_code.alphabet:
            symbols.append(symbol)
            build_symbol_to_sets(reconstructible_code, pos, symbol_to_sets, symbols)
            symbols.pop()