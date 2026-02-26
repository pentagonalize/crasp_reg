import dfa
from dfa import generate_dfa_diagram
from automata.fa.dfa import DFA
from automata.fa.nfa import NFA
import networkx as nx
from sympy import Matrix
import matplotlib.pyplot as plt

def automata_to_graph(my_dfa):
    # converts a dfa into a networkx graph
    # note that we use a MultiDiGraph
    states = my_dfa.states
    transitions = my_dfa.transitions
    G = nx.MultiDiGraph()
    state_list = list(states)
    state_mapping = {}
    for i in range(len(state_list)):
        state_mapping[state_list[i]] = i
    for q in states:
        for sym in transitions[q].keys():
            G.add_edge(state_mapping[q], state_mapping[transitions[q][sym]], key=sym, label=sym)
    return G

def show_graph(G):
    # show the graph
    # does not disambiguate multi-edges
    pos = nx.spring_layout(G)
    nx.draw_networkx(G, pos, with_labels=True, node_color='lightblue', arrows=True)
    edge_labels = {}
    for u, v, data in G.edges(data=True):
        edge_labels[(u, v)] = data.get('label', '')
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels)
    plt.show()


def get_path_edges(G, path):
    # path is a list of nodes in a multigraph
    # recursively computes all paths through the multigraph through the nodes in order
    # returns a list of lists of edges (n1, n2, symbol)
    if len(path) < 2:
        return []
    edges = []
    first = path[0]
    rest = path[1:]
    successors = list(G.successors(first))
    for succ in successors:
        if succ == rest[0]:
            for sym in G[first][succ]:
                label = G[first][succ][sym]['label']
                if isinstance(label, Matrix):
                    label = tuple(label)
                if len(rest) == 1:
                    edges.append([(first, succ, label)])
                else:
                    sub_paths = get_path_edges(G, rest)
                    for sp in sub_paths:
                        edges.append([(first, succ, label)] + sp)
    return edges

def get_simple_cycles_edges(G):
    # computes all simple cycles in G
    # returns a list of lists of edges (n1, n2, symbol)
    cycles = []
    for cycle in nx.simple_cycles(G):
        edges = get_path_edges(G, cycle + [cycle[0]])
        cycles.extend(edges)
    return cycles

def loop_equations(ordered_alphabet, cycles):
    # to maintain a balanced morphism, all loops must sum to 0
    # each loop defines a linear equation based on its Parikh image
    # we will then find the nullspace of this set of linear equations
    # this gives the basis of all balanced morphisms
    basis = []
    for cycle in cycles:
        vec = [0]*len(ordered_alphabet)
        for _,_,sym in cycle:
            index = ordered_alphabet.index(sym)
            vec[index] += 1
        basis.append(vec)
    return basis

def relabel(G, cycles, morphism):
    # morphism is a dict for relabling symbols
    # the labels will be vectors
    # choose an arbitrary start node
    # for each edge n1 -> n2 with label sym
    # pick a path from start to n1
    # the new label is the sum of morphism[sym]
    # and the sum of morphism[label] for edges in the path


    # copy the graph to avoid modifying while iterating
    start_node = list(G.nodes())[0]
    graph_copy = G.copy()
    for u, v, data in graph_copy.edges(data=True):
        label = data['label']
        label = tuple(label) if isinstance(label, Matrix) else label
        # find a path from start_node to u (they all must end up the same)
        path = nx.shortest_path(G, source=start_node, target=u)
        path_edges = get_path_edges(G, path)
        vec = morphism[label]
        if path_edges:
            # if there are paths, sum up the labels under the morphism
            # take the first (they all must end up the same)
            path_edges = path_edges[0]
            for e in path_edges:
                vec += morphism[e[2]]

        # update the label of this edge in graph_copy by column-wise concatenating vec
        for key in graph_copy[u][v]:
            edge_label = graph_copy[u][v][key]['label']
            edge_label = tuple(edge_label) if isinstance(edge_label, Matrix) else edge_label
            if edge_label == label:
                    if not isinstance(edge_label, tuple):
                        # first relabling of just symbols
                        graph_copy[u][v][key]['label'] = vec
                    else:
                        # need to keep history of labels, so we concat label vectors
                        temp_g_label = G[u][v][key]['label'].copy()
                        temp_g_label.col_join(vec)
                        graph_copy[u][v][key]['label'] = temp_g_label
    # copy back into G
    for u, v, key, data in graph_copy.edges(keys=True, data=True):
        G[u][v][key]['label'] = data['label']

def separated(G):
    # unseparated nodes are those with incoming edges having the same label
    # the graph is separated by the labels if there are no unseparated nodes
    labels = {}
    for _, v, data in G.edges(data=True):
        label = data['label']
        label = tuple(label) if isinstance(label, Matrix) else label
        if label in labels:
            if not labels[label] == v:
                return False
        else:
            labels[label] = v
    return True

def get_ordered_alphabet(G):
    # obtain an ordered list of the labels of the edges
    # the ordering is used to order the variables in linear equations
    edge_list = list(G.edges(data=True))
    hashable_labels = []
    for u, v, data in edge_list:
        label = data['label']
        if isinstance(label, Matrix):
            label = tuple(label)
        hashable_labels.append(label)
    ordered_alphabet = list(set(hashable_labels))
    return ordered_alphabet

def attack_scc(G):
    # attack an SCC (to decide membership in C-RASP)
    prev_nullspace = None
    # the algorithm converges once the null spaces don't change
    while True:
        # find all balanced morphisms using the nullspace of the loop equations
        ordered_alphabet = get_ordered_alphabet(G)
        cycles = get_simple_cycles_edges(G)
        basis = loop_equations(ordered_alphabet, cycles)
        A = Matrix(basis)
        current_nullspace = A.nullspace()
        # print("Current nullspace basis:", current_nullspace)
        if current_nullspace == prev_nullspace:
            # converged
            break
        if current_nullspace == []:
            # no more nontrivial balanced morphisms
            break
        else:
            prev_nullspace = current_nullspace
        # relabel the graph using the morphism defined by the nullspace basis
        morphism = {}
        for i in range(len(ordered_alphabet)):
            vec = []
            for nb in current_nullspace:
                vec.append(nb[i])
            morphism[ordered_alphabet[i]] = Matrix(vec)
        relabel(G, cycles, morphism)
    # after convergence, return whether the morphism separates the nodes
    return separated(G)

def decide_CRASP_membership(my_dfa):
    # decide membership in C-RASP by iterating over the SCCs in the dfa
    # the entire dfa is in C-RASP iff every SCC is in C-RASP
    G = automata_to_graph(my_dfa)
    sccs = nx.strongly_connected_components(G)
    for component in sccs:
        sg = G.subgraph(component)
        result = attack_scc(sg)
        if not result:
            return False
    return True

def decide_CRASP_membership_from_regex(regex_str):
    nfa = NFA.from_regex(regex_str)
    my_dfa = DFA.from_nfa(nfa, minify=True)
    return decide_CRASP_membership(my_dfa)

# if __name__ == "__main__":
#     # my_dfa = DFA(
#     #     states={'q0', 'q1', 'q2', 'q3', 'q4', 'q5', 'q6'},
#     #     input_symbols={'a', 'b'},
#     #     transitions={
#     #         'q0': {'a': 'q1', 'b': 'q2'},
#     #         'q1': {'a': 'q6', 'b': 'q0'},
#     #         'q2': {'a': 'q6', 'b': 'q3'},
#     #         'q3': {'a': 'q4', 'b': 'q6'},
#     #         'q4': {'a': 'q5', 'b': 'q6'},
#     #         'q5': {'a': 'q6', 'b': 'q2'},
#     #         'q6': {'a': 'q6', 'b': 'q6'}
#     #     },
#     #     initial_state='q0',
#     #     final_states={'q3'}
#     # )
#     regex_str = '(ab|aabb)*'
#     print("regex:", regex_str)
#     nfa = NFA.from_regex(regex_str)
#     my_dfa = DFA.from_nfa(nfa, minify=True)

#     generate_dfa_diagram(my_dfa, filename='drawings/my_dfa', output_format='svg', auto_open=False)
#     G = automata_to_graph(my_dfa)
#     plt.show()

#     membership = decide_CRASP_membership(my_dfa)
#     print(f"CRASP membership: {membership}")

