from .automata import Automaton
from .utils import view_graphviz, save_graphviz
import networkx as nx
import os
import tempfile
import sys
import webbrowser
from .ring import *

class monoidElement(tuple):
    def __repr__(self):
        s = ""
        for x in self:
            s += str(x)
        if s == "":
            return "1"
        else:
            return s
    def __add__(self,other):
        r = list()
        for x in self:
            r.append(x)
        for y in other:
            r.append(y)
        return monoidElement(r)
    def __getslice__(self,i,j):
        return monoidElement(tuple(self).__getslice__(i,j))
def power_set(E):
    if len(E) == 0:
        return [frozenset()]
    x = E.pop()
    l = power_set(E)
    l2 = []
    for U in l:
        l2.append(U.union([x]))
    return l+l2
def semigroup_closure(E):
    complete = False
    while not complete:
        complete = True
        F = E.copy()
        for x in F:
            for y in F:
                if x*y not in E:
                    E.add(x*y)
                    complete = False

def draw_box(box):
    s = '\\begin{tabular}{|'
    for i in range(len(box[0])):
        s = s+'c|'
    s = s + '}\n\\hline\n'
    for x in box:
        s =  s + '\n'
        l = list(x)
        s = s + '$'+str(l.pop())
        while len(l) > 0:
            s = s+'$&$'+str(l.pop())
        s = s+'$\\\\\n\\hline\n'
    s = s.replace("['","")
    s = s.replace("']","")
    s = s.replace("'","")

    return s+'\\end{tabular}'
def draw_box_dot(box,idempotents,colors_list=False):
    s = '"{'
    for L in box:
        for H in L:
            sh = ""
            for a in H:
                if a in idempotents:
                    sh= "<"+str(a)+">"+str(a)+"*,"+sh
                else:
                    sh+= str(a)+","
            s = s+sh[0:len(sh)-1]+'|'
        s = s[0:len(s)-1]+'}|{'
    return s[0:len(s)-2]+'"'

def draw_box_dot_old(box,idempotents,colors_list=False,):
    s = '<<table border="1" cellborder="1" CELLSPACING="0" cellpadding="4">'
    for L in box:
        #uns = u.replace("*","")
        s = s + "<tr>"
        for H in L:
            if not colors_list:
                s = s + "<td>"
            else:
                s = s + "<td bgcolor="+colors_list[H[0]]+">"
            for a in H:
                if a in idempotents:
                    s = s + str(a) +"*,"
                else:
                    s = s + str(a) +","

            s = s[0:len(s)-1] + "</td>"
            #uns = u.replace("*","")
        s = s+"</tr>"
    s = s + '</table>>'
    return s

class TransitionSemiGroup():
    r"""
    Return the transition semigroup of a deterministic automaton.

    INPUT:

    - ``automaton`` -- deterministic automaton

    OUTPUT:

        Transition semigroup

    EXAMPLES::

        sage: from pysemigroup import Automaton, TransitionSemiGroup
        sage: d = {(0, 'a'): [1], (1, 'a'): [0]}
        sage: A = Automaton(d,[0],[1])
        sage: S = TransitionSemiGroup(A)
        sage: S
        Transition SemiGroup of Automaton of 2 states

    """
    def __init__(self, automaton, monoid=True,compute=False, max_size= 0):
        r"""
        EXAMPLES::

            sage: from pysemigroup import Automaton, TransitionSemiGroup
            sage: d = {(0, 'a'): [1], (1, 'a'): [0]}
            sage: A = Automaton(d,[0],[1])
            sage: S = TransitionSemiGroup(A)
            sage: S
            Transition SemiGroup of Automaton of 2 states
        """
        if not isinstance(automaton, Automaton):
            raise TypeError('input(%s) must be an automaton' % automaton)
        self._monoid = monoid
        self._compute = compute
        self._idempotents = False
        self._R_class = dict()
        self._stable_comp = False
        self._J_class = dict()
        self._J_ideal = dict()
        self._L_class = dict()
        self._H_class = dict()
        self._cayley_graphs = dict()
        self._generators = [monoidElement(x) for x in automaton._alphabet]
        d = dict(automaton._transitions)
        self._automaton = Automaton(d,automaton._initial_states,automaton._final_states,aut_type=automaton._type)
        if compute:
            l = self.elements(max_size=max_size)
        self._computed_rep = dict()
    def __repr__(self):
        r"""
        EXAMPLES::

            sage: from pysemigroup import Automaton, TransitionSemiGroup
            sage: d = {(0, 'a'): [1], (1, 'a'): [0]}
            sage: A = Automaton(d,[0],[1])
            sage: S = TransitionSemiGroup(A)
            sage: S
            Transition SemiGroup of Automaton of 2 states
        """
        return "Transition SemiGroup of %s" % self._automaton

    def graphviz_string(self, arrow=True,verbose=False,unfold=True,colors_list=False,get_repr=False):
        r"""
        Return graphviz eggbox representation of self. Set arrow to False to delete the arrow.
        INPUT :
        -  ``self`` -  Automaton
        -  ``arrow`` -  boolean
        -  ``verbose`` -  boolean

        OUTPUT:

        string

        EXAMPLES::
            sage: from pysemigroup import Automaton, TransitionSemiGroup

        """

        #from sage.plot.colors import rainbow
        #col = rainbow(len(colors_list))
        #col_dic = {}

        #for x in colors_list:
        #      for y in x:
        if verbose:
            print("computing box diagramm ...")

        box = self.box_representation(verbose=verbose)
        idempotents = self.idempotents()
        if verbose:
            print("done.")
        repre = set(box)
        Gcal =  nx.condensation(self.cayley_graph(orientation="left_right"))
        if verbose:
            edge_nb = str(len(Gcal.edges()))
            print("computing global structure ...")
        count = 0
        if verbose:
            print("done.")
        Edge = []
        graph_viz = 'digraph {node [shape= record] \n'
        for x in repre:
            if x == '' or x == ():
                rx = '"1"'
            else:
                rx = '"'+str(x)+'"'
            if unfold:
               graph_viz = graph_viz + rx + ' [label='+draw_box_dot(box[x],idempotents,colors_list=colors_list)+'];\n'
            else:
               graph_viz = graph_viz + rx + ' [label='+rx+',shape="rectangle"];\n'

        if not arrow:
            graph_viz = graph_viz + 'edge [style="invis"]\n'
        if verbose:
            print("computing successor edges ...")
            count = 0
            loop_ln = len(repre)^2
        for x in repre:
            for y in repre:
                edge = (Gcal.graph["mapping"][x],Gcal.graph["mapping"][y])
                if edge  in Gcal.edges():
                    if x == '' or x == ():
                        rx = '1'
                    else:
                        rx = str(x)
                    if y == '' or y == ():
                        ry = '1'
                    else:
                        ry = str(y)
                    graph_viz = graph_viz+'"'+rx+'"->"'+ry+'";\n'
        if get_repr:
            return (graph_viz + '}',repre)
        else:
            return graph_viz + '}'


    def __call__(self,word):
        r"""
        Return the representent of word in the semigroup

        EXAMPLES::

            sage: from pysemigroup import Automaton, TransitionSemiGroup
            sage: d = {(0, 'a'): [2], (0, 'b'): [1], (1, 'a'): [0],  (1, 'b'): [2], (2, 'a'): [2], (2, 'b'): [2]}
            sage: A = Automaton(d,[1],[1])
            sage: S = TransitionSemiGroup(A)
            sage: S("abaaa")
            aa

        """
        return self.representent(word)

    def __len__(self):
        r"""

        EXAMPLE::

        """
        return len(self.elements())
    def length(self,maxsize=0):
        r"""

        EXAMPLE::

        """
        E = self.elements(maxsize=maxsize)
        if E == False:
            return maxsize+1
        return len(E)

    def __iter__(self):
        r"""
        Return an iterator of the elements of the semigroup.

        EXAMPLES::

            sage: from pysemigroup import Automaton, TransitionSemiGroup
            sage: d = {(0, 'a'): [2], (0, 'b'): [1], (1, 'a'): [0],  (1, 'b'): [2], (2, 'a'): [2], (2, 'b'): [2]}
            sage: A = Automaton(d,[1],[1])
            sage: S = TransitionSemiGroup(A)
            sage: it = iter(S)
            sage: sorted(it)
            [1, a, aa, ab, b, ba]
        """
        for x in self.elements():
            yield x


    def _compute_semigroup(self, verbose=False):
        if not self._compute:
            E = self.elements(verbose=verbose)

    def elements(self, maxsize= 0, verbose=False):
        r"""
        Compute the transition semigroup of the automaton

        INPUT :

        -  ``verbose`` - boolean
        EXAMPLES::

            sage: from pysemigroup import Automaton, TransitionSemiGroup
            sage: d = {(0, 'a'): [1], (1, 'a'): [0]}
            sage: A = Automaton(d,[0],[1])
            sage: S = TransitionSemiGroup(A)
            sage: S.elements()
            frozenset({1, a})

        """
        if self._compute == True:
            return self._elements
        A = self._generators
        Aut = self._automaton
        Sg = []
        Rep = {}
        Rep_rv = {}
        letter_quot = {}
        for x in Aut._alphabet:
            fctx = Aut.letter_to_algebra(x)
            if not(fctx  in Sg):
                Rep[fctx] = monoidElement(x)
                Rep_rv[monoidElement(x)] = fctx
                Sg.append(fctx)
            else:
                letter_quot[x] = fctx
        Gene = list(Sg)
        last = list(Sg)
        if self._monoid:
            ident = Aut.identity_on_automata_ring()
            if ident not in Sg:
                Rep[ident] = monoidElement("")
                Rep_rv[monoidElement("")] = ident
                Sg.append(ident)
        count = 0
        while len(last)>0:
            old = list(last)
            last = []
            count = count + 1
            if (maxsize>0) and (len(Sg)>maxsize):
                return False
            for u in old:
                for v in Gene:
                    v_u = u*v
                    Ru = Rep[u]
                    Rv = Rep[v]
                    if not(v_u  in Sg):
                        if verbose:
                            print("new element "+str(Ru+Rv)+" no: "+str(len(Sg)))
                        Sg.append(v_u)
                        last.append(v_u)
                        Ruv = Ru+Rv
                        Rep[v_u] = Ruv
                        Rep_rv[Ruv] = v_u

        self._Sg = Sg
        self._Representations = Rep
        self._Representations_rev = Rep_rv
        self._compute = True
        self._letter_quot=letter_quot
        self._elements = frozenset(self._Representations_rev)
        return self._elements
    def get_identity(self):
        try:
            return self._identity
        except:
            k = self.elements()
            l = list(self._Representations)
            for x in l:
                if x == self._automaton.identity_on_automata_ring():
                    self._identity = self._Representations[x]
                    return self._identity
            raise ValueError("No indentity in this semigroup")
    def representent(self, v,verbose=False):
        r"""
        Return a representent of the equivalence class of the word u in the transition semigroup
        p
        INPUT :

        - ``u``  - string or monoid element

        EXAMPLES::

            sage: from pysemigroup import Automaton, TransitionSemiGroup
            sage: d = {(0, 'a'): [1], (1, 'a'): [0]}
            sage: A = Automaton(d,[0],[1])
            sage: S = TransitionSemiGroup(A)
            sage: S.representent('aaa')
            a
        """
        if verbose:
            print(v)
        if v in self._computed_rep:
            return self._computed_rep[v]
        l = self.elements()
        u = monoidElement(v)
        if verbose:
            print(len(u))
        if u in self._Representations.values():
            fu = u
        else:
            if len(u) == 0:
                fu = self.get_identity()
            elif len(u) == 1:
                if str(u) in self._letter_quot:
                    fu = self._Representations[self._letter_quot[str(u)]]
                else:
                    raise TypeError("Letter "+str(v)+" is not in the underlying alphabet: "+str(self._generators))
            else:
                fu_raw = self._Representations_rev[self.representent((v[0],),verbose=verbose)]
                for i in v[1:len(v)]:
                    fu_raw = fu_raw*self._Representations_rev[self.representent((i,),verbose=verbose)]
                fu = self._Representations[fu_raw]
        self._computed_rep[v]=fu
        return fu

    def idempotents(self):
        r"""
        Return the idempotents of the semigroup

        EXAMPLES::

            sage: from pysemigroup import Automaton, TransitionSemiGroup
            sage: d = {(0, 'a'): [1], (1, 'a'): [0]}
            sage: A = Automaton(d,[0],[1])
            sage: S = TransitionSemiGroup(A)
            sage: S.idempotents()
            {1}

        """
        self._compute_semigroup()
        if not(self._idempotents):
            I = set()
            for u in self.elements():
                if self(u+u) == self(u):
                    I.add(u)
            self._idempotents = I
            return I
        else:
            return set(self._idempotents)


    def _relabel_idempotent(self,u):
        if u in self.idempotents():
            return u
        else:
            return u
    def cayley_graph(self, orientation="right",verbose=False):
        r"""
        Return the Cayley graph of the semigroup
        INPUT :

        - ``orientation``  -- string --  value "left", "left_right", "right".

        EXAMPLES::

            sage: from pysemigroup import Automaton, TransitionSemiGroup
            sage: d = {(0, 'a'): [2], (0, 'b'): [1], (1, 'a'): [0],  (1, 'b'): [2], (2, 'a'): [2], (2, 'b'): [2]}
            sage: A = Automaton(d,[1],[1])
            sage: S = TransitionSemiGroup(A)
            sage: G = S.cayley_graph()
            sage: G
            Looped multi-digraph on 6 vertices
        """
        if orientation in self._cayley_graphs:
            return nx.DiGraph(self._cayley_graphs[orientation])
        d = {}
        A = self._generators
        G = nx.DiGraph()
        for x in self:
            G.add_node(x)
        for x in self:
            for c in A:
                if orientation in  ["right","left_right"]:
                    if verbose:
                        print(x,"#",x+c)
                    edge = (x,self.representent(x+c,verbose=verbose))
                    if edge in G.edges():
                        G.edges[(edge[0],edge[1])]["label"].add((c,"r"))
                    else:
                        G.add_edge(x,self(x+c), label=set([(c,"r")]))
                if orientation in  ["left","left_right"]:
                    edge = (x,self(c+x))
                    if edge in G.edges():
                        G.edges[(edge[0],edge[1])]["label"].add((c,"l"))
                    else:
                        G.add_edge(x,self(c+x), label=set([(c,"l")]))
        self._cayley_graphs[orientation] = G
        return nx.DiGraph(G)

    def cayley_graphviz_string(self, orientation="left_right"):
        s = 'digraph {\n node [margin=0 shape="circle" ]\n'
        Elements = set(self)
        while (len(Elements) > 0):
            e = Elements.pop()
            if (orientation == "left"):
                E = self.L_class_of_element(e)
            elif (orientation == "right"):
                E = self.R_class_of_element(e)
            elif (orientation == "left_right"):
                E = self.J_class_of_element(e)
            else:
                raise ValueError("Orientation(={}) should be 'left', 'right', or 'left_right'".format(orientation))
            Elements = Elements.difference(E)
            s = s + 'subgraph cluster'+str(e)+'{style=filled;\ncolor=black;\nfillcolor=azure;\n'
            for x in E:
                if (x == ""):
                    s = s+'  "'+str(x)+'" [label="'+str(1)+'"'
                else:
                    s = s+'  "'+str(x)+'" [label="'+str(x)+'"'

                if x in self.idempotents():
                    s = s + ' fontcolor="red"'
                else:
                    s = s + ' fontcolor="blue"'
                s = s +'];\n'
            s = s + '}\n'
        for x in self:
            for y in self:
                edge = []
                for a in self._generators:
                    if (orientation in ["right","left_right"]) and y == self(x+a):
                        edge.append("."+str(a))
                    if (orientation in ["left","left_right"]) and y == self(a+x):
                        edge.append(str(a)+".")

                if len(edge)>0:
                    s = s+ '  "'+str(x)+'" -> "'+str(y)+'"['
                    if y == x:
                        s =  s+ 'topath="loop above"'
                    s = s + '];\n'
        s = s+'}'
        return s
    def is_Group(self):
        E = self.idempotents()
        return (len(E) == 1)



    def idempotent_power(self,u):
        r"""
        Return the idempotents power of the word u

        INPUT :

        -  ``u`` -  (word)

        EXAMPLES::

            sage: from pysemigroup import Automaton, TransitionSemiGroup
            sage: d = {(0, 'a'): [1], (1, 'a'): [0]}
            sage: A = Automaton(d,[0],[1])
            sage: S = TransitionSemiGroup(A)
            sage: S.idempotent_power('a')
            1
        """
        self._compute_semigroup()
        if isinstance(u,str) or isinstance(u,tuple):
            f = self._Representations_rev[self.representent(u)]
            g = f
            while not (g.is_idempotent()) :
                g =  f * g
            return self._Representations[g]
        else:
            g = u
            while not (g.is_idempotent) :
                g = u * g
            return g
    def _get_J_topological_sort(self):
        try :
            T = self._J_topological_sort
        except:
            G = self.cayley_graph(orientation="left_right")
            K = nx.condensation(G,nx.strongly_connected_components(G))
            T = nx.topological_sort(K)
            T = [K.nodes[i]["members"] for i in T]
            self._J_topological_sort = T
        return T
    def pop_J_maximal(self,E):
        r"""
            Output the J-maximal element of E and delete it.
        INPUT :
        -  ``E`` -  set

        EXAMPLES::
            sage: from pysemigroup import Automaton, TransitionSemiGroup
            sage: d = {(0, 'a'): [2], (0, 'b'): [1], (1, 'a'): [0],  (1, 'b'): [2], (2, 'a'): [2], (2, 'b'): [2]}
            sage: A = Automaton(d,[1],[1])
            sage: S = TransitionSemiGroup(A)
            sage: E = set(S.elements())
            sage: S.pop_J_maximal(E)
            1

        """

        if len(E) == 0 or not E.issubset(self.elements()):
            return False
        T = self._get_J_topological_sort()
        for x in T:
            F = set(x).intersection(E)
            if len(F) > 0:
                y = F.pop()
                E.remove(y)
                return y
        return E.pop()

    def is_sub_semigroup(self, S, verbose = False):
        r"""
        Return whether S is a sub semigroup of self.

        INPUT :

        -  ``S`` - iterable of words
        -  ``verbose`` - boolean


        EXAMPLES::
            sage: from pysemigroup import Automaton, TransitionSemiGroup
            sage: d = {(0, 'a'): [1], (1, 'a'): [0]}
            sage: A = Automaton(d,[0],[1])
            sage: S = TransitionSemiGroup(A)
            sage: S.is_sub_semigroup(("a"))
            False
            sage: S.is_sub_semigroup(("aa"))
            False
            sage: S.is_sub_semigroup((""))
            True
        """
        for i in S:
            for j in S:
                if verbose:
                    print(str((i,j)))
                if not self.representent(i+j) in set(S):
                    return False
        return True

    def sub_semigroup_generated(self, E, monoid=False):
        Generators = set([self(x) for x in E])
        NewElements = set(Generators)
        d = {}
        if monoid:
            for x in Generators:
                d[(self(""),(x,))] = [x]
        Reached = set(Generators)
        while len(NewElements) > 0:
            Buffer = set()
            for x in NewElements:
                for y in Generators:
                    z = self(x+y)
                    d[(x,(y,))] = [z]
                    if (not z in Reached):
                        Buffer.add(z)
                        Reached.add(z)
            NewElements = set(Buffer)
        return TransitionSemiGroup(Automaton(d,[],[]),monoid=monoid)

    def sub_monoid_generated(self,E):
        return self.sub_semigroup_generated(E,monoid=True)

    def _stable(self, verbose=False):
        if self._stable_comp:
            return self._stable_comp
        i = 1

        A = self._generators
        S = set([self.representent(a) for a in A])
        while not self.is_sub_semigroup(S):
            i = i + 1
            if verbose:
                print(i)
                print(S)
            T = set([])
            for u in S:
                for a in A:
                    T.add(self(u+a))
            S = T
        self._stable_comp = (S,i)
        return (S,i)

    def stability_index(self, verbose=False):
        r"""
        Return the stablility index

        OUTPUT :

            integer s

        EXAMPLES::

            sage: from pysemigroup import Automaton, TransitionSemiGroup
            sage: A = RegularLanguage("(a*b)^x",['a','b']) # not tested
            sage: d = {(0, 'a'): [2], (0, 'b'): [1], (1, 'a'): [0], (1, 'b'): [2], (2, 'a'): [2], (2, 'b'): [2]}
            sage: A = Automaton(d,[0],[1])
            sage: S = TransitionSemiGroup(A)
            sage: stability_index()
            2
        """

        return self._stable()[1]

    def stable_set(self, verbose=False):
        r"""
        Return the stable semigroup

        OUTPUT :

            tuple (i,S)

        EXAMPLES::

            sage: from pysemigroup import Automaton, TransitionSemiGroup
            sage: A = RegularLanguage("(a*b)^x",['a','b']) # not tested
            sage: d = {(0, 'a'): [2], (0, 'b'): [1], (1, 'a'): [0], (1, 'b'): [2], (2, 'a'): [2], (2, 'b'): [2]}
            sage: A = Automaton(d,[0],[1])
            sage: S = TransitionSemiGroup(A)
            sage: S.stable_semigroup()
            Transition SemiGroup of Automaton of 4 states

        """
        return self._stable()[0]
    def stable_semigroup(self):
        S = self._stable()[0]
        S.add(monoidElement(""))
        d = {}
        for x in S:
            for y in S:
                d[(str(x),(str(y),))] = [str(self(x+y))]
        return TransitionSemiGroup(Automaton(d,[0],[0]))
    def stabilized_automaton(self):
        r"""
        Return the automata recognizing the language enriched by modular counting
        """
        Aut = self._automaton
        d_old = Aut._transitions
        d = {}
        s = self.stability_index()
        A = Aut._alphabet
        States = Aut._states
        for a in A:
            for i in range(s):
                for x in States:
                    d[(x,(a,i))] = d_old[(x,a)]
        Enr_Aut = Automaton(d,Aut._initial_states,Aut._final_states)
        d = {}
        for a in A:
            for i in range(s):
                d[(i,(a,i))] = [(i+1)%s]
        Count_Aut = Automaton(d,[0],range(s))
        return Enr_Aut.intersection(Count_Aut).minimal_automaton()


    def J_ideal_of_element(self,x):
        xs = self.representent(x)
        if not (x in self._J_ideal):
            self._compute_J_ideal()
        return self._J_ideal[xs]
    def _compute_J_ideal(self):
        G = nx.transitive_closure(self.cayley_graph(orientation="left_right"))
        for x in self:
            self._J_ideal[x] = set(G.neighbors(x))

    def _compute_R_class(self):
        G = self.cayley_graph(orientation="right")
        for K in nx.strongly_connected_components(G):
            for x in K:
                self._R_class[x] = K
    def _compute_L_class(self):
        G = self.cayley_graph(orientation="left")
        for K in nx.strongly_connected_components(G):
            for x in K:
                self._L_class[x] = K
    def _compute_J_class(self):
        G = self.cayley_graph(orientation="left_right")
        for K in nx.strongly_connected_components(G):
            for x in K:
                self._J_class[x] = K
    def J_class_of_element(self,x):
        xs = self.representent(x)
        if xs not in self._J_class:
            self._compute_J_class()
        return set(self._J_class[xs])

    def R_class_of_element(self,x):
        xs = self.representent(x)
        if xs not in self._R_class:
            self._compute_R_class()
        return set(self._R_class[xs])

    def L_class_of_element(self,x):
        xs = self.representent(x)
        if xs not in self._L_class:
            self._compute_L_class()
        return set(self._L_class[xs])

    def H_class_of_element(self,x):
        return set(self.R_class_of_element(x).intersection(self.L_class_of_element(x)))

    def is_element_neutral(self, x):
        E = self.elements()
        for y in E:
            if not self.representent(x+y) == y:
                return False
            if not self.representent(y+x) == y:
                return False
        return True


    def _sg_reverse(self,u,v):
        S = self
        G = S.cayley_graph(orientation="left_right")
        su = S(u)
        sv = S(v)
        l = nx.shortest_path(G,su,sv)
        if l == 0:
            raise TypeError("no inverse")
        l.reverse()
        x = l.pop()
        w1 = monoidElement("")
        w2 = monoidElement("")
        while len(l)>0:
            y = l.pop()
            a = list(G.edges[(x,y)]["label"])[0]
            if a[1] == "r":
                w1 = w1 + a[0]
            else:
                w2 = a[0] + w2

            x = y
        return (w2,w1)

    def newbox(self,x,verbose=False):
        if (x not in self):
            raise TypeError("unboxable element")
        else:
            shiftR = []
            shiftL = []
            box = {}
            Jx = self.J_class_of_element(x)
            Lx = self.L_class_of_element(x)
            Rx = self.R_class_of_element(x)
            Idempotents = Jx.intersection(self.idempotents())
            cx = 0
            cy = 0
            if len(Idempotents) > 0:
                x = Idempotents.pop()
                Lx = self.L_class_of_element(x)
                Rx = self.R_class_of_element(x)
                Idempotents = Idempotents.difference(self.R_class_of_element(x))
                Idempotents = Idempotents.difference(self.L_class_of_element(x))
                while len(Idempotents) > 0:
                    cx = cx + 1
                    cy = cy + 1
                    y = Idempotents.pop()
                    Idempotents = Idempotents.difference(self.R_class_of_element(y))
                    Idempotents = Idempotents.difference(self.L_class_of_element(y))
                    z = self._sg_reverse(x,y)

                    Lx = Lx.difference(self.H_class_of_element(self(z[0]+x)))
                    Rx = Rx.difference(self.H_class_of_element(self(x+z[1])))
                    shiftR.append(z[1])
                    shiftL.append(z[0])
            while len(Lx) > 0:
                l = Lx.pop()
                cy = cy + 1
                Hl = self.H_class_of_element(l)
                Lx = Lx.difference(Hl)
                shiftL.append(self._sg_reverse(x,l)[0])
            while len(Rx) > 0:
                r = Rx.pop()
                cx = cx + 1
                Hr = self.H_class_of_element(r)
                Rx = Rx.difference(Hr)
                shiftR.append(self._sg_reverse(x,r)[1])
        for i in range(0,cx):
            for j in range(0,cy):
                sr = shiftR[i]
                sl = shiftL[j]
                y = self(sl+x+sr)
                box[(i,j)] = (y,self.H_class_of_element(y))
        box["width"] = cx
        box["height"] = cy
        return box

    def newbox_oldbox(self,x,verbose=False):
        if verbose:
            print("Oldbox dealing with:"+str(x))
        nbox = self.newbox(x,verbose=verbose)
        obox = []
        for x in range(nbox["height"]):
            L = []
            for y in range(nbox["width"]):
                H = list(nbox[(y,x)][1])

                L.append(H)
            obox.append(L)
        return obox
    def box_representation(self,verbose=False):
        E = set(self.elements())
        dic = {}
        box = {}
        if verbose:
            print("\t computing Jclass ...")
        while len(E)>0:
            x = self.pop_J_maximal(E)
            if verbose:
                    print("\t Dealing with "+str(x))

            Jclass = self.J_class_of_element(x)
            box[x] = []
            E = E.difference(Jclass)
            if verbose:
                print("\t Done.")
        for x in box:
            box[x] = self.newbox_oldbox(x,verbose=verbose)
        return box
    def box_to_file(self,filename):
        save_graphviz(self.graphviz_string(),filename)
    def view(self, save_to_file=None, arrow=True, verbose=False,unfold=True,extension="svg"):
        view_graphviz(self.graphviz_string(arrow=arrow,
                                           verbose=verbose,
                                           unfold=unfold),
                      save_to_file=save_to_file,
                      extension=extension,
                      verbose=verbose)

    def view_cayley(self, orientation="left_right", save_to_file=None, extension="svg", verbose=False):
        view_graphviz(self.cayley_graphviz_string(orientation=orientation),
                      save_to_file=save_to_file,
                      extension=extension,
                      verbose=verbose)

    def cayley_to_file(self,filename,orientation="left_right"):
        save_graphviz(self.cayley_graphviz_string(orientation=orientation),filename)

    def is_Ap(self,verbose=False):
        for e in self.idempotents():
            H = self.H_class_of_element(e)
            if len(H)>1:
                if verbose:
                    print("A non trivial H-class for element:"+str(e))
                return False
        return True
    def is_J(self,verbose=False):
        for x in self:
            J = self.H_class_of_element(x)
            if len(J)>1:
                if verbose:
                    print("A non trivial J-class for element:"+str(x))
                return False
        return True

    def is_Idempotent(self,verbose=False):
        Sg = set(self.elements())
        E = set(self.idempotents())
        if not (E == Sg) and verbose:
            print(str((Sg.difference(E)).pop()))
        return E == Sg

    def is_Commutative(self,verbose=False):
        for x in self:
            for y in self:
                if not (self(x+y) == self(y+x)):
                    if verbose:
                        print((x,y))
                    return False
        return True

    def is_Jun(self,verbose=False):
        return is_Idempotent(verbose=verbose) and is_Commutative(verbose=verbose)


    def has_a_zero(self):
        G = self.cayley_graph(orientation="left_right")
        G.remove_loops()
        if len(G.sinks()) == 1:
            return True
        else:
            return False
    def get_zero(self):
        G = self.cayley_graph(orientation="left_right")
        G.remove_loops()
        if len(G.sinks()) == 1:
            return G.sinks()[0]
        else:
            raise ValueError("The semigroup must have a zero")
    def dump_to_byte(self,file_name):
        L = [len(self)]
        M = list(self)
        for i in range(len(M)):
            for j in range(len(M)):
                 L.append(M.index(self(M[i]+M[j])))
        f = open(file_name,"wb")
        f.write(bytes(L))
        f.close()
    def load_to_byte(self, file_name):
        f = open(file_name, "rb")
        content = f.read()
        f.close()
        size = content[0]
        content = content[1:]
        t = {}
        for x in range(size):
            for y in range(size):
                t[(x,str(y))] = [content[size*x+y]]
        return TransitionSemiGroup(Automaton(t,[0],[0]))
class BuchiTransitionOmegaSG(TransitionSemiGroup):
    def __init__(self, automaton,compute=False):
        TransitionSemiGroup.__init__(self,automaton,monoid=False, compute=compute)
        self._automaton._type = "buchi"
        self._omega_compute = False
        self._omega_elements = set()

    def omega_elements(self, maxsize= 0, verbose=False):
        if len(self._omega_elements)>0:
            return set(self._omega_elements)
        elements = TransitionSemiGroup.elements(self,maxsize=maxsize,verbose=verbose)
        omega_representation= dict()
        omega_representation_rev = dict()
        set_of_omega_raw = set()
        set_of_omega = set()
        for x in elements:
            y =  monoidElement("("+str(x)+")^w")
            K = self.omega_power_raw(x)
            omega_representation[y] = K
            if K not in set_of_omega_raw:
                omega_representation_rev[K] = y
                set_of_omega.add(y)
                set_of_omega_raw.add(K)
        L = list(set_of_omega)
        for x in elements:
            for z in L:
                Mat = self._Representations_rev[x]
                D = omega_representation[z]
                K = Mat*D
                if K not in set_of_omega_raw:
                    omega_representation[x+z] = K
                    set_of_omega_raw.add(K)
                    omega_representation_rev[K] = x+z
                    set_of_omega.add(x+z)

        self._omega_representation = omega_representation
        self._omega_representation_rev = omega_representation_rev
        self._omega_compute = True
        self._omega_elements = set_of_omega
        return set_of_omega

    def omega_product(self,x,z):
        if not self._omega_compute:
            O = self.omega_elements()
        Mat = self._Representations_rev[self(x)]
        D = self._omega_representation[z]
        K = Mat*D
        return self._omega_representation_rev[K]
    def omega_power(self,x):
        if not self._omega_compute:
            O = self.omega_elements()
        O = self.omega_elements()
        return self._omega_representation_rev[self.omega_power_raw(x)]

    def omega_power_raw(self,x):
        Mat = self._Representations_rev[self.idempotent_power(x)]
        D = Mat.diagonal()
        D.projection({0:"-oo"})
        return Mat*D
    def left_omega_cayley(self,label=True):
        d = []
        A = self._generators
        O = self.omega_elements()
        for a in A:
            for x in O:
                if label:
                    d.append((x,self.omega_product(a,x),a))
                else:
                    d.append((x,self.omega_product(a,x)))
        return d

    def _omega_graphviz_string(self,give_repre=None):
        d = self.left_omega_cayley(label=False)
        G = nx.condensation(nx.DiGraph(d))
        T = nx.topological_sort(G)[0]
        sources = G.nodes[T]["members"]
        d = self.left_omega_cayley()
        J_min = self._get_J_topological_sort()
        J_min = J_min[len(J_min)-1]
        J_min = list(set(give_repre).intersection(J_min))[0]
        s = """
subgraph clusteromega{
style=filled;
color=black;
fillcolor=azure;\n"""
        O = self.omega_elements()
        for x in O:
            s +='  "'+str(x)+'"[label="'+str(x)+'",shape=rectangle,fillpcolor=white];\n'
        s += '}\n'
        for e in d:
            s += '   "'+str(e[0])+'"->"'+str(e[1])+'"[label="'+str(e[2])+'."];\n'
        if not give_repre:
            for e in self.idempotents():
                s += '   '+str(e)+'->"'+str(self.omega_power(e))+'"[color=darkgreen, weight=10];\n'
        else:
            for x in give_repre:
                I = set(self.idempotents()).intersection(self.J_class_of_element(x))
                for e in I:
                    s += '   '+str(x)+":"+str(e)+'->"'+str(self.omega_power(e))+'"[color=darkgreen, weight=10];\n'
        for k in sources:
            s +=  '   '+str(J_min)+' -> "'+str(k)+'"[style=invis];\n'
        return s
    def cayley_graphviz_string(self,orientation="left_right"):
        s = TransitionSemiGroup.cayley_graphviz_string(self,orientation=orientation)
        return s[0:len(s)-1]+self._omega_graphviz_string()+"}"

    def graphviz_string(self,arrow=True,verbose=False,unfold=True,colors_list=False):
        x = TransitionSemiGroup.graphviz_string(self,arrow=arrow,verbose=verbose,unfold=unfold,colors_list=colors_list,get_repr=True)
        repre= x[1]
        s = x[0]
        s = s.replace("node","edge[weight=50];\n node")
        return s[0:len(s)-1]+self._omega_graphviz_string(give_repre=repre)+"}"
    def representent(self,u,verbose=False):
        if verbose:
            print(u)
        return TransitionSemiGroup.representent(self,u)
