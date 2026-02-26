from sys import stdout
import os
import random
import threading

r"""
Li = Dyck of depth at most i
benchmark  
           L1 2.73 ms (monoid 6)
           L2 71.5 ms (monoid 15)
           L3 1.82  s (monoid 31)
           L4 34.6  s (monoid 56)
           L5 9   min (monoid 92)

"""
log = open("chaines/log2.txt", "w")
c_ex_tight_1 = {(0,"a"):[2],
(0,"b"):[3],
(1,"a"):[0],
(1,"b"):[1],
(2,"a"):[4],
(2,"b"):[1],
(3,"a"):[1],
(3,"b"):[4],
(4,"a"):[2],
(4,"b"):[0]}

c_ex_tight_2 = {(3,"a"):[4],
(3,"b"):[3],
(1,"a"):[4],
(4,"a"):[0],
(0,"a"):[3],
(2,"b"):[1],
(2,"a"):[2],
(0,"b"):[4],
(4,"b"):[3],
(1,"b"):[2]}

c_ex_tight_3 = {}

L={}
u = ""
for i in range(10):
    u = "(a"+u+"b)*"
    L[i+1] = simplify_regex(u)
         
def convert_transition_to_online(d):
    s = ""
    for x in d:
        s += str(x[0])+"."+str(x[1])+"->"
        for y in d[x]:
            s += str(y)+","
        s = s[0:len(s)-1]+";\n"
    return s
def convert_transition_from_online(s):
    d = {}
    dic = s.split(";")
    for x in dic:
        if len(x) >0:
            left = x.split("->")[0].split(".")    
            right = x.split("->")[1].split(",")
            d[(left[0],left[1])] = list(right)
    return d
                

     
def display_chains(Pairs, index=0):
    s = """
digraph {
ranksep=0.5;
"""
    E = set(M)
    while len(E) > 0:
        x = M.pop_J_maximal(E)
        Jx = M.J_class_of_element(x)
        E = E.difference(Jx)
        s += "subgraph cluster"+x+"{\n style=filled;\n color=black;\n fillcolor=azure;\n splines=false\n"
        for y in Jx:
            color = '"blue"'
            if y in M.idempotents():
                color = '"red"'
            s+=  '"'+y+'" [label="'+y+'" fontcolor='+color+'];\n'
        s += '}\n'
    for x in M:
        for y in Pairs[index][0][x]:
            s+='"'+x+'"->"'+y+'";\n'
    s+= "}"
    return  s
    from sage.misc.temporary_file import tmp_filename 
    from sage.misc.viewer import browser
    file_dot = tmp_filename(".",".dot")
    file_gif = tmp_filename(".",".png")
    f = file(file_dot,'w')
    f.write(u)
    f.close()
    f = file(str(i)+".dot",'w')
    f.write(u)
    f.close()

    os.system('dot -Tgif %s -o %s; %s %s  2>/dev/null 1>/dev/null '%(file_dot,file_gif,browser(),file_gif))

def next_pairs(M,Drev,verbose=False,output=stdout):    
    # D doit être un dictionnaire de paires, Drev son reverse
    d = {}
    Pairs = {}
    PairsRev = {}
    Active = []
    
    for x in M:
        Active.append((x,x))
        for a in M._generators:
            d[(x,a)] = [M(a)]
        for e in M.idempotents():
            d[(x,e)] = []
        for z in M:
            e = M.idempotent_power(z)
            for y in Drev[e]:
                u = M(e+y+e)
                if u not in d[(x,e)]:
                    d[(x,e)].append(u)
        Pairs[x] = [x]
        PairsRev[x] = [x]
    Letters = [a for a in M._generators]
    Letters.extend([e for e in M.idempotents()])

    while len(Active) > 0:
        NActive = []
        if verbose:
            print  >> output, "\r\t\tRemaning to compute:"+str(len(Active)),
            output.flush()
        for p in Active:
            for a in Letters:
                for u in d[(p[0],a)]:
                    v = (M(p[0]+a),M(p[1]+u))
                    if not v[1] in Pairs[v[0]]:
                         NActive.append(v)
                         Pairs[v[0]].append(v[1])
                         PairsRev[v[1]].append(v[0])
        Active = list(NActive)
    return (Pairs,PairsRev)



def first_pairs(M):
    d = {}
    for x in M:
        d[x] = list(M)
    return d

def compute_pairs(M,verbose=False,output=stdout):
    i = 0
    Finished = False
    Pairs = {}
    Pairs[0] = (first_pairs(M),first_pairs(M))
    while not Finished:
        if verbose:
            print  >> output, "\t Actual level: "+str(i)+"..."
            output.flush()
        Finished = True
        Pairs[i+1] = next_pairs(M,Pairs[i][1],verbose=verbose,output=output)
        print  >> output, "Done"
        for x in Pairs[i+1][0]:
            if len(Pairs[i+1][0][x]) < len(Pairs[i][0][x]):
                Finished = False
        i += 1
    return Pairs
def test_tighness(M,output=stdout,path=""):
    print  >> output, "start computing pairs ..."
    output.flush()
    Pairs = compute_pairs(M,verbose=True, output=output)
    print  >> output, "Done"
    print >> output, "Pairs:"
    if not (path==""):
        f = file(path+"/Pairs","w")
        for i in Pairs:
            f.write(str(Pairs[i])+"\n")
        f.close()
    Conflit = set()
    for i in range(1,len(Pairs)-1):
        print  >> output, "\ttest level "+str(i)
        D = Pairs[i][1]
        F = Pairs[i+1][0]
        j = 1.0
        for e in M.idempotents():
            print  >> output, "\r\t\t"+str(int(100*j/len(M.idempotents())))+"% ...", 
            output.flush()
            j += 1
            for u in F[e]:
                for v in F[e]:
                    for w in D[e]:
                        if M(u+w+v) not in F[e]:
                            Conflit.add((e,M(u+w+v),i)) 
        print  >> output, "\tDone nb of conlit:"+str(len(Conflit))
    if len(Conflit) == 0:
        return True
    else:
        return Conflit
def test_tighness_randomly(N,size=5,alphabet=["a","b"],output=stdout):

    for i in range(N):
        M = TransitionSemiGroup(random_automaton(size,alphabet))
        print  >> output, "\n"
        print  >> output, M._automaton._transitions
        print  >> output, "Size of monoid :"+str(len(M))
        if len(M) <250:
            if not test_tighness(M)==True:
                return M
            else:
                print  >> output, "Done"
        else:
            print  >> output, "Too big"
def idempotent_set(E):
    F = frozenset(E)
    while not (F == product_set(F,F)):
        F = product_set(E,F)
    return F
def product_set(E,F):
    return frozenset([M(x+y) for x in E for y in F])
def aperiodic_point_like(M):
    Pointlike = [frozenset([x]) for x in M]
    Active = list(Pointlike)
    while len(Active)>0:
        print "\r"+str(len(Active)),        
        stdout.flush()
        T = Active.pop()
        P = idempotent_set(T)
        Q = frozenset([M(x+y) for x in P for y in T])
        R = frozenset(P).union(Q)
        if R not in Pointlike:
            Pointlike.append(R)
            Active.append(R)
        for E in Pointlike:
            P = product_set(E,T)
            if P not in Pointlike:
                Pointlike.append(P)
                Active.append(P)
    print ""
    return Pointlike     
def building_one_data(prefix="",size=5,alphabet=["a","b"]):
    A = random_automaton(size,alphabet)
    M = TransitionSemiGroup(A)    
    s = "chaines/data/"+prefix+"."+stringify_aut(A)
    c_ex = file("chaines/data/c_ex.txt","w")
    if not os.path.exists(s):
        os.makedirs(s)
    else:
        return False
    A.to_svg(s+"/automata")
    A.to_dot(s+"/automata")
    f = file(s+"/transitions.txt","w")
    for x in A._transitions:
        f.write(str(x[0])+"."+str(x[1])+"->"+str(A._transitions[x][0])+";\n")
    f.close()
    log = file(s+"/log","w")
    f = file(s+"/monoid_cayley.tx","w")
    f.write(str(list(M))+"\n")        
    for x in M:
        for a in M._generators:
            f.write(x+"."+a+"->"+M(x+a)+";\n")
    f.close()
    print  >> log, "\n"
    print  >> log, A._transitions
    print  >> log, "Size of monoid :"+str(len(M))
    if len(M) <450:        
        if not test_tighness(M,output=log,path=s)==True:
            print >> c_ex,stringify_aut(A) 
        else:
            print  >> log, "Done"
    else:
        print  >> log, "Too big (size >450)"
        
def stringify_aut(A):
    s = str(A._transitions)
    return s.translate(None,"{('[, :]')}")
def building_data(prefix="",size=5,alphabet=["a","b"]):
    i = 0
    while True:
        i += 1
        print >> stdout, "\r"+prefix+str(i),
        stdout.flush()
        building_one_data(prefix=prefix+str(i),size=size,alphabet=alphabet)

def build_data_multi():
    t1 = threading.Thread(target=building_data, kwargs={"prefix":"abc5","size":5,"alphabet":["a","b","c"]})
    t2 = threading.Thread(target=building_data, kwargs={"prefix":"ab5","size":5,"alphabet":["a","b",]})
    t3 = threading.Thread(target=building_data, kwargs={"prefix":"abc4","size":4,"alphabet":["a","b","c"]})
    t4 = threading.Thread(target=building_data, kwargs={"prefix":"ab4","size":4,"alphabet":["a","b",]})
    t1.start()
    t2.start()
    t3.start()
    t4.start()

#def compare_algo_randomly(N):
#    for i in range(N):
#        M = TransitionSemiGroup(random_automaton(5,["a","b"]))
#        print  >> output, len(M)
#        if (len(M) < 50):
#            compare_algo(M)
#def compare_algo(M):
#    e = compute_pairs_new(M)
#    print  >> output, "pairs one"
#    d = compute_all_pairs(M)
#    print  >> output, "pairs two"
#    for i in range(len(e)):
#        if not (set(e[i]) == set(d[i])):
#            return (False,M)
#    print  >> output, "################"
#    
#    
#    
    #def reachable_pairs(M,E,verbose=False):
#    Reached = [("","")]
#    Todeal = [("","")]
#    while len(Todeal) > 0:
#        x = Todeal.pop()
#        if verbose:
#            print  >> output, len(Todeal)
#        for e in E:
#            a = []
#                
#            a.append(M(x[0]+e[0]))
#            a.append(M(x[1]+e[1]))
#            if tuple(a) not in Reached:
#                Reached.append(tuple(a))
#                Todeal.append(tuple(a))
#    return Reached
#def closure_pairs(M,E,verbose=False):
#    F = list(E)
#    Done = False
#    while not Done:
#        Done = True
#        G = reachable_pairs(M,F,verbose=verbose)
#        if len(G) > len(F):
#            Done = False
#            F = G
#            if verbose:
#                print  >> output, "A new loop, size="+str(len(G))
#    return F
#def next_pairs_old(M,E):
#    F = set([(x,x) for x in M])
#    for e in E:
#        a = M.idempotent_power(e[1])
#        F.add((a,M(a+e[0]+a)))
#    return F

#    
#def compute_all_pairs(M,verbose=False):
#    pairs = {}
#    E = []
#    for x in M:
#        for y in M:
#            E.append((x,y))
#    pairs[0] = list(E)
#    i = 0
#    Done = False
#    while not Done:
#        Done = True
#        nb = len(E)
#        E = reachable_pairs(M,next_pairs_old(M,E),verbose=verbose)
#        if not len(E) == nb:
#            Done = False
#        i+= 1
#        if verbose:
#            print  >> output, "Chaine of level="+str(i)
#        pairs[i] = list(E)
#    return pairs
    
#def compute_pairs_Graph(M,E,verbose=verbose):
#    F = set()
#    for x in M:
#        for y in M:
#            for e in E:
#                F.add(((x,y),(M(x+e[0]),(y+e[1]))))
#    return DiGraph(list(F))
#def compute_all_pairs_Graph(M,verbose=False):
#    pairs = {}
#    E = []
#    for x in M:
#        for y in M:
#            E.append((x,y))
#    pairs[0] = list(E)
#    i = 0
#    Done = False
#    while not Done:
#        Done = True
#        nb = len(E)
#        
#        E = compute_pairs_Graph(M,next_pairs(M,E),verbose=verbose).connected_component_containing_vertex(("",""))
#        if not len(E) == nb:
#            Done = False
#        i+= 1
#        if verbose:
#            print  >> output, "Chaine of level="+str(i)
#        pairs[i] = list(E)
#    return pairs

#def test_tightness(M,d):
#    idempotentLeft = {}
#    idempotentRight = {}
#    for i in d:
#        for e in M.idempotents():
#            idempotentLeft[(e,i)] = []
#            idempotentRight[(e,i)] = []
#        for e in d[i]:
#            if (M(e[0]+e[0]) == e[0]):
#                idempotentLeft[(e[0],i)].append(e[1]) 
#            if (M(e[1]+e[1]) == e[1]):
#                idempotentRight[(e[1],i)].append(e[0]) 
#    for e in M.idempotents():
#        for i in d:
#            for u in idempotentLeft[(e,i)]:
#                for v in idempotentLeft[(e,i)]:
#                    for w in idempotentRight[(e,i)]:
#                            if not M(u+w+v) in idempotentLeft[(e,i)]:
#                                return False
#    return True

