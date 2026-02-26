    

def process_test(n,alphabet,verbose=False):
    S = TransitionSemiGroup(random_automaton(n,alphabet))
    count = 0
    while (test_conj_semigroup(S,verbose=verbose)):
        if verbose:
            print count
        count = count+1
        f = open("status", 'a')
        f.write(str(S._automaton._transitions)+"\n")
        f.write("Taille du monoide:"+str(len(S))+"\n")
        f.write(" test n° "+str(count)+"\n")
        f.write("######################\n")
        f.close()
        S = random_aperiodic_semigroup(n,alphabet,verbose=verbose)
    f = open("status", 'a')
    f.write("Contre-exemple\n")
    f.write(str(L)+"\n")
    f.write(str(count)+"\n")
    f.close()
    
def core_js_chain(data):
    return """ <html lang="en">
<head>
    <meta charset="utf-8" />
    
<style>
.node {
    stroke: darkblue;
  stroke-width: 1.5px;
    cursor: move;
}

.group {
  stroke: #fff;
  stroke-width: 1.5px;
  cursor: move;
  opacity: 0.7;
}

.link {
  stroke: #7a4e4e;
  opacity: 0.5;
  stroke-width: 3px;
  stroke-opacity: 1;
}

.label {
    stroke: white;
    fill: white;
    font-family: Verdana;
    font-size: 25px;
    text-anchor: middle;
    cursor: move;

}

</style>
</head>
<body>

    <script src="http://marvl.infotech.monash.edu/webcola/extern/d3.v3.js"></script>
    <script src="http://marvl.infotech.monash.edu/webcola/cola.v3.min.js"></script>
<script>
    var width = 1560,
        height = 1000;
    """+data+"""
    

    var color = d3.scale.category20();
    var cola = cola.d3adaptor()
        .linkDistance(function(l){return l.length;})
        .size([width, height]);
    var cont = 1;
    var svg = d3.select("body").append("svg")
        .attr("width", width)
        .attr("height", height)
        .on("dblclick",function(){  cont = 0;});
    cola
            .nodes(nodes)
            .links(links)
            .constraints(constraints)
        .start();

        var link = svg.selectAll(".link")
            .data(links)
              .enter().append("line")
            .attr("class", "link")
            .style("stroke-width",function(d){if (d.type=="chain"){ return "8px";}})            
            .style("stroke",function(d){if (d.type=="chain"){ return color(1+d.depth*2);} else {return d.colorlink;}})
            .style("opacity",function(d){ if (d.type=="chain"){ return 1;} else {return d.opacity;}});

        link.append("title").text(function(d){ return d.type;});
        var pad = 3;
        var node = svg.selectAll(".node")
            .data(cola.nodes())
            .enter().append("g")
            .attr("class", "node")   
            .call(cola.drag);
            node.append("circle")   
            .attr("r", 5)
            .style("fill", function (d) { return color(2); });  

        node.append("title")
            .text(function (d) { return d.name; });

        cola.on("tick", function () {
            if (cont == 1){
            cola.start();
            }
            link.attr("x1", function (d) { return d.source.x; })
                .attr("y1", function (d) { return d.source.y; })
                .attr("x2", function (d) { return d.target.x; })
                .attr("y2", function (d) { return d.target.y; });

             node.attr("transform", function(d) { 
  	        return "translate(" + (d.x) + "," + (d.y) + ")"; });
         
            

        });

</script>
</body>
</html>"""
def js_chain_visu(semigroup,start,depth,verbose=False):
    from sage.misc.viewer import browser
    S = semigroup

    chain = []     
    chain.append(compute_pairs_level1(S,verbose=verbose))
    for i in range(depth-1):
        chain.append(compute_pairs(S,chain[i],verbose=verbose))
    Gcal =  DiGraph(S.caley_graph(idempotent=False, loop=False, orientation="left_right"))  
    Gd = DiGraph(Gcal)
    J = {}
    J_rev = {}
    Sp = set(S)
    while len(Sp)>0:
        x = Sp.pop()
        Jx = list(S.J_class_of_element(x))
        x = Jx[0]
        J[x] = Jx
        for y in J[x]:
            J_rev[y] = x
        
        Sp = Sp.difference(J[x])
        Gd.merge_vertices(J[x])
    Gd.allow_loops(False)
    neighbors = set(Gd.neighbors_out(""))
    layer = {"":0}
    i = 0   
    while len(neighbors)>0:
        newN = set()
        i = i+1
        for x in neighbors:
            layer[x] = i
            newN = newN.union(Gd.neighbors_out(x))  
        neighbors = newN        

    if verbose:
        print "done."
    
    A = S._automaton._alphabet
    file_html = sage.misc.temporary_file.tmp_filename(".",".html")
    #file_jsop = sage.misc.temporary_file.tmp_filename(".",".jsop")
    Ls = list(Gcal.vertices())
    
    s = "var nodes = [\n"    
    for x in Ls:
        s = s+'{"name":"'+x+'", "id":"'+str(Ls.index(x))+'"},\n'
    s = s[0:len(s)-2]+'\n];\n var links = [\n'   
    for e in Gcal.edges():
        diff = layer[J_rev[e[1]]]-layer[J_rev[e[0]]]
        if (diff > 0):
            s = s + '{"source":'+str(Ls.index(e[0]))+',"target":'+str(Ls.index(e[1]))+',length:'+str(50+50*diff)+',"chain":"False","colorlink":"blue","opacity":0.3, "type":"'+e[2]+'"},\n'
        else:
            s = s + '{"source":'+str(Ls.index(e[0]))+',"target":'+str(Ls.index(e[1]))+',length:'+str(25+40*diff)+',"colorlink":"red","opacity":1, "type":"'+e[2]+'"},\n'
    for i in range(start,depth):
        for x in chain[i]:
            for y in chain[i][x]:
                if not (x==y):
                    s = s + '{"source":'+str(Ls.index(x))+',"target":'+str(Ls.index(y))+',"depth":'+str(i)+',"colorlink":"green","opacity":'+str(float(0.1*(0.3+2*i/depth)))+', "type":"chain"},\n'

    s = s[0:len(s)-2] +'\n];\n  var constraints = [ \n'
    s = s + '{"type":"alignment","axis":"y","offsets":[  '
    for x in layer:
       s = s + '{"node":"'+str(Ls.index(x))+'","offset":"'+str(layer[x]*100)+'"},'       
    s = s[0:len(s)-1]+ ']},\n'
    s = s[0:len(s)-2] + "\n ];\n"               
    f = file(file_html,'w')
    f.write(core_js_chain(s))
    f.close()
    os.system('%s %s  2>/dev/null 1>/dev/null '%(browser(),file_html))

def compare_set_pairs(self, l,k):
    for x in l:
        if x not in k:
            return False
    for x in k:
        if x not in l:
            return False
    return True

def idempotent_set_pairs(semigroup,l):
    u2 = set()
    u = l.copy()
    ubuff= ()
    while not(compare_set_pairs(semigroup,u,u2)):
        ubuff=set()
        for x in u:
            for y in l:
                ubuff.add((semigroup(x[0]+y[0]), semigroup(x[1]+y[1])))  
        u = ubuff.copy()
        u2 = set()
        for x in u:
            for y in u:
                u2.add((semigroup(x[0]+y[0]), semigroup(x[1]+y[1])))
    

    return u


def compute_pairs_MT(semigroup, verbose=False):   
        if verbose:
            print "starting computation ... "
        A = set([semigroup(a) for a in semigroup._automaton._alphabet])
        P = power_set(A)
        P.remove(frozenset([]))
        if verbose:
            print "first set len:"+str(len(P))
            print P
        F = {}        
        T = {}
        for E in P:
                F[E] = []
                T[E] = []        
        for x in semigroup:
            if verbose:
                print "element :"+x
            for E in P:
                if set([x]) == E:
                    F[E].append(set([(x,x)]))
                if x in semigroup.sub_semigroup_generated(E):
                    T[E].append(("",x))        
             
        complete = False   
        while not complete:
            complete = True
            for E in P:
                if verbose:
                    print E
                    print len(F[E])
                    print "%%%%%%%%%%%%%%%%"    
                K = []      
                count = 0  
                for S in F[E]:
                    if verbose:
                        print str(count)+" on "+str(len(F[E]))
                        count = count +1
                    Sw = idempotent_set_pairs(semigroup,S)
                    if verbose:
                        print len(F)
                    SwTSw = set([])
                    for x in Sw:
                        for y in T[E]:
                            for z in Sw:    
                                SwTSw.add((semigroup(x[0]+y[0]+z[0]), semigroup(x[1]+y[1]+z[1])))                     
                    if SwTSw not in F[E] and SwTSw not in K:  
                        complete = False
                        K.append(SwTSw)
                F[E] = F[E]+K                          
                for G in P:
                    K = []  
                    count = 0        
                    for S1 in F[E]:
                        for S2 in F[G]:
                            if verbose:
                                print str(count)+" "+str(len(F[E])*len(F[G]))
                                count = count +1
                            
                            S1S2 = set([])                
                            for x in S1:
                                for y in S2:
                                    S1S2.add((semigroup(x[0]+y[0]), semigroup(x[1]+y[1])))                           
                            if S1S2 not in F[E.union(G)] and S1S2 not in K:
                                complete = False
                                K.append(S1S2)                                
                    F[E.union(G)] = F[E.union(G)]+K
                              
        F2 = {}
        for x in semigroup.elements():
            F2[x] = [x]
        for E in P:
            for X in F[E]:
                if len(X) > -1:
                    for y in X:
                        if y[1] not in F2[y[0]]: 
                           F2[y[0]].append(y[1])  
                         
        return F2       
    
