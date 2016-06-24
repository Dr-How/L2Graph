# Enumerate (1,0)-graphs.
# Following the algorithm of (Chein 1969)

load("Matrix.sage")

def l1(g): # Recognise graphs of level 1
    return (L1(g) and minEig(M(g))<0)

def Strict(g): # Tell if a L1 graph is strict
    d=g.order()
    M0=M(g)
    for i in range(d):
        ind=range(d)
        ind.remove(i)
        M1=M0[ind,ind]
        meig=minEig(M1)
        if meig==0:
            return false
    return true

def Init(g): #Initialize the graph, set all label 3
    for (u,v) in G.edge_iterator(labels=False):
        g.set_edge_label(u,v,3)
        
def Check(g,dict):
    if l1(g):
        if g.order() not in dict:
            dict[g.order()]=[]
        l=dict[g.order()]
        l.append(g)
        
def dchecked(g,dict):
    if g.order() not in dict:
        return false
    return lchecked(g,dict[g.order()])
    
def lchecked(g,list):
    for h in list:
        if g.is_isomorphic(h,edge_labels=true):
            return true
    return false
    
def isA1(g):
    isom,map=graphs.PathGraph(g.order()).is_isomorphic(g,certify=true)
    if not isom:
        return false,None
    for u,v,l in g.edge_iterator():
        if l!=3:
            return false,None
    return true,map

##
## Preparing level-0 graphs (affine or finite).
##

L0C=[] # affine Coxeter graphs in form of a cycle
for i in range (3,10):
    G=graphs.CycleGraph(i);Init(G)
    _circle_embedding(G,G.vertices(),radius=i)
    G.name('~A'+str(i));L0C.append(G)

L0P=[] # affine and euclidean Coxeter graphs in form of a path
L0T=[] # affine and euclidean Coxeter graphs in form of a tree
for i in range(3,10):
    G=graphs.PathGraph(i);Init(G);G.name('A'+str(i))
    pos={}
    for j in range(i):
        pos[j]=(j,j)
    G.set_pos(pos);L0P.append(G)
for i in range(3,10):
    G=graphs.PathGraph(i);Init(G);G.set_edge_label(0,1,4);G.name('B'+str(i))
    pos={}
    for j in range(i):
        pos[j]=(j,j)
    G.set_pos(pos);L0P.append(G)
for i in range(3,9):
    G=graphs.PathGraph(i);G.add_edge(1,i,3);Init(G);G.name('D'+str(i+1))
    pos={}
    for j in range(i):
        pos[j]=(j,j)
    pos[i]=(0,2)
    G.set_pos(pos);L0T.append(G)
for i in range(5,8):
    G=graphs.PathGraph(i);G.add_edge(2,i,3);Init(G);G.name('E'+str(i+1))
    pos={}
    for j in range(i):
        pos[j]=(j,j)
    pos[i]=(1,3)
    G.set_pos(pos);L0T.append(G)
G=graphs.PathGraph(4);Init(G);G.set_edge_label(1,2,4);G.name('F4');
pos={}
for j in range(4):
    pos[j]=(j,j)
G.set_pos(pos);L0P.append(G)
G=graphs.PathGraph(3);Init(G);G.set_edge_label(1,2,5);G.name('H3');
pos={}
for j in range(4):
    pos[j]=(j,j)
G.set_pos(pos);L0P.append(G)
G=graphs.PathGraph(4);Init(G);G.set_edge_label(2,3,5);G.name('H4');
pos={}
for j in range(4):
    pos[j]=(j,j)
G.set_pos(pos);L0P.append(G)

for i in range(3,9):
    G=graphs.PathGraph(i);Init(G);G.add_edge(i-2,i,3);G.set_edge_label(0,1,4);G.name('~B'+str(i+1));
    pos={}
    for j in range(i):
        pos[j]=(j,j)
    pos[i]=(i-3,i-1)
    G.set_pos(pos);L0T.append(G)
for i in range(3,10):
    G=graphs.PathGraph(i);Init(G);G.set_edge_label(i-2,i-1,4);G.set_edge_label(0,1,4);G.name('~C'+str(i));
    pos={}
    for j in range(i):
        pos[j]=(j,j)
    G.set_pos(pos);L0P.append(G)
for i in range(3,8):
    G=graphs.PathGraph(i);Init(G);G.add_edge(i-2,i+1,3);G.add_edge(1,i,3);G.name('~D'+str(i+2));
    pos={}
    for j in range(i):
        pos[j]=(j,j)
    pos[i]=(0,2);pos[i+1]=(i-3,i-1)
    G.set_pos(pos);L0T.append(G)
G=graphs.PathGraph(5);G.add_edges([(2,5),(5,6)]);Init(G);G.name('~E6');
pos={}
for j in range(5):
    pos[j]=(j,j)
pos[5]=(1,3);pos[6]=(0,4)
G.set_pos(pos);L0T.append(G)
G=graphs.PathGraph(7);G.add_edge(3,7);Init(G);G.name('~E7');
pos={}
for j in range(7):
    pos[j]=(j,j)
pos[7]=(2,4)
G.set_pos(pos);L0T.append(G)
G=graphs.PathGraph(8);G.add_edge(2,8);Init(G);G.name('~E8');
pos={}
for j in range(8):
    pos[j]=(j,j)
pos[8]=(1,3)
G.set_pos(pos);L0T.append(G)
G=graphs.PathGraph(5);Init(G);G.set_edge_label(1,2,4);G.name('~F4');
pos={}
for j in range(5):
    pos[j]=(j,j)
G.set_pos(pos);L0P.append(G)
G=graphs.PathGraph(3);Init(G);G.set_edge_label(0,1,6);G.name('~G2');
pos={}
for j in range(3):
    pos[j]=(j,j)
G.set_pos(pos);L0P.append(G)

for g in L0C:
    print g.name(),sorted(M(g).change_ring(RDF).eigenvalues())[0:2]
print len(L0C)
for g in L0P:
    print g.name(),sorted(M(g).change_ring(RDF).eigenvalues())[0:2]
print len(L0P)
for g in L0T:
    print g.name(),sorted(M(g).change_ring(RDF).eigenvalues())[0:2]
print len(L0T)

##
## List of level-0 graphs prepared.
##

L1C=[] # hyperbolic Coxeter graphs of level 1 in form of a cycle
for G in L0P:
    for l in range(3,7):
        for k in range(3,7):
            g=copy(G)
            i=g.order()
            g.add_edges([(0,i,k),(i-1,i,l)])
            if not lchecked(g,L1C):
                g.name("L1C-"+G.name()+'-'+str(l)+str(k))
                if l1(g):
                    _circle_embedding(g,g.vertices(),radius=i+1)
                    L1C.append(g)

for g in L1C:
    print g.name(),sorted(M(g).change_ring(RDF).eigenvalues())[0:2]
len(L1C)

L1T=[] # hyperbolic Coxeter graphs of level 1 in form of a tree
for G in L0P:
    for l in range(3,7):
        i=G.order()
        for j in range(0,i):
            g=copy(G)
            g.add_edge(i,j,l)
            if not lchecked(g,L1T):
                g.name("L1T-"+G.name()+'-'+str(j)+'-'+str(l))
                pos=g.get_pos(); pos[i]=(pos[j][0],pos[j][1]-1); g.set_pos(pos)
                if l1(g):
                    L1T.append(g)
for G in L0T:
    for l in range(3,7):
        i=G.order()
        for j in range(0,i):
            g=copy(G)
            g.add_edge(i,j,l)
            if not lchecked(g,L1T):
                g.name("L1T-"+G.name()+'-'+str(j)+'-'+str(l))
                pos=g.get_pos(); pos[i]=(pos[j][0],pos[j][1]-1); g.set_pos(pos)
                if l1(g):
                    L1T.append(g)

for g in L1T:
    print g.name(),sorted(M(g).change_ring(RDF).eigenvalues())[0:2]
len(L1T)
                    
Ct=[] # hyperbolic Coxeter graphs of level 1 in form of a tailed cycle

for G in L0C:
    for l in range(3,7):
        i=G.order()
        for j in range(0,i):
            g=copy(G)
            g.add_edge(i,j,l)
            pos=g.get_pos(); pos[i]=(0,0); g.set_pos(pos)
            if not lchecked(g,Ct):
                g.name("Ct-"+G.name()+'-'+str(j)+'-'+str(l))
                if l1(g):
                    Ct.append(g)

for g in Ct:
    print g.name(),sorted(M(g).change_ring(RDF).eigenvalues())[0:2]
len(Ct)
                    

types=[L1C,L1T,Ct]
total=3
for type in types:
    total+=len(type)
print "We found ",total

Level1=[] # Full list
Level1S=[] # Strict list

for t in types:
    for g in t:
        if Strict(g):
            Level1S.append(g)
        else:
            Level1.append(g)

# Special graphs

G=graphs.CompleteGraph(4)
Init(G)
G.set_pos({0:(-2,0),1:(2,0),2:(0,2),3:(0,-2)})
Level1.append(G)

G=graphs.CompleteGraph(4)
G.delete_edge(0,1)
Init(G)
G.set_pos({0:(-2,0),1:(2,0),2:(0,2),3:(0,-2)})
Level1.append(G)

G=graphs.CompleteBipartiteGraph(2,3)
Init(G)
G.set_pos({0:(2,0),1:(-2,0),2:(0,2),3:(0,0),4:(0,-2)})
Level1.append(G)

# test hinges

def Hinge(g,v):
    d=g.order()
    M0=M(g)
    ind=range(d)
    ind.remove(v)
    M1=M0[ind,ind]
    meig=minEig(M1)
    if meig==0:
        return true
    return false

for g in Level1:
    for v in g:
        g.set_vertex(v,Hinge(g,v))
    print g.get_vertices()
