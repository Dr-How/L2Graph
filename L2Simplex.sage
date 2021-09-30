import sage.graphs.graph_plot

def M(g): # Coxeter matrix from Coxeter graph
    d=g.order()
    m=matrix(SR,d)
    for i in range(d):
        m[i,i]=1 # Diagonal entries
    for u,v,l in g.edge_iterator():
        m[u,v]=-cos(pi/l) # finite order
        m[v,u]=m[u,v] # infinite order
    return m
    
def minEig(M): # Minimum eigenvalue of a matrix
    m=copy(M).change_ring(RDF)
    mineig=min(m.eigenvalues())
    if abs(mineig)<0.001: # Take error into account
        return 0
    return mineig

def L1(g): # Recognise graphs of level <=1
    d=g.order()
    M0=M(g)
    for i in range(d):
        ind=list(range(d))
        ind.remove(i) # delete one vertex
        M1=M0[ind,ind]
        if minEig(M1)<0: # has negative eigenvalue, not positive semi-definite
            return false
    return true

def l1(g): # Recognise graphs of level 1
    return (L1(g) and minEig(M(g))<0)

def L2(g): # Recognise graphs of level 2
    d=g.order()
    if L1(g):
        return false
    M0=M(g)
    for i in range(d):
        for j in range(i+1,d):
            ind=list(range(d))
            ind.remove(i) # delete two vertices
            ind.remove(j)
            M2=M0[ind,ind]
            if minEig(M2)<0: # has negative eigenvalue
                return false
    return true

def Strict(g): # Tell if a L2 graph is strict
    d=g.order()
    M0=M(g)
    for i in range(d):
        for j in range(i+1,d):
            ind=list(range(d))
            ind.remove(i)
            ind.remove(j)
            M2=M0[ind,ind]
            meig=minEig(M2)
            if meig==0: # has eigenvalue 0, not positive definite
                return false
    return true
    
def Init(g): #Initialize the graph, set all label 3
    for (u,v) in G.edge_iterator(labels=False):
        g.set_edge_label(u,v,3)
        
def Check(g,dict): # Check the graph g; if level 2, add it into the dictionary
    if L2(g):
        if g.order() not in dict:
            dict[g.order()]=[]
        list=dict[g.order()]
        list.append(g)
        
def dchecked(g,dict): # If g was already in the dictionary
    if g.order() not in dict:
        return false
    return lchecked(g,dict[g.order()])
    
def lchecked(g,list): # If g was already in the list
    for h in list:
        if g.is_isomorphic(h,edge_labels=true):
            return true
    return false
    
def coloring(g): # Distinguish real vertices by colors
    color_dict={'black':[],'lightgray':[],'white':[]}
    B=M(g).change_ring(RDF)
    W=B.inverse()
    for v in range(g.order()):
        norm=W.row(v)*B*W.column(v)
        if abs(norm-1)<0.0001: # Surreal vertices
            color_dict['lightgray'].append(v)
        elif norm>0.0001: # Real vertices
            color_dict['white'].append(v)
        else: #Imaginary vertices
            color_dict['black'].append(v)
    return color_dict

def output(type,s): # Output images
    for i in type:
        images=[]
        graphstr=""
        for g in type[i]:
            h=copy(g)
            for u,v,l in g.edges(): # Remove label 3 on the edges
                if l==3:
                    h.set_edge_label(u,v,'')
            p=h.plot(vertex_labels=false,vertex_colors=coloring(g),vertex_size=30,edge_labels=true,edge_color='red',graph_border=Strict(g),edge_labels_background='transparent')
            images.append(p)
            graphstr=graphstr+h.graphviz_string(edge_labels=True)+"\n"
        n=floor(30/i) # Control the size of the figure array
        m=ceil(len(type[i])/n)
        Image=graphics_array(images,m,n) # Create an array of figures
        filename=s+'-'+str(g.order())+".eps" # Output eps file
        if type==T:
            Image.save(filename,figsize=[8,i*m/4])
        else:
            Image.save(filename,figsize=[6,i*m/5])
        # Image.show()
        print(filename,"output")
        filename=s+'-'+str(g.order())+".txt" # Output txt file
        with open(filename, "w") as f:
            f.write(graphstr)
            f.close()

print("George Maxwell found "+str(186+66+36+13+10+8+4)) # For comparison

L0C=[] # affine Coxeter graphs in form of a cycle
for i in range (3,10): # type ~A
    G=graphs.CycleGraph(i);Init(G)
    pos=G.layout_circular()
    for jj in range(len(pos)):
        pos[jj] = (pos[jj][1]*i, -pos[jj][0]*i)
    G.set_pos(pos)
    G.name('~A'+str(i));L0C.append(G)

L0P=[] # affine and euclidean Coxeter graphs in form of a path
L0T=[] # affine and euclidean Coxeter graphs in form of a tree

# Manually input finite graphs
for i in range(3,10): # type A
    G=graphs.PathGraph(i);Init(G);G.name('A'+str(i))
    pos={}
    for j in range(i):
        pos[j]=(j,j)
    G.set_pos(pos);L0P.append(G)
for i in range(3,10): # type B
    G=graphs.PathGraph(i);Init(G);G.set_edge_label(0,1,4);G.name('B'+str(i))
    pos={}
    for j in range(i):
        pos[j]=(j,j)
    G.set_pos(pos);L0P.append(G)
for i in range(3,9): # type D
    G=graphs.PathGraph(i);G.add_edge(1,i,3);Init(G);G.name('D'+str(i+1))
    pos={}
    for j in range(i):
        pos[j]=(j,j)
    pos[i]=(0,2)
    G.set_pos(pos);L0T.append(G)
for i in range(5,8): # type E
    G=graphs.PathGraph(i);G.add_edge(2,i,3);Init(G);G.name('E'+str(i+1))
    pos={}
    for j in range(i):
        pos[j]=(j,j)
    pos[i]=(1,3)
    G.set_pos(pos);L0T.append(G)
G=graphs.PathGraph(4);Init(G);G.set_edge_label(1,2,4);G.name('F4'); # type F4
pos={}
for j in range(4):
    pos[j]=(j,j)
G.set_pos(pos);L0P.append(G)
G=graphs.PathGraph(3);Init(G);G.set_edge_label(1,2,5);G.name('H3'); # type H3
pos={}
for j in range(4):
    pos[j]=(j,j)
G.set_pos(pos);L0P.append(G)
G=graphs.PathGraph(4);Init(G);G.set_edge_label(2,3,5);G.name('H4'); # type H4
pos={}
for j in range(4):
    pos[j]=(j,j)
G.set_pos(pos);L0P.append(G)

# Manually input affine Coxeter 
for i in range(3,9): # type ~B
    G=graphs.PathGraph(i);Init(G);G.add_edge(i-2,i,3);G.set_edge_label(0,1,4);G.name('~B'+str(i+1));
    pos={}
    for j in range(i):
        pos[j]=(j,j)
    pos[i]=(i-3,i-1)
    G.set_pos(pos);L0T.append(G)
for i in range(3,10): # type ~C
    G=graphs.PathGraph(i);Init(G);G.set_edge_label(i-2,i-1,4);G.set_edge_label(0,1,4);G.name('~C'+str(i));
    pos={}
    for j in range(i):
        pos[j]=(j,j)
    G.set_pos(pos);L0P.append(G)
for i in range(3,8): # type ~D
    G=graphs.PathGraph(i);Init(G);G.add_edge(i-2,i+1,3);G.add_edge(1,i,3);G.name('~D'+str(i+2));
    pos={}
    for j in range(i):
        pos[j]=(j,j)
    pos[i]=(0,2);pos[i+1]=(i-1,i-3)
    G.set_pos(pos);L0T.append(G)
G=graphs.PathGraph(5);G.add_edges([(2,5),(5,6)]);Init(G);G.name('~E6'); # type ~E6
pos={}
for j in range(5):
    pos[j]=(j,j)
pos[5]=(1,3);pos[6]=(0,4)
G.set_pos(pos);L0T.append(G)
G=graphs.PathGraph(7);G.add_edge(3,7);Init(G);G.name('~E7'); # type ~E7
pos={}
for j in range(7):
    pos[j]=(j,j)
pos[7]=(2,4)
G.set_pos(pos);L0T.append(G)
G=graphs.PathGraph(8);G.add_edge(2,8);Init(G);G.name('~E8'); # type ~E8
pos={}
for j in range(8):
    pos[j]=(j,j)
pos[8]=(1,3)
G.set_pos(pos);L0T.append(G)
G=graphs.PathGraph(5);Init(G);G.set_edge_label(1,2,4);G.name('~F4'); # type ~F4
pos={}
for j in range(5):
    pos[j]=(j,j)
G.set_pos(pos);L0P.append(G)
G=graphs.PathGraph(3);Init(G);G.set_edge_label(0,1,6);G.name('~G2'); # type ~G2
pos={}
for j in range(3):
    pos[j]=(j,j)
G.set_pos(pos);L0P.append(G)

# Verify the eigenvalues of level 0 graphs
for g in L0C:
    print(g.name(),sorted(M(g).change_ring(RDF).eigenvalues())[0:2])
print(len(L0C))
for g in L0P:
    print(g.name(),sorted(M(g).change_ring(RDF).eigenvalues())[0:2])
print(len(L0P))
for g in L0T:
    print(g.name(),sorted(M(g).change_ring(RDF).eigenvalues())[0:2])
print(len(L0T))

L1C=[] # hyperbolic Coxeter graphs of level 1 in the form of a cycle
for G in L0P: # Construct from a path
    for l in range(3,7): # Take two labels
        for k in range(3,7):
            g=copy(G)
            i=g.order()
            g.add_edges([(0,i,k),(i-1,i,l)]) # Connect a vertex to the two ends
            if not lchecked(g,L1C):
                g.name("L1C-"+G.name()+'-'+str(l)+str(k))
                if l1(g):
                    pos=g.layout_circular()
                    for jj in range(len(pos)):
                        pos[jj] = (pos[jj][1]*(i+1), -pos[jj][0]*(i+1))
                    g.set_pos(pos)
                    L1C.append(g)

for g in L1C: # Verify eigenvalues
    print(g.name(),sorted(M(g).change_ring(RDF).eigenvalues())[0:2])
len(L1C)

L1T=[] # hyperbolic Coxeter graphs of level 1 in form of a cycle
for G in L0P: # Construct from a path
    for l in range(3,7): # label
        i=G.order()
        for j in range(0,i):
            g=copy(G)
            g.add_edge(i,j,l) # Attach an edge
            if not lchecked(g,L1T):
                g.name("L1T-"+G.name()+'-'+str(j)+'-'+str(l))
                pos=g.get_pos(); pos[i]=(pos[j][0],pos[j][1]-1); g.set_pos(pos)
                if l1(g):
                    L1T.append(g)
for G in L0T: # Construct from a tree
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

for g in L1T: # Verify eigenvalues
    print(g.name(),sorted(M(g).change_ring(RDF).eigenvalues())[0:2])
len(L1T)

C={} # Coxeter graphs of level 2 in form of a cycle
T={} # Coxeter graphs of level 2 in form of a tree
for G in L1T: # Construct from a tree
    i=G.order()
    for l in range(3,7):
        for j in range(0,i):
            g=copy(G)
            g.add_edge(i,j,l) # Attach an edge and obtain another tree
            if not dchecked(g,T):
                g.name("T-"+G.name()+'-'+str(j)+'-'+str(l))
                pos=g.get_pos() 
                if (pos[j][0]+1,pos[j][1]) in pos.values():
                    pos[i]=(pos[j][0]-1,pos[j][1])
                else:
                    pos[i]=(pos[j][0]+1,pos[j][1])
                g.set_pos(pos)
                Check(g,T)
                    
    isom,map=G.is_isomorphic(graphs.PathGraph(i),certificate=true)
    if isom: # If the tree is actually a path
        for l in range(3,7):
            for k in range(3,7):
                g=copy(G)
                g.relabel(map)
                g.add_edges([(0,i,k),(i-1,i,l)]) # Create a cycle
                if not dchecked(g,C):
                    g.name("C-"+G.name()+'-'+str(l)+str(k))
                    pos=g.layout_circular()
                    for jj in range(len(pos)):
                        pos[jj] = (pos[jj][1]*(i+1), -pos[jj][0]*(i+1))
                    g.set_pos(pos)
                    Check(g,C)
                    
output(C,"C")
output(T,"T")

Ct={} # hyperbolic Coxeter graphs of level 2 in form of a tailed cycle

for G in L1C: # Construct from a cycle
    for l in range(3,7):
        i=G.order()
        for j in range(0,i):
            g=copy(G)
            g.add_edge(i,j,l) # Attach an edge
            pos=g.get_pos(); pos[i]=(0,0); g.set_pos(pos)
            if not dchecked(g,Ct):
                g.name("Ct-"+G.name()+'-'+str(j)+'-'+str(l))
                Check(g,Ct)
                                
for G in L1T: # Construct from a tree
    i=G.order()
    H=G.degree_histogram()
    if H[1]==3: # If there are three leaves
        deg=G.degree()
        center=deg.index(3) # Branching point
        leaves=[]
        for j in range(i): # Find all the leaves
            if deg[j]==1:
                leaves.append(j)
        for j in range(3):
            if G.has_edge(leaves[j],center): # This leaf is connected to the center
                for l in range(3,7):
                    for k in range(3,7):
                        g=copy(G)
                        g.add_edges([(leaves[(j+1)%3],i,l),(leaves[(j+2)%3],i,k)]) # Attach the other two leaves to a new vertex
                        if not dchecked(g,Ct):
                            h=graphs.CycleGraph(i); h.add_edge(0,i)
                            isom,map=g.is_isomorphic(h,certificate=true)
                            g.relabel(map)
                            pos=g.layout_circular()
                            for jj in range(len(pos)):
                                pos[jj] = (pos[jj][1]*(i+1), -pos[jj][0]*(i+1))
                            g.set_pos(pos)
                            pos=g.get_pos(); pos[i]=(0,0); g.set_pos(pos)
                            g.name("Ct-"+G.name()+'-'+str(leaves[(j+1)%3])+str(leaves[(j+2)%3])+'-'+str(l)+str(k))
                            Check(g,Ct)
    if H[1]==2: # There are two leaves, i.e. a path
        isom,map=G.is_isomorphic(graphs.PathGraph(i),certificate=true) # Number the vertices in order
        for l in range(3,7):
            for k in range(3,7):
                g=copy(G)
                g.relabel(map)
                g.add_edges([(1,i,k),(i-1,i,l)]) # Connect the second and the last to a new vertex
                if not dchecked(g,Ct):
                    h=graphs.CycleGraph(i); h.add_edge(0,i)
                    isom,map1=g.is_isomorphic(h,certificate=true)
                    g.relabel(map1)
                    pos=g.layout_circular()
                    for jj in range(len(pos)):
                        pos[jj] = (pos[jj][1]*(i+1), -pos[jj][0]*(i+1))
                    g.set_pos(pos)
                    pos=g.get_pos(); pos[i]=(0,0); g.set_pos(pos)
                    g.name("Ct-sigma"+G.name()+'-'+str(l)+str(k))
                    Check(g,Ct)
                g=copy(G)
                g.relabel(map)
                g.add_edges([(0,i,k),(i-2,i,l)]) # Connect the first and the second last to a new vertex
                if not dchecked(g,Ct):
                    h=graphs.CycleGraph(i); h.add_edge(0,i)
                    isom,map1=g.is_isomorphic(h,certificate=true)
                    g.relabel(map1)
                    pos=g.layout_circular()
                    for jj in range(len(pos)):
                        pos[jj] = (pos[jj][1]*(i+1), -pos[jj][0]*(i+1))
                    g.set_pos(pos)
                    pos=g.get_pos(); pos[i]=(0,0); g.set_pos(pos)
                    g.name("Ct\'-"+G.name()+'-'+str(l)+str(k))
                    Check(g,Ct)
                                     
output(Ct,"Ct")

K4={} # Special graphs from K4
G=graphs.CompleteGraph(4)
Init(G)
pos={0:(-2,0),1:(2,0),2:(-1,2),3:(1,2),4:(0,-1)}
for l in range(3,7):
    g=copy(G)
    g.name("K4+1-"+str(l))
    g.add_edge(0,4,l)
    g.set_pos(pos)
    Check(g,K4)
for i in range(2,5):
    g=copy(G)
    g.name("K4+"+str(i))
    g.add_edges([(j,4,3) for j in range(i)])
    g.set_pos(pos)
    Check(g,K4)
    
output(K4,"K4")

K4dK2={} # Special graphs from K4\K2
G=graphs.CompleteGraph(4)
G.delete_edge(0,1)
Init(G)
Pos={0:(0,2),1:(0,-2),2:(1,0),3:(-1,0)}
for i in range(3,7): # Attach new vertex to one vertex
    g=copy(G)
    g.name("K4/K2+A-"+str(i)) # 'A' represents vertex 0 or 1
    g.add_edge(0,4,i)
    pos=copy(Pos);pos[4]=(2,2);g.set_pos(pos)
    Check(g,K4dK2)
    g=copy(G)
    g.name("K4/K2+B-"+str(i)) # 'B' represents vertex 2 or 3
    g.add_edge(2,4,i)
    pos=copy(Pos);pos[4]=(2,2);g.set_pos(pos)
    Check(g,K4dK2)
for l in range(3,7): # Attach new vertex to two vertices
    for k in range(3,5):
        g=copy(G)
        g.add_edges([(1,4,l),(0,4,k)])
        g.name("K4/K2+AA-"+str(l)+str(k))
        if not dchecked(g,K4dK2):
            pos=copy(Pos);pos[4]=(2,0);g.set_pos(pos)
            Check(g,K4dK2)
        g=copy(G)
        g.add_edges([(2,4,l),(3,4,k)])
        g.name("K4/K2+BB-"+str(l)+str(k))
        if not dchecked(g,K4dK2):
            pos=copy(Pos);pos[4]=(0,1);g.set_pos(pos)
            Check(g,K4dK2)

g=copy(G);g.name("K4/K2+AB");g.add_edges([(0,4,3),(2,4,3)])
pos=copy(Pos);pos[4]=(2,2);g.set_pos(pos);Check(g,K4dK2)
g=copy(G);g.name("K4/K2+AAB");g.add_edges([(0,4,3),(1,4,3),(2,4,3)]) # Attach new vertex to three vertices
pos=copy(Pos);pos[4]=(2,0);g.set_pos(pos);Check(g,K4dK2)
# The following case is not possible
# g=copy(G);g.name("K4/K2+ABB");g.add_edges([(1,4,3),(2,4,3),(3,4,3)]);Check(g,K4dK2)

output(K4dK2,"K4-K2")

K23={} # Special graph group C, built from K23
G=graphs.CompleteBipartiteGraph(2,3)
Pos={0:(2,0),1:(-2,0),2:(0,2),3:(0,0),4:(0,-2)}
Init(G)
for i in range(3,7): # New vertex attach to one vertex
    g=copy(G)
    g.name("K23+A-"+str(i)) # 'A' represents the part with two vertices
    g.add_edge(0,5,i)
    pos=copy(Pos);pos[5]=(2,2);g.set_pos(pos)
    Check(g,K23)
    g=copy(G)
    g.name("K23+B-"+str(i)) # 'B' represents the part with three vertices
    g.add_edge(2,5,i)
    pos=copy(Pos);pos[5]=(2,2);g.set_pos(pos)
    Check(g,K23)

for l in range(3,7): # New vertex attach to two vertices
    for k in range(3,5):
        g=copy(G)
        g.add_edges([(0,5,3),(1,5,3)])
        g.name("K23+AA-"+str(l)+str(k))
        if not dchecked(g,K23):
            pos=copy(Pos);pos[5]=(0,1);g.set_pos(pos)
            Check(g,K23)
        g=copy(G)
        g.add_edges([(2,5,3),(4,5,3)])
        g.name("K23+BB-"+str(l)+str(k))
        if not dchecked(g,K23):
            pos=copy(Pos);pos[5]=(1/2,1);g.set_pos(pos)
            Check(g,K23)
            
g=copy(G);g.name("K23+AB");g.add_edges([(1,5,3),(2,5,3)]);Check(g,K23)

for l in range(3,7): # New vertex attach to three vertices
    for m in range(3,5):
        for n in range(3,7):
            g=copy(G)
            g.add_edges([(2,5,l),(3,5,m),(4,5,n)])
            g.name("K23+BBB-"+str(l)+str(m)+str(n))
            if not dchecked(g,K23):
                pos=copy(Pos);pos[5]=(1/2,1);g.set_pos(pos)
                Check(g,K23)
            g=copy(G)
            g.add_edges([(0,5,l),(1,5,m),(3,5,n)])
            g.name("K23+AAB-"+str(l)+str(m)+str(n))
            if not dchecked(g,K23):
                Check(g,K23)
            g=copy(G)
            g.add_edges([(0,5,l),(2,5,m),(3,5,n)])
            g.name("K23+ABB-"+str(l)+str(m)+str(n))
            if not dchecked(g,K23):
                Check(g,K23)

output(K23,"K23")

Ckt2={} # Cycles with a tail of length 2
for G in L0C: # Construct from cycle
    i=G.order()
    for l in range(3,7):
        for k in range(3,7):
            g=copy(G)
            g.add_edges([(0,i,l),(i,i+1,k)]) # attach the tail
            pos=g.get_pos(); pos[i]=(2*i,0); pos[i+1]=(3*i,0); g.set_pos(pos)
            g.name('C'+str(i)+'t2-'+str(l)+str(k))
            if not dchecked(g,Ckt2):
                Check(g,Ckt2)
                
output(Ckt2,"Ct2")

CC={} # Graphs with two cycles
g=Graph([(0,1,3),(1,2,3),(2,0,3),(0,3,3),(3,4,3),(4,0,3)],name='Butterfly'); # Butterfly graph
g.set_pos({0:(0,0),1:(0,1),2:(1,0),3:(0,-1),4:(-1,0)}); Check(g,CC)

G=graphs.CompleteGraph(4) # From kite graphs
G.delete_edge(2,3)
Init(G)
Pos={0:(-3,0),1:(3,0),2:(0,3),3:(0,-3)}
# In the following 'CCabc' means that on the three paths, there are a,b,c vertices respectively

g=copy(G);g.subdivide_edge(0,1,2);g.subdivide_edges([(0,2),(0,3)],1);g.name('CC222')
pos=copy(Pos);pos[4]=(-1,0);pos[5]=(1,0);pos[6]=(-1,3);pos[2]=(1,3);pos[7]=(-1,-3);pos[3]=(1,-3);g.set_pos(pos); Check(g,CC)
g=copy(G);g.subdivide_edges([(0,1),(0,2),(0,3)],1);g.name('CC122')
pos=copy(Pos);pos[4]=(0,0);pos[5]=(-1,3);pos[2]=(1,3);pos[6]=[-1,-3];pos[3]=(1,-3);g.set_pos(pos); Check(g,CC)
g=copy(G);g.subdivide_edges([(0,2),(0,3)],1);g.name('CC022')
pos=copy(Pos);pos[4]=(-1,3);pos[2]=(1,3);pos[5]=(-1,-3);pos[3]=(1,-3);g.set_pos(pos);Check(g,CC)

g=copy(G);g.subdivide_edges([(0,1),(0,2)],1)
pos=copy(Pos);pos[4]=(0,0);pos[5]=(-1,3);pos[2]=(1,3);g.set_pos(pos);Check(g,CC)
for l in range(3,7):
    for n in range(l,7):
        for m in range (3,7):
            h=copy(g)
            h.set_edge_label(0,5,l)
            h.set_edge_label(5,2,m)
            h.set_edge_label(2,1,n)
            h.name("CC112-"+str(l)+str(m)+str(n))
            if not dchecked(h,CC):
                Check(h,CC)

g=copy(G);g.subdivide_edge(0,2,1)
pos=copy(Pos);pos[4]=(-1,3);pos[2]=(1,3);g.set_pos(pos)
for l in range(3,7):
    for n in range(l,7):
        for m in range (3,7):
            h=copy(g)
            h.set_edge_label(0,4,l)
            h.set_edge_label(4,2,m)
            h.set_edge_label(2,1,n)
            h.name("CC012-"+str(l)+str(m)+str(n))
            if not dchecked(h,CC):
                Check(h,CC)

g=copy(G);g.subdivide_edge(0,1,1)
pos=copy(Pos);pos[4]=(0,0);g.set_pos(pos)
for l04 in range(3,5):
    for l14 in range(3,5):
        for l02 in range(3,5):
            for l03 in range(3,7):
                for l12 in range(3,7):
                    for l13 in range(3,5):
                        h=copy(g)
                        h.set_edge_label(0,2,l02)
                        h.set_edge_label(0,3,l03)
                        h.set_edge_label(0,4,l04)
                        h.set_edge_label(1,2,l12)
                        h.set_edge_label(1,3,l13)
                        h.set_edge_label(1,4,l14)
                        h.name("CC111-"+str(l02)+str(l03)+str(l04)+str(l12)+str(l13)+str(l14))
                        if not dchecked(h,CC):
                            Check(h,CC)

output(CC,"CC")

Cktt={} # Graphs with one cycle and two tails
for G in L0C: # Construct from cycle
    i=G.order()
    for j in range(0,floor(i/2+1)):
        for l in range(3,7):
            for k in range(3,7):
                g=copy(G)
                g.add_edge(0,i,l)
                g.add_edge(j,i+1,k)
                pos=g.get_pos(); pos[i]=(2*i,0); pos[i+1]=(0,0); g.set_pos(pos)
                g.name('C'+str(i)+'t-'+str(l)+str(k))
                if not dchecked(g,Cktt):
                    Check(g,Cktt)
              
output(Cktt,"Ctt")

# Conclusion
types=[K4,K4dK2,K23,C,T,Ct,Cktt,Ckt2,CC]
total=0
for i in range(5,12):
    number=0
    for type in types:
        if i in type:
            number+=len(type[i])
    total+=number
print("We found ",total)

Level2=[] # Output data for the use of other programs
Level2S=[]
for type in types:
    for i in type:
        for g in type[i]:
            if Strict(g):
                Level2S.append(g)
            else:
                Level2.append(g)

###########################################################################
# The following implement the in Proposition 3.1 & 3.2 of (Maxwell, 1998) #
# which find a subgroup of finite index, as a double-check.               #
###########################################################################

def isA(g): # Tell if the graph g is a A_n graph.  If yes, number the vertices in order.
    isom,map=graphs.PathGraph(g.order()).is_isomorphic(g,certificate=true)
    if not isom:
        return false,None
    for u,v,l in g.edge_iterator():
        if l!=3:
            return false,None
    return true,map

def subgroup(g): #subgroups
    sub=[]
    for u,v,l in g.edges():
        g1=copy(g)
        g1.delete_edge(u,v)
        if (l==4 or l==6) and not g1.is_connected():
            gu=g.subgraph(g1.connected_component_containing_vertex(u))
            gv=g.subgraph(g1.connected_component_containing_vertex(v))
            h=copy(g)
            vuv=h.add_vertex()
            uvu=h.add_vertex()
            uinA1,mapu=isA(gu)
            if gu.degree(u)<=1 and uinA1: 
                if l==6:
                    h.add_edge(uvu,v,3)   
                for w in gv.neighbor_iterator(v):
                    h.add_edge(uvu,w,g.edge_label(v,w))
                if gu.degree(u)==1:
                    #for i in range(1):#maxwell1982
                    for i in range(gu.order()-1): # maxwell, 1998
                        if mapu[0]==u:
                            j=i+1
                        else:
                            j=gu.order()-i-2
                        h1=copy(h)
                        h1.add_edge(uvu,mapu[j],l) 
                        h1.delete_vertices([u,vuv])
                        sub.append(h1) 
                else:
                    h1=copy(h); h1.delete_vertices([u,vuv])
                    sub.append(h1);                 
            vinA1,mapv=isA(gv)
            if gv.degree(v)<=1 and vinA1:
                if l==6:
                    h.add_edge(vuv,u,3)
                for w in gu.neighbor_iterator(u):
                    h.add_edge(vuv,w,g.edge_label(u,w))
                if gv.degree(v)==1:
                    #for i in range(1):#maxwell1982
                    for i in range(gv.order()-1): # maxwell, 1998
                        if mapv[0]==v:
                            j=i+1
                        else:
                            j=gv.order()-i-2
                        h2=copy(h)
                        h2.add_edge(vuv,mapv[j],l)
                        h2.delete_vertices([v,uvu]);
                        sub.append(h2)
                else:
                    h2=copy(h); h2.delete_vertices([v,uvu])
                    sub.append(h2)
            if h.degree(uvu)>0 and h.degree(vuv)>0 and l==6:
                h3=copy(h)
                if gu.degree(u)==1:
                    h3.add_edge(gu.neighbors(u)[0],uvu,l)
                if gv.degree(v)==1:
                    h3.add_edge(gv.neighbors(v)[0],vuv,l)
                h3.delete_vertices([u,v])
                sub.append(h3)     
    return sub

def maxgroup(list):
    maxl=copy(list)
    for g in list:
        subg=subgroup(g)
        for h in subg:
            i=0
            while i<len(maxl) and not h.is_isomorphic(maxl[i],edge_labels=true):
                i=i+1
            if i<len(maxl):
                h1=maxl.pop(i)
    return maxl
    
# Conclusion and detailed comparison
result={}
for k in range(5,12):
    all=[]
    for type in types:
        if k in type:
            all=all+type[k]
    result[k]=maxgroup(all)
    print("dimension:",k)
    print("number of groups:",len(all))
    print("number of packings:",len(result[k]))
    images=[]
    graphstr=""
    for g in result[k]:
        h=copy(g)
        for u,v,l in h.edges():
            if l==3:
                h.set_edge_label(u,v,'')
        p=h.plot(vertex_labels=false,vertex_colors=coloring(g),vertex_size=30,edge_labels=true,edge_color='red',edge_labels_background='transparent')
        images.append(p)
        graphstr=graphstr+h.graphviz_string(edge_labels=True)+"\n"
    n=floor(30/k)
    m=ceil(len(result[k])/n)     
    Image=graphics_array(images,m,n)
    filename="max"+'-'+str(g.order())+".eps"
    Image.save(filename,figsize=[8,k*m/4])
    #Image.show()
    print(filename,"output")
    filename="max"+'-'+str(g.order())+".txt" # Output txt file
    with open(filename, "w") as f:
        f.write(graphstr)
        f.close()
