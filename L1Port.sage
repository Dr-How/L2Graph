# A vertex u of a (1,0)-graph H is called a port of H if there is a (2,0)-graph
# in the form of H + v in which v is the only real vertex and u is a neighbor
# of v.
# This program finds ports for (1,0)-graphs.

load("Matrix.sage")

def OnlyReal(g,v): # Tell if v is the only real vertex of g.
    B=M(g).change_ring(RDF)
    W=B.I
    for u in g:
        norm=W.row(u)*B*W.column(u)
        if (u!=v and norm>0.0001) or (u==v and norm<0.0001):
            return false
    return true

def Port(g,v): # Tell if v is a port of g.  The idea follows (Esselmann 1996), Lemma 4.2.
    h=copy(g)
    h.add_edge(v,g.order(),3)
    while L2(h) and OnlyReal(h,g.order()): # Confirm that v is a port
        h.set_edge_label(v,g.order(),h.edge_label(v,g.order())+1) # Calculate the capacity of v
    cap=h.edge_label(v,g.order())-1
    if cap>2:
        return cap # return the capacity of v.  This information is important when connecting g to other graphs.
    else:
        return 0

StrictL1=load("L1SimplexR.sobj") # Load cocompact (1,0)-graphs of >= 4 vertices.

pos={}
for j in range(3):
    pos[j]=(j,j)

# Triangle groups of label <= 7, with one label 2.
for i in range(3,7):
    for j in range(i,7):
        if (1/i+1/j<1/2):
            g=Graph([(0,1,i),(1,2,j)])
            g.set_pos(pos)
            StrictL1.append(g)

# Triangle groups of label 3~7
for i in range(3,7):
    for j in range(i,7):
        for k in range(j,7):
            if (1/i+1/j+1/k<1):
                g=Graph([(0,1,i),(1,2,j),(2,0,k)])
                _circle_embedding(g,g.vertices(),radius=1)
                StrictL1.append(g)

# Triangle groups with higher labels, to be eliminated
g=Graph([(0,1,3),(1,2,8)])
g.set_pos(pos)
StrictL1.append(g)

g=Graph([(0,1,3),(1,2,10)])
g.set_pos(pos)
StrictL1.append(g)

g=Graph([(0,1,4),(1,2,8)])
g.set_pos(pos)
StrictL1.append(g)

g=Graph([(0,1,4),(1,2,10)])
g.set_pos(pos)
StrictL1.append(g)

# the list of cocompact level-1 groups completed

Ports=[]
Images=[]

for g in StrictL1:
    have_port=false
    colors={"black":[], "white":[]}
    for v in g:
        cap=Port(g,v)
        g.set_vertex(v,cap)
        if cap>0:
            colors["white"].append(v)
            have_port=true
        else:
            colors["black"].append(v)
    for u in g:
        for v in range(u+1,g.order()): # test if it is possible to connect two vertices to the same vertex of another graph.
            capu=g.get_vertex(u)
            capv=g.get_vertex(v)
            if (capu!=0 and capv!=0):
                h=copy(g)
                h.add_edge(u,g.order(),3)
                h.add_edge(v,g.order(),3)
                if L2(h) and OnlyReal(h,g.order()):
                    g.set_vertices({u:-abs(capu),v:-abs(capv)}) # If yes, record it with negative cap
    if have_port:
        Ports.append(g)
        h=copy(g)
        for u,v,l in g.edges():
            if l==3:
                h.set_edge_label(u,v,'')
        image=h.plot(vertex_labels=true, vertex_size=110, vertex_colors=colors, edge_labels=true)
        Images.append(image)

for g in Ports:
    print g.edges()
    print g.get_vertices()
    print g.get_pos()

Image=graphics_array(Images,8,5)
Image.save("L1Ports.eps", axes=false, figsize=[10,16])
