load("Matrix.sage")

def OnlyReal(g,v):
    B=M(g).change_ring(RDF)
    W=B.I
    for u in g:
        norm=W.row(u)*B*W.column(u)
        if (u!=v and norm>0.0001) or (u==v and norm<0.0001):
            return false
    return true

def Port(g,v):
    h=copy(g)
    h.add_edge(v,g.order(),3)
    while L2(h) and OnlyReal(h,g.order()): # Confirm that v is a port
        h.set_edge_label(v,g.order(),h.edge_label(v,g.order())+1) # Calculate the capacity of v
    return h.edge_label(v,g.order())-1

Components=[]

g=Graph([(0,1,4),(1,2,5)])
g.set_pos({0:(1.5,1), 1:(0.5,0), 2:(1.5,-1)})
Components.append(g)
g=Graph([(0,1,5),(1,2,5)])
g.set_pos({0:(0.5,1), 1:(1.5,0), 2:(0.5,-1)})
Components.append(g)
g=Graph([(0,1,3),(1,2,8)])
g.set_pos({0:(0.5,0), 1:(1.5,0), 2:(2.5,0)})
Components.append(g)
g=Graph([(0,1,3),(1,2,10)])
g.set_pos({0:(0.5,0), 1:(1.5,0), 2:(2.5,0)})
Components.append(g)
g=Graph([(0,1,3),(1,2,3),(0,2,4)])
g.set_pos({0:(1.5,1), 1:(0.5,0), 2:(1.5,-1)})
Components.append(g)
g=Graph([(0,1,3),(1,2,3),(0,2,5)])
g.set_pos({0:(1.5,1), 1:(0.5,0), 2:(1.5,-1)})
Components.append(g)
g=Graph([(0,1,4),(1,2,4),(0,2,3)])
g.set_pos({0:(0.5,1), 1:(1.5,0), 2:(0.5,-1)})
Components.append(g)

for g in Components:
    for v in g:
        g.set_vertex(v,Port(g,v))

Bases=[]

Bases.append([0,3,[((0,1),(1,0),3)]])
Bases.append([0,5,[((0,1),(1,1),3)]])
Bases.append([1,3,[((0,0),(1,0),3),((0,2),(1,0),3)]])
Bases.append([1,5,[((0,0),(1,1),3),((0,2),(1,1),3)]])
Bases.append([2,4,[((0,0),(1,1),4)]])
Bases.append([4,4,[((0,1),(1,1),4)]])
Bases.append([2,2,[((0,0),(1,0),4)]])
Bases.append([6,6,[((0,0),(1,0),3),((0,2),(1,2),3)]])

def ConnectTo(g):
    c=[{}]
    for v in g:
        cap=g.get_vertex(v)
        for l in range(3,cap+1):
            c.append({v:l})
    for u in g:
        for v in range(u+1,g.order()):
            capu=g.get_vertex(u)
            capv=g.get_vertex(v)
            if (capu>0 and capv>0):
                h=copy(g)
                h.add_edge(u,g.order(),3)
                h.add_edge(v,g.order(),3)
                if L2(h) and OnlyReal(h,g.order()):
                    c.append({u:3,v:3})
    return c

Result=[]
Images=[]

for bases in Bases:
    g0=Components[bases[0]]
    g1=Components[bases[1]]
    g=g0.disjoint_union(g1)
    g.add_edges(bases[2])

    pos={}
    colors={"black":g.vertices(), "white":[0]}
    # gpos=copy(g)
    for v in g0:
        pos[(0,v)]=(-g0.get_pos()[v][0],g0.get_pos()[v][1])
    for v in g1:
        pos[(1,v)]=g1.get_pos()[v]
    
    apex=g.add_vertex()
    pos[apex]=(0,.5)
    g.set_pos(pos)
    C0=ConnectTo(g0)
    C1=ConnectTo(g1)
    for c0 in C0:
        for c1 in C1:
            G=copy(g)
            for v in c0:
                G.add_edge((0,v),apex,c0[v])
            for v in c1:
                G.add_edge((1,v),apex,c1[v])
            Level2=True
            h=copy(G)
            h.relabel()
            r=M(h).change_ring(RDF).eigenvalues()
            r.sort()
            if abs(r[1])>0.0001:
                Level2=False
            for v0 in g0:
                for v1 in g1:
                    h=copy(G)
                    h.delete_vertices([(0,v0),(1,v1)])
                    h.relabel()
                    if not L1(h):
                        Level2=False
            if Level2 and (c0!={} or c1!={}):
                checked=false
                for H in Result:
                    if G.is_isomorphic(H,edge_labels=true):
                        checked=true
                if not checked:
                    print r[1],bases,c0,c1
                    Result.append(G)
                    H=copy(G)
                    for u,v,l in G.edges():
                        if l==3:
                            H.set_edge_label(u,v,'')
                    image=H.plot(vertex_labels=false, vertex_size=30, vertex_colors=colors, edge_labels=true)
                    Images.append(image)

Image=graphics_array(Images,1,3)
Image.save("Pyr1-1-nn.eps", axes=false, figsize=[6,2])
Image.show()
