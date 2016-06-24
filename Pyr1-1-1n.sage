load("Matrix.sage")

StrictL1=load("L1SimplexR.sobj")

pos={}
pos[0]=(-1,0)
pos[1]=(0,.5)
pos[2]=(1,0)
# pos[3]=(-1,2)
# pos[4]=(0,1.5)
# pos[5]=(1,2)

# for i in range(3,7):
#     for j in range(i,7):
#         if (1/i+1/j<1/2):
#             g=Graph([(0,1,i),(0,2,j)])
#             g.set_pos(pos)
#             StrictL1.append(g)
# 
# for i in range(3,7):
#     for j in range(i,7):
#         for k in range(j,7):
#             if (1/i+1/j+1/k<1):
#                 g=Graph([(0,1,i),(1,2,j),(2,0,k)])
#                 g.set_pos(pos)
#                 StrictL1.append(g)

def OnlyReal(g,v):
    B=M(g).change_ring(RDF)
    W=B.I
    for u in g:
        norm=W.row(u)*B*W.column(u)
        if (u!=v and norm>0.0000001) or (u==v and norm<0.0000001):
            return false
    return true

def Port(g,v):
    h=copy(g)
    h.add_edge(v,g.order(),3)
    while L2(h) and OnlyReal(h,g.order()): # Confirm that v is a port
        h.set_edge_label(v,g.order(),h.edge_label(v,g.order())+1) # Calculate the capacity of v
    return h.edge_label(v,g.order())-1

for g in StrictL1:
    for v in g:
        g.set_vertex(v,Port(g,v))

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
            if (capu>2 and capv>2):
                h=copy(g)
                h.add_edge(u,g.order(),3)
                h.add_edge(v,g.order(),3)
                if L2(h) and OnlyReal(h,g.order()):
                    c.append({u:3,v:3})
    return c

def LabelPrism(G,cv,cw):
    if cv=={} and cw=={}:
        return 0
    g=copy(G)
    v=g.add_vertex()
    w=g.add_vertex()
    for port in cv:
        g.add_edge(port,v,cv[port])
    for port in cw:
        g.add_edge(port,w,cw[port])
    M2=M(g).change_ring(SR)
    x=var('x')
    V=copy(M2[v]); V[w]=x; M2[v]=V
    W=copy(M2[w]); W[v]=x; M2[w]=W
    if not M2.determinant().is_zero():
        label=solve(M2.determinant()==0,x)
        # print label
        for l in label:
            if (l.rhs().n()+1)<-0.000001:
                return l.rhs()
    return 0

def LabelL3(G,cv,cw):
    result=[]
    g=copy(G)
    v=g.add_vertex()
    w=g.add_vertex()
    for port in cv:
        g.add_edge(port,v,cv[port])
    for port in cw:
        g.add_edge(port,w,cw[port])
    for l in range(2,7):
        h=copy(g)
        h.add_edge(v,w,l)
        Level3=True
        for u in G:
            H=copy(h)
            H.delete_vertex(u)
            H.relabel()
            if not L1(H):
                Level3=False
        if Level3:
            result.append(l)
    return result

Result=[]
Images=[]
Graphs=[]

for g in StrictL1:
    print g.edges()
    C=ConnectTo(g)
    for i in range(len(C)):
        print C[i]
        for j in range(len(C)):
            for k in range(j,len(C)):
                G=copy(g)
                u=G.add_vertex()
                for port in C[i]:
                    G.add_edge(port,u,C[i][port])
                v=G.add_vertex()
                for port in C[j]:
                    G.add_edge(port,v,C[j][port])
                w=G.add_vertex()
                for port in C[k]:
                    G.add_edge(port,w,C[k][port])

                Lv=LabelL3(g,C[i],C[j])
                Lw=LabelL3(g,C[i],C[k])
                l=LabelPrism(g,C[j],C[k])
                for lv in Lv:
                    for lw in Lw:
                        if not l.is_zero() and not (lv==2 and lw==2 and C[i]=={}):
                            H=copy(G)
                            if lv>2:
                                H.add_edge(u,v,lv)
                            if lw>2:
                                H.add_edge(u,w,lw)
                            H.add_edge(v,w,l.n())
                            H.relabel()
                            r=M(H).change_ring(RDF).eigenvalues()
                            r.sort()
                            if r[0]<-0.00001 and abs(r[1])<1e-10 and r[2]>0.00001:
                                checked=False
                                for h in Graphs:
                                    if H.is_isomorphic(h,edge_labels=true):
                                        checked=True
                                if not checked:
                                    Result.append((g,C[i],C[j],C[k],lv,lw))
                                    Graphs.append(H)
                                    print C[j], C[k], lv, lw, l, r[1]
                                    if g.order()==5:
                                        Pos=copy(g.get_pos())
                                        Pos[5]=(-2,2)
                                        Pos[6]=(2,-2)
                                        Pos[7]=(0,2)
                                        colors={"black":[0,1,2,3,4], "white":[5,6,7]}
                                        H.set_pos(Pos)
                                    if g.order()==4:
                                        Pos={0:(1,1),1:(-1,1),2:(-1,-1),3:(1,-1),4:(-.5,0),5:(.5,.5),6:(.5,-.5)}
                                        colors={"black":[0,1,2,3], "white":[4,5,6]}
                                        H.set_pos(Pos)
                                    if g.order()==3:
                                        Pos={0:(-1,-1),1:(-.5,0),2:(-1,1),3:(.5,0),4:(1,-1),5:(1,1)}
                                        colors={"black":[0,1,2], "white":[3,4,5]}
                                        H.set_pos(Pos)
                                    for u,v,l in H.edges():
                                        if l==3:
                                            H.set_edge_label(u,v,'')
                                        if l<0:
                                            H.set_edge_label(u,v,l.n(digits=4))
                                    image=H.plot(vertex_labels=false, vertex_size=30, vertex_colors=colors, edge_labels=true)
                                    Images.append(image)

Image=graphics_array(Images,3,6)
Image.save("Pyr1-1-1n.eps", axes=false, figsize=[12,6])
# Image=graphics_array(Images,3,6)
# Image.save("Pyr1-1-12.eps", axes=false, figsize=[13,6])
Image.show()
            # if l>2:
            #     G.add_edge(v0,v1,l)
            #     # print G.edges()
