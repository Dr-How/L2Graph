load("Matrix.sage")

L2Hinges=load("L2Hinge.sobj")
Total=0
Result=[]
Graphs=[]

for m in range(len(L2Hinges)):
    g=L2Hinges[m][0]
    u=L2Hinges[m][1]
    v=L2Hinges[m][2]
    print g.edges(),u,v
    if g.has_edge(u,v):
        l=g.edge_label(u,v)
    else:
        l=2
    L3=[]
    L3Graph=[]
    for lu in range(2,7):
        for lv in range(2,7):
            h=copy(g)
            w=h.add_vertex()
            h.add_edge(u,w,lu)
            h.add_edge(v,w,lv)
            Level3=True
            if 1/l+1/lu+1/lv<1:
                Level3=False
            for x in g:
                if x!=u and x!=v:
                    H=copy(h)
                    H.delete_vertex(x)
                    H.relabel()
                    if not L1(H):
                        Level3=False
            if Level3:
                checked=False
                for G in L3Graph:
                    if H.is_isomorphic(G,edge_labels=true):
                        checked=True
                if not checked:
                    G=copy(h)
                    G.relabel()
                    L3.append((lu,lv))
                    # Eig=M(G).change_ring(RDF).eigenvalues()
                    # Eig.sort()
                    # if Eig[0]>0 or Eig[1]<0:
                    #     print "Ouch!!!!!!!!!!!!"
                    # print lu,lv
                    L3Graph.append(H)
    print len(L3Graph)
    total=0
    for i in range(len(L3)):
        for j in range(i,len(L3)):
            if not (L3[i][0]==L3[j][0]==2 or L3[i][1]==L3[j][1]==2):
                h=copy(g)
                vi=h.add_vertex()
                vj=h.add_vertex()
                h.add_edge(vi,u,L3[i][0])
                h.add_edge(vi,v,L3[i][1])
                h.add_edge(vj,u,L3[j][0])
                h.add_edge(vj,v,L3[j][1])
                h.add_edge(vi,vj,-1)
                Eig=M(h).change_ring(RDF).eigenvalues()
                Eig.sort()
                # print Eig
                if Eig[0]<0.0001 and abs(Eig[1])<0.00000001:
                    chekced=False
                    for H in Graphs:
                        if H.is_isomorphic(h,edge_labels=True):
                            checked=True
                    if not checked:
                        total+=1
                        Result.append([m,L3[i],L3[j]])
                        print m,L3[i],L3[j]
                        Graphs.append(h)
    print total
    Total+=total
print Total
