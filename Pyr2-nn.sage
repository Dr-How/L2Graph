load("Matrix.sage")

def newL2(g0,g1,u0,u1,v0,v1):
    g=g0.disjoint_union(g1)
    g.merge_vertices([(0,u0),(1,u1)])
    g.merge_vertices([(0,v0),(1,v1)])
    G=copy(g)
    G.relabel()
    Eig=M(G).change_ring(RDF).eigenvalues()
    Eig.sort()
    if Eig[0]>-0.00001 or abs(Eig[1])>1e-10 or Eig[2]<0.00001:
        # print g0.edges(),g1.edges()
        # print u0,u1,v0,v1
        # print Eig[0],Eig[1]
        return None
    # level1=true

    h=copy(g)
    h.delete_vertex((0,u0))
    h.relabel()
    Eig=M(h).change_ring(RDF).eigenvalues()
    Eig.sort()
    if Eig[0]>-0.00001 or abs(Eig[1])>1e-10 or Eig[2]<0.00001:
        # print "Ouch!"
        return None

    h=copy(g)
    h.delete_vertex((0,v0))
    h.relabel()
    Eig=M(h).change_ring(RDF).eigenvalues()
    Eig.sort()
    if Eig[0]>-0.00001 or abs(Eig[1])>1e-10 or Eig[2]<0.00001:
        # print "Ouch!"
        return None

    for w0 in g0:
        for w1 in g1:
            if w0!=u0 and w0!=v0 and w1!=u1 and w1!=v1:
                h=copy(g)
                h.delete_vertices([(0,w0),(1,w1)])
                h.relabel()
                if not L1(h):
                    return None
                # if minEig(M(h))<0:
                #     level1=false

    # if level1:
    #     print "Level 1:"
    #     print g0.edges()
    #     print g1.edges()
    #     return None
    # else:
    return g

Result=[]
Output=[]
Level2=load("L2Hinge.sobj")

for i in range(len(Level2)):
    print i
    for j in range(i,len(Level2)):
        g0,u0,v0=Level2[i]
        g1,u1,v1=Level2[j]
        if (not g0.has_edge(u0,v0) and not g1.has_edge(u1,v1)) or (g0.has_edge(u0,v0) and g1.has_edge(u1,v1) and g0.edge_label(u0,v0)==g1.edge_label(u1,v1)):
            g=newL2(g0,g1,u0,u1,v0,v1)
            if g:
                checked=false
                for h in Result:
                    if g.is_isomorphic(h,edge_labels=true):
                        checked=true
                if not checked:
                    print i,j
                    # h=copy(g)
                    # h.relabel()
                    # print h.order()-M(h).rank()
                    print g0.edges()
                    print g1.edges()
                    Result.append(g)
                    Output.append([i,j,u0,u1,v0,v1])
            g=newL2(g0,g1,u0,v1,v0,u1)
            if g:
                checked=false
                for h in Result:
                    if g.is_isomorphic(h,edge_labels=true):
                        checked=true
                if not checked:
                    print i,j
                    # h=copy(g)
                    # h.relabel()
                    # print h.order()-M(h).rank()
                    print g0.edges()
                    print g1.edges()
                    Result.append(g)
                    Output.append([i,j,u0,v1,v0,u1])

Relation=Graph()
Relation.allow_loops(True)
for data in Output:
    Relation.add_edge(data[0],data[1])
    if data[0]!=data[1]:
        Relation.set_vertex(data[0],[data[2],data[4]])
        Relation.set_vertex(data[1],[data[3],data[5]])
# for g in Relation:
    # print "Hinge used by graph",g
    # for data in Output:
    #     if data[0]==g:
    #         print data[2],data[4]
    #     if data[1]==g:
    #         print data[3],data[5]
