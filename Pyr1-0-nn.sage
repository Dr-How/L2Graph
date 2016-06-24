load("Matrix.sage")

def newL2(g0,g1,u0,u1):
    g=g0.disjoint_union(g1)
    g.merge_vertices([(0,u0),(1,u1)])
    # G=copy(g)
    # G.relabel()
    # if M(G).rank()!=G.order()-1:
    #     print "Hey!"
    #     return None
    for v in g:
        if v!=(0,u0):
            G=copy(g)
            G.delete_vertex(v)
            G.relabel()
            if M(G).rank()!=G.order():
                print g0,g1,u0,u1
                return false
    level1=true
    for v0 in g0:
        for v1 in g1:
            if v0!=u0 and v1!=u1:
                h=copy(g)
                h.delete_vertices([(0,v0),(1,v1)])
                h.relabel()
                if not L1(h):
                    return None
                if minEig(M(h))<0:
                    level1=false
    if level1:
        # print "Level 1:"
        # print u0,g0.edges()
        # print u1,g1.edges()
        return None
    else:
        return g

Result=[]
Output=[]
Level1=load("L1Hinge.sobj")

for i in range(len(Level1)):
    print i
    for j in range(i,len(Level1)):
        g0,u0=Level1[i]
        g1,u1=Level1[j]
        if g0.get_vertex(u0) and g1.get_vertex(u1):
            g=newL2(g0,g1,u0,u1)
            if g:
                checked=false
                for h in Result:
                    if g.is_isomorphic(h,edge_labels=true):
                        checked=true
                if not checked:
                    print i,j
                    print u0,g0.edges()
                    print u1,g1.edges()
                    Result.append(g)
                    Output.append([i,j,u0,u1])
Relation=Graph()
Relation.allow_loops(True)
for data in Output:
    Relation.add_edge(data[0],data[1])
