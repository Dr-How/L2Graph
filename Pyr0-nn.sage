# Enumerate products of two simplices of level 2.

load("Matrix.sage")

def newL2(g0,g1,cabel):
    g=g0.disjoint_union(g1)
    g.add_edges(cabel)
    G=copy(g)
    G.relabel()
    if M(G).rank()!=G.order()-1:
        return false
    for v in G:
        H=copy(G)
        H.delete_vertex(v)
        H.relabel()
        if M(H).rank()!=G.order()-1:
            print g0,g1,cabel
            return false
    level1=true
    for v0 in g0:
        for v1 in g1:
            h=copy(g)
            h.delete_vertices([(0,v0),(1,v1)])
            h.relabel()
            if not L1(h):
                return false
            if minEig(M(h))<0:
                level1=false
    if level1:
        print "Level 1:"
        print g0.edges(),g1.edges(),cabel
    return not level1

def Check(Res,res):
    for G in res:
        g=G[0]
        checked=false
        for H in Res:
            h=H[0]
            if g.is_isomorphic(h,edge_labels=true):
                checked=true
        if not checked:
            Res.append(G)
    
def listL2(g0,g1,cabel):
    result=[]
    g=g0.disjoint_union(g1)

    pos={}
    pos0=g0.get_pos()
    pos1=g1.get_pos()
    for u in g0:
        pos[(0,u)]=(pos0[u][0]+1,pos0[u][1])
    for u in g1:
        pos[(1,u)]=(pos1[u][0]-1,pos1[u][1])
    # print pos
    g.set_pos(pos)

    g.add_edges(cabel)

    if newL2(g0,g1,cabel):
        print "Level 2:"
        print g0.edges(),g1.edges(),cabel
        result.append((g,cabel))
    for v0 in g0:
        for v1 in g1:
            cap0=g0.get_vertex(v0)
            cap1=g1.get_vertex(v1)
            if (cap0!=0 and cap1!=0 and (cabel==[] or ((0,v0),(1,v1))>max(cabel))):
                new_cabel=copy(cabel)
                h0=copy(g0)
                h1=copy(g1)
                n0=[c for c in cabel if c[0]==(0,v0)]
                n1=[c for c in cabel if c[1]==(1,v1)]
                if n0==[] and n1==[]:
                    for l in range(3,min(abs(cap0),abs(cap1))+1):
                        new_cabel.append(((0,v0),(1,v1),l))
                        if cap0>0:
                            h1.set_vertex(v1,0)
                        if cap1>0:
                            h0.set_vertex(v0,0)
                        Check(result,listL2(h0,h1,new_cabel))
                if len(n0)==1 and n1==[]:
                    neighbor=n0[0][1][1]
                    if g1.get_vertex(neighbor)<0 and cap1<0:
                        new_cabel.append(((0,v0),(1,v1),3))
                        h0.set_vertex(v0,0)
                        if cap0>0:
                            h1.set_vertex(v1,0)
                        Check(result,listL2(h0,h1,new_cabel))
                if len(n1)==1 and n0==[]:
                    neighbor=n1[0][0][1]
                    if g0.get_vertex(neighbor)<0 and cap0<0:
                        new_cabel.append(((0,v0),(1,v1),3))
                        h1.set_vertex(v1,0)
                        if cap1>0:
                            h0.set_vertex(v0,0)
                        Check(result,listL2(h0,h1,new_cabel))
                if len(n0)==1 and len(n1)==1:
                    neighbor0=n0[0][1][1]
                    neighbor1=n0[0][1][1]
                    if g0.get_vertex(neighbor0)<0 and cap0<0:
                        if g1.get_vertex(neighbor1)<0 and cap1<0:
                            new_cabel.append(((0,v0),(1,v1),3))
                            h1.set_vertex(v1,0)
                            h0.set_vertex(v0,0)
                            Check(result,listL2(h0,h1,new_cabel))
    return result

Result=[]
Output=[]
Ports=load("L1Port.sobj")
for i in range(len(Ports)):
    for j in range(i,len(Ports)):
        print i,j
        ListL2=listL2(Ports[i],Ports[j],[])
        if ListL2!=[]:
            Result.extend(ListL2)
            for G in ListL2:
                Output.append([i,j,G[1]])
                print G[0].edges(),G[0].order()
