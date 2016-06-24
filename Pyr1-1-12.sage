load("Matrix.sage")

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
        L=[]
        for l in label:
            if not (l.rhs()+1).is_zero() and not (l.rhs()+1).is_positive():
                L.append(l.rhs().simplify_rational())
    return L

k=var('k')
for i in range(2,7):
    for j in range(max(i,3),7):
        g=Graph([(0,1,i),(1,2,j),(2,0,k)])
        if i==2:
            g.delete_edge(0,1)
        print i,j
        L=3
        while (1/i+1/L>=1/2) and (1/j+1/L>=1/2):
            L+=1
        for l in range(2,L):
            for m in range(2,L):
                for n in range(m,L):
                    G=copy(g)
                    u=G.add_vertex()
                    if l>2:
                        G.add_edge(1,u,l)
                    v=G.add_vertex()
                    if m>2:
                        G.add_edge(1,v,m)
                    w=G.add_vertex()
                    if n>2:
                        G.add_edge(1,w,n)

                    print " ", l,m,n

                    if not(l==m==2):
                        Lv=range(2,ceil(1/(1-1/m-1/l)))
                    else:
                        Lv=range(2,6)
                    if not(l==n==2):
                        Lw=range(2,ceil(1/(1-1/n-1/l)))
                    else:
                        Lw=range(2,6)
                    Lu=LabelPrism(g,{1:m},{1:n})
                    for lv in Lv:
                        for lw in Lw:
                            for lu in Lu:
                                if not (lv==lw==l==2) and not (m==n==2):
                                    H=copy(G)
                                    if lv>2:
                                        H.add_edge(u,v,lv)
                                    if lw>2:
                                        H.add_edge(u,w,lw)
                                    H.add_edge(v,w,lu)
                                    H.relabel()
                                    try:
                                        r=M(H).determinant().find_root(-1,-0.9,xtol=1e-16)
                                        print "  ", lv,lw,lu
                                        print "  ", r, pi.n()/arccos(-r)
                                    except:
                                        r=None
