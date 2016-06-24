# This program contains functions that is used by other programs.

# from sage.graphs.graph_plot import _line_embedding
from sage.graphs.graph_plot import _circle_embedding

def M(g): # Coxeter matrix from graph
    d=g.order()
    m=matrix(SR,d)
    for i in range(d):
        m[i,i]=1
    for u,v,l in g.edge_iterator():
        if l>0:
            m[u,v]=-cos(pi/l)
        else:
            m[u,v]=l
        m[v,u]=m[u,v]
    return m

def minEig(Mat): # Minimum eigenvalue
    m=copy(Mat).change_ring(RDF)
    mineig=min(m.eigenvalues())
    if abs(mineig)<0.000001:
        #print log(abs(mineig),10)
        return 0
    return mineig

def L1(g): # Recognise graphs of level <=1
    d=g.order()
    M0=M(g)
    for i in range(d):
        ind=range(d)
        ind.remove(i)
        M1=M0[ind,ind]
        if minEig(M1)<0:
            return false
    return true

def L2(g): # Recognise graphs of level 2
    d=g.order()
    if L1(g):
        return false
    M0=M(g)
    for i in range(d):
        for j in range(i+1,d):
            ind=range(d)
            ind.remove(i)
            ind.remove(j)
            M2=M0[ind,ind]
            if minEig(M2)<0:
                return false
    return true

