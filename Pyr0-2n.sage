# Enumerate products of a simplex and a triangle, of level 2.
# The objective is to eliminate labels 7, 9 and >= 11
# Following (Esselmann 2004), Section 4.1.

load("Matrix.sage")

def DetTriangle(k,l,m):
    return 1-(cos(pi/k)^2+cos(pi/l)^2+2*cos(pi/k)*cos(pi/l)*cos(pi/m))/sin(pi/m)^2
# There is a 2 before the triple product!!
# Esselmann made a mistake / typo.

Ports=load("L1Port.sobj")

for g in Ports:
    # print g.edges()
    G=copy(g)
    G.relabel()
    A1=M(G).determinant()
    for v in g:
        h=copy(g)
        h.delete_vertex(v)
        h.relabel()
        A0=M(h).determinant()
        r=(A1/A0)
        for i in range(3,7):
            Target=cos(pi/i)^2
            for k in range(2,7):
                for l in range(max(k,3),7):
                    m=7
                    while (Target-r*DetTriangle(k,l,m)).n()>0.0000001:
                        m=m+1
                    if abs((r*DetTriangle(k,l,m)-Target).n())<0.0000001:
                        print g.edges()
                        print v,i,k,l,m
print("-----------")
for g in Ports:
    # print g.edges()
    G=copy(g)
    G.relabel()
    A1=M(G).determinant()
    for u in g:
        for v in range(u+1,g.order()):
            capu=g.get_vertex(u)
            capv=g.get_vertex(v)
            if (capu<0 and capv<0):
                h=copy(g)
                h.add_edge(u,g.order(),3)
                h.add_edge(v,g.order(),3)
                h.relabel()
                A2=M(h).determinant()
                r=(A2/A1).n()
                for k in range(2,7):
                    for l in range(max(k,3),7):
                        m=7
                        while (r+DetTriangle(k,l,m)-1).n()>0.0000001:
                            m=m+1
                        if abs((r+DetTriangle(k,l,m)-1).n())<0.0000001:
                            print g.edges()
                            print u,v,k,l,m
