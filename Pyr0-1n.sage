load("Matrix.sage")
L2=load("L2SimplexR.sobj")
L2+=load("L2Simplex.sobj")

Total=0

for g in L2:
    B=M(g).change_ring(RDF)
    W=B.I
    G,O=g.automorphism_group(edge_labels=true, orbits=true)
    prism=false
    total=0
    for o in O:
        v=o[0]
        norm=W.row(v)*B*W.column(v)
        if norm>0.00000001:
            total+=len(o)
            h=copy(g)
            d1=M(h).determinant()
            h.delete_vertex(v)
            h.relabel()
            d2=M(h).determinant()
            x=d1.n()/d2.n()
            if x>1.00000001:
                prism=true
                Total+=1
            # elif x<0.9999999:
            #     print g.edges(),v
    if total==1 and prism:
        print g.edges()
        print O
        Total-=1

print Total
