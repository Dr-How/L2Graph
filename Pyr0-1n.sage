# Enumerate simplicial prisms of level 2.

load("Matrix.sage")
L2=load("L2SimplexR.sobj")
L2+=load("L2Simplex.sobj") # Load all (2,0)-graphs

Total=0 # Number of level-1 prisms

for g in L2:
    B=M(g).change_ring(RDF)
    W=B.I
    G,O=g.automorphism_group(edge_labels=true, orbits=true) # Avoid repeated computation
    prism=false
    total=0 # Number of space-like vertices
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
            x=d1.n()/d2.n() # Label of the dashed edge
            if x>1.00000001: # Otherwise, it's a pyramid
                prism=true
                Total+=1
            # elif x<0.9999999:
            #     print g.edges(),v
    if total==1 and prism: # If not the only space-like vertex, then level 1 (not 2)
        print g.edges()
        print O
        Total-=1

print Total
