load("Matrix.sage")

def Ideal(g,u):
    B=M(g).change_ring(RDF)
    W=B.I
    return abs(W.row(u)*B*W.column(u))<0.000001

def OnlyReal(g,v):
    B=M(g).change_ring(RDF)
    W=B.I
    for u in g:
        norm=W.row(u)*B*W.column(u)
        if (u!=v and norm>0.0001) or (u==v and norm<0.0001):
            return false
    return true

def Capacity(g,v):
    h=copy(g)
    cap=2
    h.add_edge(v,g.order(),cap+1)
    cap0=2
    while L2(h): # Confirm that v is a port
        cap+=1
        if OnlyReal(h,g.order()):
            cap0=cap
        h.set_edge_label(v,g.order(),cap+1) # Calculate the capacity of v
    return cap0,cap

Level1=load("L1Simplex.sobj")
# Ideals=[]
# Images=[]

Total=0
total=0

for g in Level1:
    has_Ideal=false
    # colors={"black":[], "white":[]}
    for v in g:
        cap0,cap=Capacity(g,v)
        if Ideal(g,v):
            # g.set_vertex(v,cap)
            # colors["white"].append(v)
            has_Ideal=true
            # Total+=(cap-2)*(cap-1)/2+cap-2
        # else:
            # g.set_vertex(v,2)
            # colors["black"].append(v)
    if has_Ideal:
        # Ideals.append(g)
        G,O=g.automorphism_group(edge_labels=true, orbits=true)
        # print O
        for o in O:
            if Ideal(g,o[0]):
                cap0,cap=Capacity(g,o[0])
                # print cap0,cap
                # cap=g.get_vertex(o[0])
                Total+=(cap-2)*(cap-1)/2+cap-2
                Total-=(cap0-2)*(cap0-1)/2+cap0-2
                total+=cap-2

print Total
print total
