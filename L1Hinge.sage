# If a (1,0)-graph has a unique ideal vertex, we call the vertex a hinge.
# This program finds (1,0)-graphs with a hinge.

load("Matrix.sage")

def Hinge(g,u):  # Tell if the vertex u is a hinge of g.
    B=M(g).change_ring(RDF)
    W=B.I
    for w in g:
        norm=W.row(w)*B*W.column(w)
        if w==u:
            if abs(norm)>0.000001: # If w is not light-like
                return false
        else:
            if norm>-0.000001: # If another vertex is light- or space-like (in the latter case, the graph is not of level 1).
                return false
    return true

L1Simplex=load("L1Simplex.sobj") # read list of (1,0)-graphs
Hinges=[]
Images=[]

for g in L1Simplex:
    for u in range(g.order()):
        if Hinge(g,u):
            Hinges.append((g,u))
            V=g.vertices()
            V.remove(u)
            colors={"black":V, "white":[u]}
            h=copy(g)
            print u,h.edges()
            for v,w,l in g.edges():
                if l==3:
                    h.set_edge_label(v,w,'')
            image=h.plot(vertex_labels=false, vertex_size=30, vertex_colors=colors, edge_labels=true)
            Images.append(image)
            u=g.order()

Image=graphics_array(Images,5,6)
# Image.show(axes=false, figsize=[7.5,9])
Image.save("L1Hinges.eps", axes=false, figsize=[7.5,9]) # Output graph
