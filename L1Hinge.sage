load("Matrix.sage")

def Hinge(g,u):
    B=M(g).change_ring(RDF)
    W=B.I
    for w in g:
        norm=W.row(w)*B*W.column(w)
        if w==u:
            if abs(norm)>0.000001:
                return false
        else:
            if norm>-0.000001:
                return false
    return true

L1Simplex=load("L1Simplex.sobj")
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
Image.save("L1Hinges.eps", axes=false, figsize=[7.5,9])
