load("Matrix.sage")

def OnlyTwoReal(g,u,v):
    B=M(g).change_ring(RDF)
    W=B.I
    for w in g:
        norm=W.row(w)*B*W.column(w)
        if w==v or w==u:
            if norm<0.0001:
                return false
        else:
            if norm>0.0001:
                return false
    return true

def Hinge(g,u,v):
    if OnlyTwoReal(g,u,v):
        M0=M(g)
        ind=range(g.order())
        ind.remove(u)
        ind.remove(v)
        M2=M0[ind,ind]
        if minEig(M2)==0:
            return true
    return false

L2Simplex=load("L2Simplex.sobj")
Hinges=[]
Images=[]

for g in L2Simplex:
    for u in range(g.order()):
        for v in range(u+1,g.order()):
            if Hinge(g,u,v):
                print u,v,g.edges(),g
                Hinges.append((g,u,v))
                V=g.vertices()
                V.remove(u)
                V.remove(v)
                colors={"black":V, "lightgray":[u], "white":[v]}
                h=copy(g)
                print u,h.edges()
                for x,y,l in g.edges():
                    if l==3:
                        h.set_edge_label(x,y,'')
                image=h.plot(vertex_labels=false, vertex_size=20, vertex_colors=colors, edge_labels=true)
                Images.append(image)
                u=g.order()
                v=g.order()

# n=0
# Images2=[]
# pos={0:(1,1),1:(1,-1),2:(-1,1),3:(-1,-1)}
# colors={"black":[2,3],"white":[0,1]}
# for ab in range(2,7):
#     for ac in range(2,7):
#         for bc in range(max(ac,3),7):
#             for ad in range(ac,7):
#                 for bd in range(bc,7):
#                     if (1/ab+1/ac+1/bc)>=1 and (1/ab+1/ad+1/bd)>=1:
#                         g=Graph([(0,1,ab),(0,2,ac),(0,3,ad),(1,2,bc),(1,3,bd),(2,3,-1)])
#                         g.set_pos(pos)
#                         for x,y,l in g.edges():
#                             if l==3:
#                                 g.set_edge_label(x,y,'')
#                             if l==2:
#                                 g.delete_edge(x,y)
#                         image=g.plot(vertex_labels=false, vertex_size=15, vertex_colors=colors, edge_labels=true)
#                         Images2.append(image)
#                         n+=1
# 
# print n

# Image2=graphics_array(Images2,9,6)
# Image2.show(axes=false, figsize=[6,9])

Image=graphics_array(Images,11,7)
# Image.show(axes=false, figsize=[10,8])
Image.save("L2Hinges.eps", axes=false, figsize=[7,11])
