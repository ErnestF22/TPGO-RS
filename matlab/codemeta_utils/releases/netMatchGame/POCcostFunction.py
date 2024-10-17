#!/bin/env python

import ipdb
import copy
import numpy as np
import matplotlib.pyplot as plt
import geometry
import graph2
import random
from scipy.optimize import minimize
from scipy.optimize import linprog


plt.ion()

def euclidean_distance_squared_matrix(a,b=None):
    if not b:
        b=a
        
    #get number of points
    nb_a=a.shape[1]
    nb_b=b.shape[1]

    #note: the vectors below are column for a and row for b
    
    #compute vectors of squared column norms
    a_sq=np.reshape(sum(a*a),(-1,1))
    b_sq=np.reshape(sum(b*b),(1,-1))

    #all-1 vectors
    u_a=np.ones((nb_a,1))
    u_b=np.ones((1,nb_b))

    #matrix of inner products
    a_b=a.T.dot(b)

    #distance matrix
    return a_sq.dot(u_b)-2*a_b+u_a.dot(b_sq)

    
def softmax(a,approximation_factor=20.0):
    """ 
    Compute the soft max of an array.
    The approximation_factor optional argument (default is 20.0) controls the "softness" (the higher the factor, the closer the approximation)
    """
    a=np.array(a)
    return np.log(sum(np.exp(a*approximation_factor)))/approximation_factor

def softmin(a,approximation_factor=20.0):
    """ 
    Compute the soft min of an array.
    The approximation_factor optional argument (default is 20.0) controls the "softness" (the higher the factor, the closer the approximation)
    """
    a=np.array(a)
    return -softmax(-a,approximation_factor=approximation_factor)

def softmax_test():
    t=np.arange(-2.0,2.0,0.1)
    s=[softmax([1,a]) for a in t]
    plt.clf()
    plt.plot(t,s)
    plt.show()

def softmin_test():
    t=np.arange(-2.0,2.0,0.1)
    s=[softmin([1,a]) for a in t]
    plt.clf()
    plt.plot(t,s)
    plt.show()

def cost_nearest_single(samples,i_sample,polygon):
    """ Compute the nearest-structure cost for a the i-th sample """
    #sample in question
    s=samples[:,[i_sample]]
    #concatenate distances to polygon's lines and other points
    d_polygon=polygon.lines_distance(s)
    d_samples=np.linalg.norm(np.delete(samples,i_sample,axis=1)-s,axis=0)
    d=np.concatenate((d_polygon,d_samples))
    
    return -softmin(d)

def cost_nearest(samples,polygon):
    cost=0
    samples=np.reshape(samples,(2,-1))
    for i_sample in range(0,samples.shape[1]):
       cost+=cost_nearest_single(samples,i_sample,polygon)
    return cost

def node_start(polygon):
    """ Generate starting locations equispaced on the edges of the polygon """
    v=np.zeros((2,0))
    for e in polygon.edges:
        v=np.hstack((v,e.linspace(4,endpoints=True,delta_fraction=0.125)))
    return v
        

def node_packing(nb_nodes=41):
    """ Optimizes position of the nodes to approximately maximize the minimum distance """
    p=geometry.polygon_regular(1,12)
    s0=p.rand_cvx(nb_nodes)
    plt.clf()
    p.plot()
    plt.plot(s0[0,:],s0[1,:],'*')
    plt.show()
    print cost_nearest(s0,p)
    res=minimize(cost_nearest,np.reshape(s0,(-1,)),args=(p),method='BFGS',options={'disp':True})
    s=np.reshape(res.x,(2,-1))
    plt.plot(s[0,:],s[1,:],'o')
    print cost_nearest(s,p)
    np.save('nodes',s)
    return s

def board_graph(nb_starts,nb_nodes):
    """ Create the graph for the board """
    g=graph2.Graph()

    starts=['s{}'.format(i) for i in range(0,nb_starts)]
    nodes=['n{}'.format(i) for i in range(0,nb_nodes)]
    g.add_vertex_list(starts)
    g.add_vertex_list(nodes)

    #add edges between starts and nodes
    g.add_edge_from_vertex_list(starts,nodes)
    #add edges between each node and the nodes with higher index
    for i_node in range(0,len(nodes)-1):
        g.add_edge_from_vertex_list([nodes[i_node]],nodes[i_node+1:])
    g.vnames=starts+nodes
    return g

def board_graph_weights(g,d):
    return np.array([[d[e[0],e[1]]] for e in g.edges_indeces()])

def board_solve_edges(g,v_starts,v_nodes,nb_edges_starts=2,nb_edges_nodes=3):
    #compute distances
    nb_nodes=v_nodes.shape[1]
    nb_starts=v_starts.shape[1]
    v=np.hstack((v_starts,v_nodes))
    d=euclidean_distance_squared_matrix(v)
    #get parameters for the problem
    A_ub=-g.incidence_matrix()
    b_ub=-np.vstack((nb_edges_starts*np.ones((nb_starts,1)),nb_edges_nodes*np.ones((nb_nodes,1))))
    w=board_graph_weights(g,d)
    #add box constraints
    nb_edges=w.shape[0]
    A_ub=np.vstack((A_ub,np.eye(nb_edges),-np.eye(nb_edges)))
    b_ub=np.vstack((b_ub,np.ones((nb_edges,1)),np.zeros((nb_edges,1))))

    #xtrue=np.ones((nb_edges,1))
    
    b_ub=np.reshape(b_ub,(-1,))
    w=np.reshape(w,(-1,))

    #solve the problem
    res=linprog(w, A_ub=A_ub, b_ub=b_ub,method='interior-point')
    xres=np.reshape(res.x,(-1,1))
    b_ub=np.reshape(b_ub,(-1,1))
    print np.hstack((A_ub.dot(xres),b_ub))

    res.A_ub=A_ub
    res.b_ub=b_ub
    res.xres=xres
    return res
    
def board_plot(g,v_starts,v_nodes):
    plt.plot(v_nodes[0,:],v_nodes[1,:],'o')
    plt.plot(v_starts[0,:],v_starts[1,:],'o')
    v=np.hstack((v_starts,v_nodes));
    for e in g.edges_indeces():
        e=list(e)
        plt.plot(v[0,e],v[1,e],'b-')
    plt.show()

def find_edges():
    p=geometry.polygon_regular(1,12)
    v_starts=node_start(p)
    v_nodes=np.load('nodes.npy')

    #v_starts=np.array([[0,1,0,1],[0,0,1,1]])
    #v_nodes=np.array([[0.5,0.5,0.25,0.75],[0.25,0.75,0.5,0.5]])

    nb_starts=v_starts.shape[1]
    nb_nodes=v_nodes.shape[1]
    g=board_graph(nb_starts,nb_nodes)

    res=board_solve_edges(g,v_starts,v_nodes,nb_edges_nodes=6)
    print res

    idx=list(np.where(res.x>0.5)[0])
    v=np.hstack((v_starts,v_nodes));
    e_masked=[g.edges_indeces()[i] for i in idx]

    plt.clf()
    p.plot()
    ax=plt.gca()
    vn=g.vertices()
    for e in e_masked:
        e=list(e)
        v0=v[:,e[0]]
        v1=v[:,e[1]]
        plt.plot([v0[0], v1[0]],[v0[1], v1[1]],'b-')
        #ax.text(v0[0],v0[1],'{}'.format(e[0]))
        #ax.text(v1[0],v1[1],'{}'.format(e[1]))

    np.save('starts',v_starts)
    np.save('edges',e_masked)
    
    #plt.clf()
    #board_plot(g_masked,v_starts,v_nodes)
    plt.show()
    
    #x=board_solve_edges(g,v_starts,v_nodes)
    #print x

def tikz_output():
    import os.path as path
    scale=8
    nb_sides=12
    nb_players=4
    path_tikz='latex/tikz'
    v_starts=np.load('starts.npy')
    v_nodes=np.load('nodes.npy')
    names_starts=name_coordinates('s',v_starts)
    names_nodes=name_coordinates('n',v_nodes)
    edges=np.load('edges.npy')
    v_polygon=geometry.polygon_regular(1,nb_sides).vertices
    names_polygon=name_coordinates('p',v_polygon)

    #scaling
    v_starts*=scale
    v_nodes*=scale
    v_polygon*=scale

    #constants
    eol=';\n'
    eol_table='\\\\\n'
    cycle=' -- cycle'+eol

    #board
    with open(path.join(path_tikz,'board.tex'),'w') as f:
        f.write('%This file has been generated by {}\n\n'.format(__file__))
        f.write('%Player sector colors\n')
        for player in range(0,nb_players):
            for i_sector in range(player,nb_sides,nb_players):
                v_sector=sector_coordinates(v_polygon,i_sector,1.1)
                f.write('\path[fill=player{}]'.format(player+1)+join_coordinates(v_sector)+cycle)
        f.write('%Polygon\n')
        f.write('\draw[polygon] '+join_coordinates(v_polygon)+cycle)
        f.write('%Nodes and links\n')
        f.write('\path[start] '+join_coordinates(v_starts,names_starts,' ',' coordinate[start] ')+eol)
        f.write('\path '+join_coordinates(v_nodes,names_nodes,' ',' coordinate[netnode] ')+eol)
        f.write('\draw[link] '+join_edges_names(names_starts+names_nodes,edges)+eol)
        
        f.write('%Named coordinates for starting locations\n')
        cnt=0
        for i_sector in range(0,3):
            for i_player in range(0,4):
                s='\path '
                for i_place in range(0,4):
                    s+='coordinate (Spl{}se{}pa{}) at ({},{}) '.format(i_player+1,i_sector,i_place,v_starts[0,cnt],v_starts[1,cnt])
                    cnt+=1
                s+=eol
                f.write(s)
        
    #card boxes    
    with open(path.join(path_tikz,'card-places.tex'),'w') as f:
        f.write('%This file has been generated by {}\n\n'.format(__file__))
        f.write('%Places for cards\n')
        f.write('\matrix[column sep=5mm]{\n')
        s_card_places=' & '.join(['\\node[card place,draw=player{0}] {{Player {0}}};'.format(i_player+1) for i_player in range(0,nb_players)])
        f.write(s_card_places+eol_table)
        f.write('}'+eol)

    #word cards
    language='english'
    nb_words_cols=5
    nb_words_rows=5
    nb_words_max=75
    with open(path.join(path_tikz,'word-card-table-{}.tex').format(language),'w') as f:
        with open('wordlist-{}-{}l.txt'.format(language,nb_players),'r') as fw:
             words=fw.read().splitlines()
             if len(words)>nb_words_max:
                 words=random.sample(words,nb_words_max)
                 
             nb_words=len(words)
             idx_rows=chunks(range(0,nb_words),nb_words_cols)
             
             #join elements into rows
             s_rows=[' & '.join([ words[idx].upper() for idx in row])+eol_table for row in idx_rows]
             #split rows into tables
             s_tables=['\\begin{tikzpicture}\n\\matrix[table]{\n'+''.join(sr)+'}'+eol+'\end{tikzpicture}\n' for sr in chunks(s_rows,nb_words_rows)]
             #write all tables
             for st in s_tables:
                 f.write(st)

def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in xrange(0, len(l), n):
        yield l[i:i + n]
        
def sector_coordinates(v_polygon,i_sector,scale1,scale2=1.0):
    """ Returns four coordinates that define an anular sector for a closed polygon """
    nb_vertices=v_polygon.shape[1]
    v_segment=v_polygon[:,[i_sector,(i_sector+1)%nb_vertices]]
    return np.hstack((scale1*v_segment,scale2*np.fliplr(v_segment)))
        
def name_coordinates(prefix,v):
    """ Return a list of strings of the form prefix+number corresponding to each column of v """
    names=[prefix+str(i) for i in range(0,v.shape[1])]
    return names
    
def join_coordinates(vertices,names=None,s=' -- ',sc=' coordinate '):
    """ 
    Assemble a string of the form (x,y) -- ... -- (x,y) with the coordinates in vertices. 
    If names is provided, then add also a named coordinate for each point.
    """
    
    if not names:
        r=s.join(['({},{})'.format(v[0],v[1]) for v in vertices.T])
    else:
        r=s.join([sc.join(['({},{})'.format(v[0],v[1]),'('+n+')']) for v,n in zip(vertices.T,names)])
    return r

def join_edges_names(name_vertices,edges,s=' edge '):
    #transform edge indeces into pairs of names
    edges_vertices=[(name_vertices[e[0]],name_vertices[e[1]]) for e in edges]
    #build string with sequence of edges
    return ' '.join(['({}) {} ({})'.format(vn1,s,vn2) for vn1,vn2 in edges_vertices])


def main():
    pass
    
    
if __name__=="__main__":
    main()
