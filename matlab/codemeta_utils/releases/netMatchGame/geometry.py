import numpy as np
from math import pi
import matplotlib.pyplot as plt
import ipdb


def orthogonal_complement_2d(v):
    """
    Returns the orthogonal complement to the vector v, obtained by rotating the vector of pi/2 counterclockwise.
    Works only for 2-D.
    """
    return cnormalize(np.vstack((-v[1],v[0])))

def orthogonal_complement_2d_test():
    v=np.random.randn(2,1)
    print v.T.dot(orthogonal_complement_2d(v)), v.T.dot(v)

def cnormalize(a):
    """ 
    Returns array of same dimension with each column having norm one.
    If a column has zero norm (more precisely, less than eps), it is left as it is.
    """
    #compute norms
    n=[np.linalg.norm(c) for c in a.T]
    eps=np.finfo('float').eps

    #do not touch vectors with norm zero (they will be divided by one instead of the norm)
    n=[x if x>eps else 1.0 for x in n]

    #return the original array with columns normalized (uses broadcasting)
    return a/n

class Line2D:
    def __init__(self,p1,p2):
        #Stores lines as origin, pointing vector, and a basis for the orthogonal complement
        self.origin=p1
        self.vector=p2-p1
        self.orthogonal=orthogonal_complement_2d(self.vector)

    def signed_distance(self,p):
        return np.asscalar(self.orthogonal.T.dot(p-self.origin))


def line_signed_distance_test():
    p1=np.array([[1],[0]])
    p2=np.array([[0],[1]])
    p=np.array([[2],[1]])
    o=np.zeros((2,1))
    print Line2D(o,p1).signed_distance(p)
    print Line2D(p1,o).signed_distance(p)
    print Line2D(o,p2).signed_distance(p)
    print Line2D(p2,0).signed_distance(p)

class Segment:
    def __init__(self,p1,p2):
        self.p1=p1
        self.p2=p2
        self.vector=p2-p1

    def linspace(self,num=10,endpoints=True,delta_fraction=0):
        """ 
        Generate equispaced points on the segment. If endpoints=false, then the endpoints of the segment do not appear in the returned array.
        """
        if endpoints:
            t=np.linspace(delta_fraction,1-delta_fraction,num)
        else:
            t=np.linspace(delta_fraction,1-delta_fraction,num+2)
            t=np.delete(t,[0,num+1])
        return self.p1+self.vector*t

class Polygon:
    def __init__(self,p):
        """ Stores the polygon mainly as an array of vertices"""
        self.vertices=p
        self.nb_vertices=p.shape[1]

    def plot(self):
        vertices=self.vertices_extended
        plt.plot(vertices[0,:],vertices[1,:])

    @property
    def lines(self):
        """ Return the list of lines corresponding to the polygon """
        return [Line2D(self.vertex(i),self.vertex_next(i)) for i in range(0,self.nb_vertices)]
    
    def lines_distance(self,p):
        """ Return signed distances from all lines corresponding to the polygon as a 1-D NumPy array"""
        return np.array([l.signed_distance(p) for l in self.lines])

    @property
    def edges(self):
        """ Return the list of segments corresponding to the vertices """
        return [Segment(self.vertex(i),self.vertex_next(i)) for i in range(0,self.nb_vertices)]
        
    @property
    def vertices_extended(self):
        """ Returns the array of vertices padded by repeating the first vertex at the end """
        return np.hstack((self.vertices, self.vertices[:,[0]]))

    def vertex(self,i):
        """ Returns the coordinates of the i-th vertex"""
        return self.vertices[:,[i]]

    def vertex_next(self,i):
        """ Returns the coordinates of the vertex following the i-th (using modulo aritmetic) """
        i_next=(i+1)%self.nb_vertices
        return self.vertices[:,[i_next]]

    def distance_cvx_inside(self,p):
        """ 
        Compute the distance between a point and the polygon, but 
        WARNING: gives correct results only for convex polygons, and if p is inside the polygon
        """
        if not self.is_inside_cvx(p):
            r=float('inf')
        else:
            r=min(self.lines_distance(p))
        return r
        
    def is_inside_cvx(self,p):
        """
        Returns true if the point is inside the polygon, but
        WARNING: gives correct results only for solid (counterclockwise-oriented), convex polygons
        """
        return min(self.lines_distance(p))>0

    def bounding_box(self):
        """
        Returns an array with the lower left and upper right coordinates of the bounding box containing
        all the points of the polygon
        """
        return np.array([[min(self.vertices[0,:]), max(self.vertices[0,:])],[min(self.vertices[1,:]), max(self.vertices[1,:])]])

    def rand_cvx(self,nb_samples=1,max_tries=100):
        """
        Use rejection sampling to sample points from inside the polygon, but
        WARNING: relies on is_inside_cvx, so it works only for convex polygons.
        Each sample is subject to max_tries tries. If the maximum number of tries is exceeded for some samples,
        the total number of samples returned might be less than nb_samples.
        """

        #Samples will be accumulated in this array
        samples=np.zeros((2,0))

        #Compute parameters of transformations that adapt the range of np.random.rand to the bounding box
        box=self.bounding_box()
        M=np.diag([box[0,1]-box[0,0],box[1,1]-box[1,0]])
        m=box[:,[0]]
        
        #Generate samples
        for i_sample in range(0,nb_samples):
            for i_try in range(0,max_tries):
                s=M.dot(np.random.rand(2,1))+m
                if self.is_inside_cvx(s):
                    samples=np.hstack((samples,s))
                    break
                
        return samples

def polygon_regular(radius,nb_sides):
    t=np.linspace(0.0,2*pi,num=nb_sides+1,endpoint=True)
    vertices=radius*np.vstack((np.cos(t),np.sin(t)))
    vertices=np.delete(vertices,nb_sides,1)
    polygon=Polygon(vertices)
    return polygon

def polygon_regular_test():
    p=polygon_regular(1,12)
    plt.clf()
    p.plot()
    plt.show()
