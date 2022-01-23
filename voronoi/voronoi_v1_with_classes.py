import numpy as np
import matplotlib.pyplot as plt
from random import random

# Author: Anthony Maltsev
# Date: 01/06/2022
#
# Description: Generate and display a randomly generated voronoi diagram 
# by finding a Delaunay triangulation of a random set of points and then 
# finding the corresponding/dual Voronoi diagram.

################################################################### define helper functions and classes below

class Edge() :
    def __init__(self, v1, v2) :
        self.v1 = v1
        self.v2 = v2
        self.length = euclidean_dist(self.v1, self.v2)
        self.inc_polygons = []
        self.ref()
    
    def equals(self, other) :
        if (self.v1.equals(other.v1) and self.v2.equals(other.v2)) or (self.v1.equals(other.v2) and self.v2.equals(other.v1)) :
            return True
        else :
            return False

    def ref(self) :
        if not self in self.v1.inc_edges :
            self.v1.inc_edges.append(self)
        if not self in self.v2.inc_edges :
            self.v2.inc_edges.append(self)
    def deref(self) :
        self.v1.inc_edges.remove(self)
        self.v2.inc_edges.remove(self)

    def replace_vertex(self, old, new) :
        self.deref()
        if old.equals(self.v1) :
            self.v1 = new
        else :
            self.v2 = new
        self.length = euclidean_dist(self.v1, self.v2)
        self.ref()
    
    def __str__(self) :
        return "(("+str(self.v1.x)+", "+str(self.v1.y)+"), ("+str(self.v2.x)+", "+str(self.v2.y)+"))"

class Vertex() :
    def __init__(self, x, y) :
        self.x = x
        self.y = y
        self.inc_edges = []
    
    def equals(self, other) :
        if (self.x == other.x and self.y == other.y) : # abs(x1 - x2 ) < 0.001
            return True
        else :
            return False

class Polygon() :
    def __init__(self, edges, vertices) :
        self.edges = edges
        self.vertices = vertices
        self.ref()
    def ref(self) :
        for edge in self.edges :
            edge.inc_polygons.append(self)
    def deref(self) :
        for edge in self.edges :
            edge.inc_polygons.remove(self)
    def deref_edges(self) :
        for edge in self.edges :
            edge.deref()


class Triangle(Polygon) : # pass vertices in and this creates a triangle that stores all info inc. circumcircle
    def __init__(self, vertices, edges=None, first_edge=None) :
        if (edges == None) and (first_edge == None) :
            self.vert__init__(vertices)
        elif (edges == None) :
            self.vert_and_firstedge__init__(vertices, first_edge)
        else :
            self.vert_and_edges__init__(edges, vertices)
    def vert__init__(self, vertices) :
        e1 = Edge(vertices[0], vertices[1])
        e2 = Edge(vertices[0], vertices[2])
        e3 = Edge(vertices[1], vertices[2])
        self.edges = [e1, e2, e3]
        self.vertices = vertices
        temp = Polygon(self.edges, self.vertices)
        cc_info = circumcenter_and_radius(temp)
        self.circumcenter = cc_info[0]
        self.circumradius = cc_info[1]
        self.ref()
    def vert_and_edges__init__(self, edges, vertices) :
        self.edges = edges
        self.vertices = vertices
        temp = Polygon(self.edges, self.vertices)
        cc_info = circumcenter_and_radius(temp)
        self.circumcenter = cc_info[0]
        self.circumradius = cc_info[1]
        self.ref()
    def vert_and_firstedge__init__(self, vertices, first_edge) : # first edge MUST be the edge between vert[0] and vert[1]
        e1 = first_edge
        e2 = Edge(vertices[0], vertices[2])
        e3 = Edge(vertices[1], vertices[2])
        self.edges = [e1, e2, e3]
        self.vertices = vertices
        temp = Polygon(self.edges, self.vertices)
        cc_info = circumcenter_and_radius(temp)
        self.circumcenter = cc_info[0]
        self.circumradius = cc_info[1]
        self.ref()


# find the euclidean distance between two vertices
def euclidean_dist(v1, v2) :
    return ((v1.x - v2.x)**2 + (v1.y - v2.y)**2)**0.5

#check whether two points are far enough apart
def distance_okay(v1, v2):
    if euclidean_dist(v1, v2) >= min_dist : #euclidean distance >= max_dist
        return True
    else :
        return False

#return a tuple containing the circumcenter and circumradius of a given triangle
def circumcenter_and_radius(polygon) :
    assert len(polygon.edges) == 3 # must be triangle
    point_a = polygon.vertices[0]
    point_b = polygon.vertices[1]
    point_c = polygon.vertices[2]
    dist_ab = euclidean_dist(point_a, point_b) # dont use pre-calculated distances from edges to ensure edge is correct, e.g. edges[0] = ab
    dist_ac = euclidean_dist(point_a, point_c) 
    dist_bc = euclidean_dist(point_b, point_c) 
    angle_a = np.arccos((dist_ab**2 + dist_ac**2 - dist_bc**2)/(2*dist_ab*dist_ac)) # law of cosines
    angle_b = np.arccos((dist_ab**2 + dist_bc**2 - dist_ac**2)/(2*dist_ab*dist_bc))
    angle_c = np.arccos((dist_bc**2 + dist_ac**2 - dist_ab**2)/(2*dist_bc*dist_ac))
    circumcenter = Vertex(
        (point_a.x*np.sin(2*angle_a) + point_b.x*np.sin(2*angle_b) + point_c.x*np.sin(2*angle_c))/(np.sin(2*angle_a) + np.sin(2*angle_b) + np.sin(2*angle_c)),
        (point_a.y*np.sin(2*angle_a) + point_b.y*np.sin(2*angle_b) + point_c.y*np.sin(2*angle_c))/(np.sin(2*angle_a) + np.sin(2*angle_b) + np.sin(2*angle_c))
    ) # formula for circumcenter found on google
    radius = euclidean_dist(circumcenter, point_a)
    return (circumcenter, radius)


# if the vertex is out of the range [0,1] for either x or y then return true
def out_of_bounds(vert) :
    x = vert.x
    y = vert.y
    if (x < 0 or x > 1) :
        return True
    elif (y < 0 or y > 1) :
        return True
    else :
        return False

#################################################################################### end helper functions

#### DEBUG
#v1 = Vertex(0.7623, 0.43221)
#v2 = Vertex(0.2367, 0.8173)
#e1 = Edge(v1, v2)
#e2 = Edge(v2, v1)
#assert e1.equals(e2)
#### END DEBUG

# generate a random set of points
def generate_points():
    point_set = []
    i = 0
    while i < num_points :
        new_point = Vertex(random(), random())
        if all([distance_okay(point, new_point) for point in point_set]) :
            point_set.append(new_point)
            i += 1
    return point_set

## Now do the Bowyer Watson triangulation algorithm on this set of points
# yields a set of triangles whose circumcenters are vertices on the voronoi diagram
# the idea of the algorithm : add points inductively and adjust our triangulation so that the bowyer-watson condition
# (defined as all of the vertices of any triangle are not inside the circumcenter of any other triangle)
# is maintained
def bowyer_watson_triangulation(point_set) :

    super_triangle = Triangle([Vertex(-0.5, -0.5), Vertex(-0.5, 3), Vertex(3, -0.5)]) # super triangle, triangle part
    triangle_set = [super_triangle] # store as a tuple of 3 items: the tuple of vertices(tuples), the tuple containing the circumcenter position, the radius of the circumcircle

    for point in point_set :
        #check all existing triangle to see if they are violated by the addition of point
        bad_triangle_set = []
        for triangle in triangle_set :
            # check if the point is within the circumcircle (subtract a bit to allow for floating point error)
            if (euclidean_dist(point, triangle.circumcenter) <= triangle.circumradius - 0.0001) : 
                bad_triangle_set.append(triangle)
        polygon_edge_set = []
        for bad_triangle in bad_triangle_set :
            # construct a polygon hole by getting all the edges that were on the outside of the violated region
            # i.e. all edges that occur only once in the bad_triangle_set / not shared by multiple bad triangles
            
            for bad_edge in bad_triangle.edges :
                if all([not bad_edge.equals(edge) for edge in polygon_edge_set]) : # check edge 1 unique so far
                    polygon_edge_set.append(bad_edge)
                else : # if not unique, remove from set
                    [polygon_edge_set.remove(edge) for edge in polygon_edge_set if edge.equals(bad_edge)]
                    bad_edge.deref() # make sure vertices don't reference this edge anymore since it doesn't exist
            triangle_set.remove(bad_triangle)
            bad_triangle.deref()
        # now add triangles constructed by connecting the point to each edge of the polygon hole, like slices of a pie
        for edge in polygon_edge_set :
            new_triangle = Triangle([edge.v1, edge.v2, point], first_edge=edge)
            triangle_set.append(new_triangle)
    # now we are done adding all points to our triangulation
    return triangle_set

def prune_triangles(triangle_set) :
    # remove any triangle that is connected to the original super_triangle (as that is extra leftover from construction process)
    tri_set_copy = triangle_set.copy()
    for triangle in tri_set_copy :
        if any([out_of_bounds(vert) for vert in triangle.vertices]) :
            triangle_set.remove(triangle)
            triangle.deref_edges()
            triangle.deref()
            
    return triangle_set


# now, get the dual voronoi map by looking at the circumcenters
def voronoi_diagram_from_triangulation(triangle_set) :
    voronoi_edge_set = []
    for triangle1 in triangle_set :
        for triangle2 in triangle_set :
            if not (triangle1.circumcenter.equals(triangle2.circumcenter)) : # check that 2 unique triangles by comparing circumcenters and hoping diff cc means diff triangle
                # if the two triangle have any shared edge, connect their circumcenters and thats an edge on the voronoi diagram
                if any([any([t1_edge.equals(t2_edge) for t2_edge in triangle2.edges]) for t1_edge in triangle1.edges]) :
                    voronoi_edge_set.append(Edge(triangle1.circumcenter, triangle2.circumcenter))
    return voronoi_edge_set
    
#prune edges that are out of bounds
def prune_bounds(edge_set) :
    edge_copy = edge_set.copy()
    for edge in edge_copy :
        if (out_of_bounds(edge.v1) or out_of_bounds(edge.v2)) :
            edge_set.remove(edge)
            edge.deref()
    return edge_set

def prune_small_edges(edge_set, min_size) :
    any_small = True
    while any_small :
        any_small = False
        for edge in edge_set :
            if (edge.length < min_size) :
                any_small = True
                edge_set.remove(edge)
                edge.deref()
                # shrink edge to one point (midpoint of edge) and combine
                mid_point = Vertex((edge.v1.x + edge.v2.x)/2, (edge.v1.y + edge.v2.y)/2)
                # make any edge that connected to any of the previous vertices connect to this midpoint instead
                v1_edges_copy = edge.v1.inc_edges.copy()
                v2_edges_copy = edge.v2.inc_edges.copy()
                for v1_edge in v1_edges_copy :
                    if (not v1_edge.equals(edge)) :
                        v1_edge.replace_vertex(edge.v1, mid_point)
                for v2_edge in v2_edges_copy :
                    if (not v2_edge.equals(edge)) :
                        v2_edge.replace_vertex(edge.v2, mid_point)
                break
    return edge_set

def display_edge_set(edge_set) :
    for edge in edge_set :
        #print(edge)
        plt.plot((edge.v1.x, edge.v2.x), (edge.v1.y, edge.v2.y), 'b')
    plt.show()

def display_triangle_set(triangle_set) :
    #visualise the triangles
    for triangle in triangle_set :
        plt.plot((triangle.edges[0].v1.x, triangle.edges[0].v2.x), (triangle.edges[0].v1.y, triangle.edges[0].v2.y), 'g')
        plt.plot((triangle.edges[1].v1.x, triangle.edges[1].v2.x), (triangle.edges[1].v1.y, triangle.edges[1].v2.y), 'g')
        plt.plot((triangle.edges[2].v1.x, triangle.edges[2].v2.x), (triangle.edges[2].v1.y, triangle.edges[2].v2.y), 'g')
    plt.show()

# now run all of the functions in a row!

#### PARAMETERS, CHANGE THESE IF YOU WANT
num_points = 100
min_dist = 0.02
min_edge = 0.05
####

point_set = generate_points()
triangle_set = bowyer_watson_triangulation(point_set)
#display_triangle_set(triangle_set)
triangle_set = prune_triangles(triangle_set)
#display_triangle_set(triangle_set)
voronoi_edges = voronoi_diagram_from_triangulation(triangle_set)
voronoi_edges = prune_bounds(voronoi_edges)
display_edge_set(voronoi_edges)
voronoi_edges = prune_small_edges(voronoi_edges, min_edge)
display_edge_set(voronoi_edges)
