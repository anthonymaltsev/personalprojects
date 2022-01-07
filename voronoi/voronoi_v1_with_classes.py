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
        v1.inc_edges.append(self)
        v2.inc_edges.append(self)
    
    def equals(self, other) :
        if (self.v1.equals(other.v1) and self.v2.equals(other.v2)) or (self.v1.equals(other.v2) and self.v2.equals(other.v1)) :
            return True
        else :
            return False


class Vertex() :
    def __init__(self, x, y) :
        self.x = x
        self.y = y
        self.inc_edges = []
    
    def equals(self, other) :
        if (self.x == other.x and self.y == other.y) :
            return True
        else :
            return False

class Polygon() :
    def __init__(self, edges, vertices) :
        self.edges = edges
        self.vertices = vertices
        for edge in self.edges :
            edge.inc_polygons.append(self)

class Triangle(Polygon) : # pass vertices in and this creates a triangle that stores all info inc. circumcircle
    def __init__(self, edges, vertices, circumcenter, circumradius) :
        super(edges, vertices)
        self.circumcenter = circumcenter
        self.circumradius = circumradius
    def __init__(self, vertices) :
        e1 = Edge(vertices[0], vertices[1])
        e2 = Edge(vertices[0], vertices[2])
        e3 = Edge(vertices[1], vertices[2])
        self.edges = [e1, e2, e3]
        self.vertices = vertices
        temp = Polygon(self.edges, self.vertices)
        cc_info = circumcenter_and_radius(temp)
        self.circumcenter = cc_info[0]
        self.circumradius = cc_info[1]
    def __init__(self, edges, vertices) :
        self.edges = edges
        self.vertices = vertices
        temp = Polygon(self.edges, self.vertices)
        cc_info = circumcenter_and_radius(temp)
        self.circumcenter = cc_info[0]
        self.circumradius = cc_info[1]



# find the euclidean distance between two vertices
def euclidean_dist(v1, v2) :
    return ((v1.x - v2.x)**2 + (v1.y - v2.y)**2)**0.5

#check whether two points are far enough apart
def distance_okay(v1, v2):
    if euclidean_dist(v1, v2) >= max_dist : #euclidean distance >= max_dist
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


def add_to_triangle_set(triangle_set, triangle) :
    triangle_set.append(triangle)

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

# generate a random set of points
def generate_points():
    point_set = []
    i = 0
    while i < num_points :
        new_point = (random(), random())
        if all([distance_okay(point, new_point) for point in point_set]) :
            point_set.append(new_point)
            i += 1
    return point_set

## Now do the Bowyer Watson triangulation algorithm on this set of points
# yields a set of triangles whose circumcenters are vertices on the voronoi diagram
# the idea of the algorithm : add points inductively and adjust our triangulation so that the bowyer-watson condition
# (defined as meaning all of the vertices of any triangle are not inside the circumcenter of any other triangle)
# is maintained
def bowyer_watson_triangulation(point_set) :

    super_triangle = ((-0.5, -0.5), (-0.5, 3), (3, -0.5)) # super triangle, triangle part
    triangle_set = [] # store as a tuple of 3 items: the tuple of vertices(tuples), the tuple containing the circumcenter position, the radius of the circumcircle
    add_to_triangle_set(triangle_set, super_triangle)

    for point in point_set :
        #check all existing triangle to see if they are violated by the addition of point
        bad_triangle_set = []
        for triangle in triangle_set :
            # check if the point is within the circumcircle (subtract a bit to allow for floating point error)
            if (euclidean_dist(point, triangle[1]) <= triangle[2] - 0.0001) : 
                bad_triangle_set.append(triangle)
        polygon_edge_set = []
        for bad_triangle in bad_triangle_set :
            # construct a polygon hole by getting all the edges that were on the outside of the violated region
            # i.e. all edges that occur only once in the bad_triangle_set / not shared by multiple bad triangles
            edge_1 = (bad_triangle[0][0], bad_triangle[0][1])
            edge_2 = (bad_triangle[0][0], bad_triangle[0][2])
            edge_3 = (bad_triangle[0][1], bad_triangle[0][2])
            if all([not edges_equal(edge_1, edge) for edge in polygon_edge_set]) : # check edge 1 unique so far
                polygon_edge_set.append(edge_1)
            else : # if not unique, remove from set
                polygon_edge_set.remove(edge_1)
            if all([not edges_equal(edge_2, edge) for edge in polygon_edge_set]) :
                polygon_edge_set.append(edge_2)
            else :
                polygon_edge_set.remove(edge_2)
            if all([not edges_equal(edge_3, edge) for edge in polygon_edge_set]) :
                polygon_edge_set.append(edge_3)
            else :
                polygon_edge_set.remove(edge_3)
            triangle_set.remove(bad_triangle)
        # now add triangles constructed by connecting the point to each edge of the polygon hole, like slices of a pie
        for edge in polygon_edge_set :
            new_triangle = (point, edge[0], edge[1])
            add_to_triangle_set(triangle_set, new_triangle)
    # now we are done adding all points to our triangulation
    return triangle_set

def prune_triangles(triangle_set) :
    # remove any triangle that is connected to the original super_triangle (as that is extra leftover from construction process)
    tri_set_copy = triangle_set.copy()
    for triangle in triangle_set.copy() :
        if any([out_of_bounds(vert) for vert in triangle[0]]) :
            triangle_set.remove(triangle)
    return triangle_set


# now, get the dual voronoi map by looking at the circumcenters
def voronoi_diagram_from_triangulation(triangle_set) :
    voronoi_edge_set = []
    for triangle1 in triangle_set :
        for triangle2 in triangle_set :
            if not (triangle1 == triangle2) : # check that 2 unique triangles
                triangle1_e1 = (triangle1[0][0], triangle1[0][1])
                triangle1_e2 = (triangle1[0][0], triangle1[0][2])
                triangle1_e3 = (triangle1[0][1], triangle1[0][2])
                triangle1_edges = [triangle1_e1, triangle1_e2, triangle1_e3]
                triangle2_e1 = (triangle2[0][0], triangle2[0][1])
                triangle2_e2 = (triangle2[0][0], triangle2[0][2])
                triangle2_e3 = (triangle2[0][1], triangle2[0][2])
                triangle2_edges = [triangle2_e1, triangle2_e2, triangle2_e3]
                # if the two triangle have any shared edge, connect their circumcenters and thats an edge on the voronoi diagram
                if any([any([edges_equal(t1_edge, t2_edge) for t2_edge in triangle2_edges]) for t1_edge in triangle1_edges]) :
                    voronoi_edge_set.append((triangle1[1], triangle2[1]))
    return voronoi_edge_set
    
#prune edges that are out of bounds
def prune_bounds(edge_set) :
    edge_copy = edge_set.copy()
    for edge in edge_copy :
        if (out_of_bounds(edge[0]) or out_of_bounds(edge[1])) :
            edge_set.remove(edge)
    return edge_set

def prune_small_edges(edge_set, min_size) :
    return edge_set

def display_edge_set(edge_set) :
    for edge in edge_set :
        plt.plot((edge[0][0], edge[1][0]), (edge[0][1], edge[1][1]), 'b')
    plt.show()

def display_triangle_set(triangle_set) :
    #visualise the triangles
    for triangle in triangle_set :
        plt.plot((triangle[0][0][0], triangle[0][1][0]), (triangle[0][0][1], triangle[0][1][1]), 'g')
        plt.plot((triangle[0][0][0], triangle[0][2][0]), (triangle[0][0][1], triangle[0][2][1]), 'g')
        plt.plot((triangle[0][1][0], triangle[0][2][0]), (triangle[0][1][1], triangle[0][2][1]), 'g')
    plt.show()

# now run all of the functions in a row!

#### PARAMETERS, CHANGE THESE IF YOU WANT
num_points = 100
max_dist = 0.05
min_edge = 0.05
####

point_set = generate_points()
triangle_set = bowyer_watson_triangulation(point_set)
#display_triangle_set(triangle_set)
triangle_set = prune_triangles(triangle_set)
#display_triangle_set(triangle_set)
voronoi_edges = voronoi_diagram_from_triangulation(triangle_set)
voronoi_edges = prune_bounds(voronoi_edges)
voronoi_edges = prune_small_edges(voronoi_edges)
display_edge_set(voronoi_edges)
