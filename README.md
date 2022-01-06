# Personal Projects
A Collection of Some Non-School Projects

### Voronoi Diagrams
A python program for generating random Voronoi Diagrams, which given a collection of points on the plane will divide that plane into sections depending on which point it is closest to. The purpose of this project is to potentially be used as a grid type for Slitherlink puzzles (see https://ejelta.com/android/slitherlink/).  
To generate Voronoi Diagrams, I implemented the Bowyer-Watson algorithm to construct a triangulation of a set of randomly generated points. Interestingly, there is a duality between this triangulation and the Voronoi Diagram for those given points: connecting the circumcenters of adjacent triangles yields the edges on the Voronoi diagram.  
Read more at: https://en.wikipedia.org/wiki/Delaunay_triangulation  
