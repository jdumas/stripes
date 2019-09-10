StripePatterns - Keenan Crane (2015)
====================================

This code implements the algorithm described in

   > Knöppel, Crane, Pinkall, Schröder, "Stripe Patterns on Surfaces" (ACM Transactions on Graphics 2015)

It is not particularly well-written, well-documented, or well-optimized, but it does the right thing.
The basic function of this  program is to compute a pattern of evenly-spaced stripes aligned with a
given direction field.  In principal the direction field could be computed using any method, but for
convenience the code implements a couple techniques that make it easy to generate a nice direction
field on pretty much any surface.  In particular, it implements variants of the algorithms described in

   1. Knöppel, Crane, Pinkall, Schröder, "Globally Optimal Direction Fields" (ACM Transactions on Graphics 2013)
   2. Crane, Desbrun, and Schröder, "Trivial Connections on Discrete Surfaces" (SGP 2010)

The first algorithm is used to initialize the field with either (a) the smoothest possible line field on
the surface or (b) a field aligned with principal curvature directions.  The second method allows one to
manually change the placement of singularities in the initial field, if desired.  The stripe pattern is
updated interactively with changes in the field.


Building
--------

The code has two dependencies:

   1. the [SuiteSparse](http://faculty.cse.tamu.edu/davis/suitesparse.html) library for numerical linear algebra, and
   2. OpenGL/GLUT.

Both of these libraries are quite standard, but building and linking them can vary quite a bit on different platforms. Alternatively, if you want to use your own linear algebra library (like Eigen), you can simply change the implementation of methods like solvePositiveDefinite() in DDG::SparseMatrix and DDG::DenseMatrix, which currently serve as wrappers around SuiteSparse. OpenGL and GLUT should be available by default on many platforms.

Once the dependencies have been installed/built, the code can be built by simply typing

   ```
   mkdir build; cd build
   cmake ..
   build -j 4
   ```

in the root directory.  The result will be an executable called `StripePatterns` (or `StripePatterns` on Windows).


Running
-------

To run the code from the command line, type

   ```
   ./StripePatterns input.obj
   ```

where “input.obj” is a path to a Wavefront OBJ file. The mesh must be connected and manifold, possibly with boundary. You will then be presented with a window displaying your mesh; hitting the spacebar will automatically generate a stripe pattern. Further commands are listed below.


User Interface
--------------

```
SPACE - compute a stripe pattern aligned with the globally smoothest direction field
c     - compute a stripe pattern aligned with the minimum principal curvature direction
d     - toggle display of input direction field
s     - set draw mode to smooth shaded
f     - set draw mode to wireframe
w     - write mesh to the file "out.obj"; the stripe pattern is stored in the texture coordinates
1     - compute just a 1D stripe pattern
2     - compute two orthogonal coordinates to get a 2D parameterization (not visualized, but will be saved to disk)
e     - edit the input direction field
*     - toggle singularities
-     - decrease stripe frequency
+     - increase stripe frequency
(     - rotate field left
)     - rotate field right
TAB   - animate rotating stripe pattern
r     - reload the input mesh from disk
\     - write a screenshot to the frames/ directory
ESC   - exit

ALT-CLICK: in edit mode, alt-clicking on a pair of triangles will adjust the singularities; in particular:
   --clicking on a pair of nonsingular triangles will create an equal and opposite pair
   --clicking on an equal and opposite pair will remove both singularities
   --clicking on a singularity and a nonsingular triangle will move the singularity
```
