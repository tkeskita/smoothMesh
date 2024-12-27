# smoothMesh

OpenFOAM mesh smoothing tool to improve mesh quality. Moves internal
mesh points by using the Centroidal smoothing algorithm (a version of the
[Laplacian smoothing algorithm](https://en.wikipedia.org/wiki/Laplacian_smoothing),
which uses surrounding cell centers instead of the neighbour point
locations to calculate the new point position). Optional heuristic
quality constraint options exist to constrain the smoothing, to avoid
self-intersections. No changes to mesh topology are made.

[smoothMesh demo video](https://vimeo.com/1023687267) illustrates mesh
smoothing on an example mesh, cross-section from which is shown
below. Without quality constraints, centroidal smoothing would move the
point highlighted in blue to the location highlighted with green,
which is outside of the domain, and therefore smoothing tends to
create self-intersecting mesh in the concave part of the geometry in
this case. Self-intersections like this can be avoided by using
additional quality constraints, which restrict the movement of vertices.

[![smoothMesh demo video](images/base_mesh_with_problematic_vertex.png)](https://vimeo.com/1023687267)

## Current restrictions

- Works on 3D polyhedron meshes
- Requires a consistent (not self-intersecting or tangled) initial mesh
- Smoothes internal mesh points only (boundary points are frozen)
- Developed on OpenFOAM.com v2312

## Compilation instructions

```
. /usr/lib/openfoam/openfoam2312/etc/bashrc
cd smoothMesh/src
wclean; wmake
```

## Command line options

- `-centroidalIters` specifies the number of smoothing iterations (default 20).

- `-maxStepLength` is the maximum length (in metres) for moving a point in one iteration (default 0.01). Adjust this value for your case. Smoothing process seems to be stable when this value is smaller than about one tenth of cell side length.

- `-orthogonalBlendingFraction` is the fraction (0 <= value <= 1) by which the edges touching the boundary faces are forced towards orthogonal direction. Zero value causes no orthogonal direction for boundary edges (default 0).

- `-qualityControl true` enables quality control constraints for smoothing. The constraints limit the tendency of smoothing to squish cells and create self-intersecting cells near concave geometry features. Without constraints, self-intersections can be created if the mesh contains skewed faces or low determinant cells (`cellDeterminant` according to `checkMesh`). When this option is enabled (default is true), the following options affect the results:

  - `-minEdgeLength` defines edge length below which edge points are fully frozen at their current location, but only if edge length would decrease during smoothing (default 0.05). Adjust this value for your case.

  - `-totalMinFreeze` option makes `-minEdgeLength` an absolute requirement, freezing short edges, even if edge length would increase during smoothing (default false). This option is useful to keep boundary layers in the mesh unmodified, and smooth the large cells only.

  - `-minAngle` defines the minimum edge-edge angle for face corners as well as the minimum of the face-face angles of the surrounding edges (in degrees, 0 < value < 180, default 45) below which points are fully frozen in their current location, if either the edge-edge angle or the face-face angle would **decrease** in smoothing. Points are allowed to move if both minimum angle values **increase** with smoothing, regardless of this value.

## Description of the algorithm

Please see [the algorithm description document](algorithm_description.md).

## Usage examples

Adjust `-maxStepLength` and `-minEdgeLength` according to your mesh cell size.

- Parallel run example: `mpirun -np 3 smoothMesh -centroidalIters 20 -maxStepLength 0.01 -minEdgeLength 0.05 -parallel`
- Serial run example: `smoothMesh -centroidalIters 20 -maxStepLength 0.01 -minEdgeLength 0.05`

## Getting help and feedback

Please use Github issues section for asking help. If you like this
tool, please star this repository in Github!
