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

## Current features and restrictions

- Works on 3D polyhedron meshes
- Requires a consistent (not self-intersecting or tangled) initial mesh with "good enough" quality
- Smoothes internal mesh points only (boundary points are frozen)
- Developed on OpenFOAM.com v2312

**Warning: This tool is still under active development!**

## Compilation instructions

```
. /usr/lib/openfoam/openfoam2312/etc/bashrc
cd smoothMesh/src
wclean; wmake
```

## Command line options

- `-centroidalIters` specifies the number of smoothing iterations (default 20).

- `-maxStepLength` is the maximum length (in metres) for moving a point in one iteration (default 0.01). Adjust this value for your case. Smoothing process seems to be stable when this value is smaller than about one tenth of minimum cell side length.

- `-minEdgeLength` defines edge length below which edge points are fully frozen at their current location, but only if edge length would decrease during smoothing (default 0.05). Adjust this value for your case.

- `-totalMinFreeze` option causes mesh points on all edges shorter than `-minEdgeLength` to freeze (default false). This option is useful to keep boundary layers in the mesh unmodified, and smooth the large cells only.

- `-orthogonalBlendingFraction` is the fraction (0 <= value <= 1) by which the edges touching the boundary faces are forced towards orthogonal direction. Zero value causes no orthogonal direction for boundary edges (default 0). Warning: This is an experimental option (WIP) which currently works correctly only for serial run!

The following options are related to additional **heuristic quality control constraints for smoothing**. The constraints work by disallowing movement of point (freezing of points) if the movement would cause quality of the mesh would suffer too much. Without constraints, centroidal smoothing may squish cells and create self-intersecting cells e.g. near concave geometry features, depending on the mesh details. Have a look at [the algorithm description document](algorithm_description.md) for details.

**Note:** The old `-qualityControl` option has been superceded by the options below.

- `-faceAngleConstraint` boolean option enables an additional quality control which restricts decrease of smallest and largest face-face angle (default is true). When this option is enabled, the `-minAngle` option defines the minimum angle (in degrees, default value 35), and the `-maxAngle` option specifies the maximum angle (in degrees, default value 170).

- `-minAngle` value causes point freezing *only* if the angle is below this value and if the angle would *decrease* in smoothing. Points are allowed to move if the angle value *increases* with smoothing, regardless of this value. The same applies for the `-maxAngle` option: Freezing takes place only if angle is above the specified value and if the angle would *increase* in smoothing.

## Description of the algorithm

Please view [the algorithm description document](algorithm_description.md).

## Usage examples

Adjust at least the `-maxStepLength` and `-minEdgeLength` options according to your mesh cell size.

- Parallel run example: `mpirun -np 3 smoothMesh -centroidalIters 20 -maxStepLength 0.01 -minEdgeLength 0.05 -parallel`
- Serial run example: `smoothMesh -centroidalIters 20 -maxStepLength 0.01 -minEdgeLength 0.05`

## Test case

The folder `testcase` contains an artificial test case which contains
skewed and non-orthogonal cells, as well as variance in geometric
shell shapes and topology. This is meant to be a challenging (but not
impossible) task for centroidal smoothing.

## Getting help and feedback

Please use Github issues section for asking help. If you like this
tool, please star this repository in Github!
