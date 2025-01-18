# smoothMesh

<img src="images/smooth_mesh.png" width="600"/>

OpenFOAM mesh smoothing tool to improve mesh quality. Moves internal
mesh points by using the Centroidal smoothing algorithm (a version of the
[Laplacian smoothing algorithm](https://en.wikipedia.org/wiki/Laplacian_smoothing),
which uses surrounding cell centers instead of the neighbour point
locations to calculate the new point position). Optional heuristic
quality constraint options exist to constrain the smoothing, to avoid
self-intersections. No changes to mesh topology are made.

Image below illustrates the need for restricting centroidal
smoothing. Without quality constraints, centroidal smoothing would
move the point highlighted in blue to the location highlighted with
green, which is outside of the domain. Thereby, centroidal smoothing
can create self-intersecting cells, depending on the geometry and
topology of the mesh. Self-intersections can be avoided by using
additional quality constraints, which restrict the movement of
vertices.

<img src="images/base_mesh_with_problematic_vertex.png" width="600"/>

## Current features and restrictions

- Works on 3D polyhedron meshes
- Requires a consistent (not self-intersecting or tangled) initial mesh with "good enough" quality
- Smoothes internal mesh points only (boundary points are frozen)
- Optional handling of prismatic boundary layers
- Developed on OpenFOAM.com v2312, tested on v2412


## Compilation instructions

```
. /usr/lib/openfoam/openfoam2412/etc/bashrc
cd smoothMesh/src
wclean; wmake
```


## Command line options

### Basic options

- `-centroidalIters` specifies the number of smoothing iterations (default 20).

- `-maxStepLength` is the maximum length (in metres) for moving a point in one iteration (default 0.01). Adjust this value for your case. Smoothing process seems to be stable when this value is in the range 10% - 50% of the minimum cell side length.

- `-minEdgeLength` defines edge length below which edge points are fully frozen at their current location. Freezing happens only if edge length would decrease during smoothing (default 0.05). Edge length is allowed to increase regardless of this value. Adjust this value for your case.

- `-totalMinFreeze` option causes mesh points on all edges shorter than `-minEdgeLength` to freeze, even if edge length would increase in smoothing (default false). This option is useful to keep boundary layers in the mesh unmodified, and smooth the large cells only, if the special boundary layer related options below are not used.

### Quality constraint options

The following options are related to additional **heuristic quality control constraints for smoothing**. The constraints work by disallowing movement of point (freezing of points) if the movement would cause quality of the mesh to suffer too much. Without constraining, centroidal smoothing may squish cells and create self-intersecting cells e.g. near concave geometry features, depending on the mesh details. Have a look at [the algorithm description document](algorithm_description.md) for details.

**Note:** The old `-qualityControl` option has been superceded by the options below.

- `-faceAngleConstraint` boolean option enables an additional quality control which restricts decrease of smallest and largest face-face angle (default is true).

- `-minAngle` option defines the value for minimum angle (in degrees, default value 35).

- `-maxAngle` option specifies the value for maximum angle (in degrees, default value 160).

- Note: `-minAngle` value causes point freezing *only* if the minimum angle is below this value and if the minimum angle would *decrease* in smoothing. Points are allowed to move if the minimum angle value *increases* with smoothing, regardless of this value. The same applies for the `-maxAngle` option: Freezing takes place only if maximum angle is above the specified value and if the maximum angle would *increase* in smoothing.

### Boundary layer related options

The options below are related to handling of prismatic cells near mesh boundaries, to either preserve or improve the orthogonality and the thickness of boundary layer cells in the mesh. If the mesh contains prismatic boundary layers, the unconstrained centroidal smoothing will tend to bloat the boundary layer cells into normal size. That can be avoided using the options below. These options affect only the prismatic cell edges near the mesh boundaries (see [the algorithm description document](algorithm_description.md) for details).

Warning: This is an experimental feature!

- `-boundaryMaxBlendingFraction` is the maximum fraction (0 <= value <= 1) by which boundary layer edge length and edge direction are blended with the centroidal smoothing locations. Zero value disables the effect of all other boundary related variables below (default 0). Value 0.8 seems to produce good results in practice.

- `-boundaryEdgeLength` specifies the target thickness for the first boundary layer cells (prismatic side edge length) (default: 0.05).

- `-boundaryExpansionRatio` specifies the thickness ratio by which the boundary edge length is assumed to increase (default: 1.3).

- `-boundaryMinLayers` is an integer value specifying the number of boundary layers which experience a full force of boundary blending specified with the `-boundaryMaxBlendingFraction` option (default: 1).

- `-boundaryMaxLayers` specifies the number of boundary cell layers beyond which boundary blending options above ceases to affect smoothing, and only centroidal smoothing is applied (default: 4).

- `-patches` option can be used to limit the boundary layer treatment to specified patches only. You can specify one or several patches, optionally with wild cards. For example `-patches 'walls'` or `-patches '( stator "rotor.*" )'`.


## Description of the algorithm

Please view [the algorithm description document](algorithm_description.md).


## Basic usage examples

Adjust at least the `-centroidalIters`, `-maxStepLength` and `-minEdgeLength` options according to your case.

- Parallel run example: `mpirun -np 3 smoothMesh -centroidalIters 20 -maxStepLength 0.01 -minEdgeLength 0.05 -parallel`

- Serial run example: `smoothMesh -centroidalIters 20 -maxStepLength 0.01 -minEdgeLength 0.05`


## Test case

The folder `testcase` contains an artificial test case which contains
skewed and non-orthogonal cells, as well as variance in geometric
cell shapes and topology. This is meant to be a challenging (but not
impossible) task for centroidal smoothing.

Figures below illustrate how smoothing can go wrong (or right),
depending on the parameters applied.

This is the starting mesh:

<img src="images/testcase_00_original_mesh.png" width="600"/>

Test 1 (bad results): No boundary layers, without faceAngleConstraint, produces self-intersections.
Full command: `smoothMesh -centroidalIters 100 -minEdgeLength 0.01 -maxStepLength 0.004 -faceAngleConstraint false`

<img src="images/testcase_01_self_intersections.png" width="600"/>

Test 2 (bad results): No boundary layers, too large minAngle does not allow much smoothing.
Full command: `smoothMesh -centroidalIters 100 -minEdgeLength 0.01 -maxStepLength 0.004 -minAngle 45 -maxAngle 160 -faceAngleConstraint true`

<img src="images/testcase_02_large_minAngle.png" width="600"/>

Test 3 (good results): No boundary layers, faceAngleConstraint creates mesh without self-intersections.
Full command: `smoothMesh -centroidalIters 100 -minEdgeLength 0.01 -maxStepLength 0.004 -minAngle 15 -maxAngle 160 -faceAngleConstraint true`

<img src="images/testcase_03_small_minAngle.png" width="600"/>

Test 4 (bad results): Boundary layers without patches specification creates boundary layers on outer walls.
Full command: `smoothMesh -centroidalIters 100 -minEdgeLength 0.01 -maxStepLength 0.004 -minAngle 15 -maxAngle 160 -faceAngleConstraint true -boundaryMaxBlendingFraction 0.8 -boundaryEdgeLength 0.01`

<img src="images/testcase_04_layers_without_patches_specified.png" width="600"/>

Test 5 (good results): Boundary layers for patch named default, best result.
Full command: `smoothMesh -centroidalIters 100 -minEdgeLength 0.01 -maxStepLength 0.004 -minAngle 15 -maxAngle 160 -faceAngleConstraint true -boundaryMaxBlendingFraction 0.8 -boundaryEdgeLength 0.01 -patches '("def.*")'`

<img src="images/testcase_05_layers_with_patches.png" width="600"/>


## Getting help and feedback

Please use Github issues section for asking help. If you like this
tool, please star this repository in Github!
