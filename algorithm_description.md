# SmoothMesh Algorithm Description

SmoothMesh relies on two algorithms for mesh smoothing: **Centroidal
smoothing** to *provide new target coordinates for moving each mesh
point*, and **heuristic quality constraints** to *disallow movement if
mesh quality would suffer too much*.

Centroidal smoothing is a version of the
[Laplacian smoothing algorithm](https://en.wikipedia.org/wiki/Laplacian_smoothing),
which uses surrounding cell centers instead of the neighbouring point
locations to calculate the new target position for mesh points. The
jump vectors (from current internal point positions to new
coordinates) calculated with the centroidal smoothing algorithm are
taken as a starting point. The absolute length of the jump done in
each iteration is limited to a user given maximum step length value
(`-maxStepLength` option), to allow stable iteration during the
smoothing process. Boundary points are not allowed to move, to keep
the mesh geometry intact.

Unconstrained centroidal smoothing works well as such for cases where
the initial mesh quality w.r.t face non-orthogonality, face skewness
and cell aspect ratio is moderate. However, if the mesh contains cells
with low quality, then centroidal smoothing can produce either
self-intersecting or boundary-intersecting cells. Therefore SmoothMesh
applies optional **heuristic quality constraints** which constrain
point movement, to avoid mesh issues. All of the quality constraints
work by prohibiting point movement, effectively "freezing" a point to
it's current coordinates, if quality would be decreased too much with
moving of the point.

The quality constraints include the following, and they are by default
all evaluated in sequence:

## 1. Avoid shortening of short edge length

This constraint considers only **the minimum length of edges connected
to an internal point**. For each moving point, the minimum length of
all edges of the point are calculated first from current mesh, and
then assuming the point is moved to new coordinates. If new length is
below allowed minimum edge length (specified with the `-minEdgeLength`
option) and if the new minimum length is smaller than the current
minimum length, then the point is frozen.

This constraint is useful in the case where cells contain skewed faces
near boundaries, which tend to shrink the edge length in centroidal
smoothing, like in the figure below.

<p align="left"><img src="images/base_mesh_with_problematic_vertex.png"></p>


## 2. Restrict decrease of smallest edge-edge angle

This constraint deals only with **the minimum angle between edges
connected to an internal point**. For each point, the minimum
angle is calculated for all edge pairs which share a face and the
point. First, the current minimum angle is calculated using the
current mesh point locations. Then, a minimum new angle is calculated
using new coordinates for points. Three variations are tested in the
calculation: New position of the center point, new position for one
edge end point, and new position for the other edge position. Minimum
of those is selected as the new minimum angle. If the new minimum angle
is below allowed minimum edge angle (specified with the
`-minEdgeAngle` option), and if the new minimum angle is smaller than the
current minimum angle, then the point is frozen.

TODO: Add the neighbour point freezing functionality to edge-edge
angle checking, like done for face-face angles?

While the smallest edge-edge angle is helpful in catching many bad
cell shapes, it is not enough in the cases where cells are flattened in
such a way that edges don't collapse on top of each other, like the
example cell in the figure below (top view on left, front view on
right).

<p align="left"><img src="images/base_mesh_with_problematic_vertex.png"></p>


## 3. Restrict deterioration of face-face angles

This constraint focuses only on **the minimum and maximum angle
between _cell faces_, which are connected to the all of the edges,
which surround a single mesh point**. For all mesh points (not just
internal points!), each point is processed as follows, please refer to
the example figure below:

- Minimum and maximum face angle is calculated for the point *p0* edges
  using the current mesh point locations in the
  calculation. Calculation procedure is described below in more
  detail.

- The main point *p0* is hypothetically moved to it's new
  coordinates. The effect of the move to the minimum and maximum face
  angles at *p0* is recalculated using new *p0* coordinates. If minimum
  or maximum face angle is deteriorated compared to current values,
  then *p0* is frozen.

- Each of the *p0* neighbour points (*p1*, *p2* and *p6*) are
  hypothetically moved to their new coordinates, one at a time. The
  effect of the move to the minimum and maximum face angles at *p0* is
  calculated using new coordinates of the neighbour points. If minimum
  or maximum face angle at *p0* is deteriorated compared to current
  values, then *that neighbour point* is frozen.

The deterioration of face angles is defined as: Minimum face angle is
below a threshold value (currently 45 deg), and the minimum angle
would decrease if the move would be allowed. Similary for maximum
angle: Maximum angle is above threshold (currently 170 deg), and the
maximum angle would increase if the move would be allowed.

### The calculation method for face angles

Shortly put, the calculation of face angles at an edge is done by
**first projecting face and cell centers to an edge normal plane**,
and then **calculating and summing the angles from edge center to each
projected point**.

Figure below illustrates an example for calculation of face angles for
a simple two cell case. For the sake of simplicity, the front and back
faces are not considered. 

- Point *p0* in the figure is the point for which face angles are to
  be calculated, and the figure illustrates the calculation of face
  angles for the edge *p0-p1*. Edge center point is *d0*.

- There are two faces for this edge: First face *p0-p6-p7-p1* with
  face center at *c4*, and a second face *p0-p1-p3-p2* with face center
  at *c3*.

- The cell centers for the upper cell is *c2* and for the lower cell
  *c1*.

- Edge normal plane is defined as a plane located at *d0* with normal
  direction *p1-p0*. Each of the center points *c* are projected onto
  this plane, resulting in corresponding *d* points.

- Face angles are calculated from the projected *d* points by summing
  the angles between projected face and cell center coordinates
  with the edge center coordinates *d0*:

  - The face angle for lower cell is the sum of angles *a1* (angle from
    *d3-d0-d1*) and *a2* (angle from *d1-d0-d4*).

  - Similarly, the face angle for upper cell is the sum of angles *a3* (angle from
    *d4-d0-d2*) and *a4* (angle from *d2-d0-d3*).

The reason for using cell centers in angle calculation is that it
allows better evaluation of concave angles (>180 deg face angles),
like *a1+a2* in the figure.

<p align="left"><img src="images/face_angle.png"></p>
