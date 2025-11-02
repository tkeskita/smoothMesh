/*---------------------------------------------------------------------------*\
Library
    Boundary Point Smoothing

Description
    Projection of boundary points to feature edges and boundary surfaces
\*---------------------------------------------------------------------------*/

#include "fvMesh.H"
#include "smoothMeshCommon.H"

using namespace Foam;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Help function to check the sanity of the edge mesh. Checks minimum
// edge length and perimeter size in comparison with polyMesh.

int checkEdgeMeshSanity
(
    const edgeMesh& em,
    const double meshMinEdgeLength,
    const double meshPerimeter
)
{
    double minEdgeLength = VGREAT;
    double bbMinX = VGREAT;
    double bbMaxX = -VGREAT;
    double bbMinY = VGREAT;
    double bbMaxY = -VGREAT;
    double bbMinZ = VGREAT;
    double bbMaxZ = -VGREAT;

    forAll(em.edges(), edgeI)
    {
        const label startPointI = em.edges()[edgeI][0];
        const label endPointI = em.edges()[edgeI][1];

        const point startPoint = em.points()[startPointI];
        const point endPoint = em.points()[endPointI];

        const double edgeLength = mag(endPoint - startPoint);

        if (edgeLength < minEdgeLength)
        {
            minEdgeLength = edgeLength;
        }

        // Bounding box update
        if (startPoint[0] < bbMinX) { bbMinX = startPoint[0]; }
        if (startPoint[1] < bbMinY) { bbMinY = startPoint[1]; }
        if (startPoint[2] < bbMinZ) { bbMinZ = startPoint[2]; }
        if (startPoint[0] > bbMaxX) { bbMaxX = startPoint[0]; }
        if (startPoint[1] > bbMaxY) { bbMaxY = startPoint[1]; }
        if (startPoint[2] > bbMaxZ) { bbMaxZ = startPoint[2]; }

        if (endPoint[0] < bbMinX) { bbMinX = endPoint[0]; }
        if (endPoint[1] < bbMinY) { bbMinY = endPoint[1]; }
        if (endPoint[2] < bbMinZ) { bbMinZ = endPoint[2]; }
        if (endPoint[0] > bbMaxX) { bbMaxX = endPoint[0]; }
        if (endPoint[1] > bbMaxY) { bbMaxY = endPoint[1]; }
        if (endPoint[2] > bbMaxZ) { bbMaxZ = endPoint[2]; }
    }

    if (minEdgeLength < REL_TOL * meshMinEdgeLength)
    {
        FatalError << "Minimum edge length in edge mesh " << minEdgeLength << " is too small in comparison to minimum edge length in polyMesh " << meshMinEdgeLength << endl << abort(FatalError);
    }

    const double emPerimeter = bbMaxX - bbMinX + bbMaxY - bbMinY + bbMaxZ + bbMinZ;
    const double PERIMETER_TOLERANCE = 0.5;

    if (abs((emPerimeter / meshPerimeter) - 1.0) > PERIMETER_TOLERANCE)
    {
        FatalError << "Perimeter (sum of bounding box side lengths) of edge mesh " << emPerimeter << " is too different in comparison to perimeter of polyMesh " << meshPerimeter << endl << abort(FatalError);
    }

    return 0;
}


// Help function to project point coordinates to a given edge in an
// edge mesh. Does not allow extrapolation over the edge end points
// (clips projection at edge end points). Additionally saves edge mesh
// point index to edgePointI, but only if the freely projected point
// (without clipping at edge ends) coincides with the edge mesh point.

int projectPointToEdge
(
    const point& pt,
    const edgeMesh& em,
    const label edgeI,
    const double distanceTolerance,
    point& projPoint,
    label& edgePointI
)
{
    edgePointI = UNDEF_LABEL;

    const label startPointI = em.edges()[edgeI][0];
    const label endPointI = em.edges()[edgeI][1];

    const point startPoint = em.points()[startPointI];
    const point endPoint = em.points()[endPointI];

    const double edgeLength = mag(endPoint - startPoint);

    // Project on the edge
    const point c2pt = pt - startPoint;
    const point edgeVec = endPoint - startPoint;
    const double normalizedDotProd = (c2pt & edgeVec) / sqr(edgeLength);
    const point testProjPoint = startPoint + normalizedDotProd * edgeVec;

    //if (mag(pt - vector(1,1,1)) < 1e-6)
    //    Info << " edgeI " << edgeI << " startPoint " << startPoint <<  " endPoint " << endPoint << " c2pt " << c2pt << " edgeVec " << edgeVec << " normalizedDotProd " << normalizedDotProd << " projPoint " << startPoint + normalizedDotProd * edgeVec << endl;

    // Clip at start point
    if (normalizedDotProd <= ABS_TOL)
    {
        projPoint = startPoint;
        if (mag(testProjPoint - startPoint) <= distanceTolerance)
        {
            edgePointI = startPointI;
        }
    }

    // Clip at end point
    else if (normalizedDotProd >= (1.0 - ABS_TOL))
    {
        projPoint = endPoint;
        if (mag(testProjPoint - endPoint) <= distanceTolerance)
        {
            edgePointI = endPointI;
        }
    }

    // Point is on the edge
    else
    {
        projPoint = testProjPoint;
    }

    return 0;
}


// Help function to find the index of the closest corner point in an
// edge mesh.

label findClosestEdgeMeshCornerPointIndex
(
    const point pt,
    const edgeMesh& em
)
{
    scalar distance(GREAT);
    label closestPointI = UNDEF_LABEL;

    forAll(em.points(), pointI)
    {
        // Skip non-corner points
        if (em.pointEdges()[pointI].size() == 2)
        {
            continue;
        }

        const point c = em.points()[pointI];
        const double testDistance = mag(pt - c);

        if (testDistance < distance)
        {
            distance = testDistance;
            closestPointI = pointI;
        }
    }

    if (closestPointI == UNDEF_LABEL)
    {
        FatalError << "Did not find any eligible corner points in edge mesh" << endl << abort(FatalError);
    }

    return closestPointI;
}

// Help function to project an arbitrary point pt onto the closest
// edge in an edge mesh. Returns projected point coordinates and
// various information from the closest edge.
//
// Search of edges eligible for projection can be optionally limited
// to a required edge string index number requiredStringI. Value -1
// for requiredStringI means that all edges are included in the
// search.
//
// Output variables are:
//
// - projPoint: the point location on the closest edge
// - closestEdgeI: the closest edge mesh edge index
// - closestEdgeStringI: the string index of the closest edge.
//   If projected point is exactly on an edge end (e.g. at a
//   corner point), then largest string index is returned.
// - closestEdgePointI: index of edge mesh point index, if the
//   freely projected point (not clipped at end) coincides with an
//   edge mesh point

int findClosestEdgeInfo
(
    const point& pt,
    const edgeMesh& em,
    const label requiredStringI,
    const labelList& targetEdgeStrings,
    const double distanceTolerance,
    point& projPoint,
    label& closestEdgeI,
    label& closestEdgeStringI,
    label& closestEdgePointI
)
{
    scalar distance(GREAT);
    projPoint = UNDEF_VECTOR;
    closestEdgeI = UNDEF_LABEL;
    closestEdgeStringI = UNDEF_LABEL;
    closestEdgePointI = UNDEF_LABEL;

    forAll(em.edges(), edgeI)
    {
        // Skip edge if string index is wrong
        if ((requiredStringI >= 0) and (targetEdgeStrings[edgeI] != requiredStringI))
        {
            continue;
        }

        // Project to edge
        point testProjPoint;
        label edgePointI;
        projectPointToEdge(pt, em, edgeI, distanceTolerance, testProjPoint, edgePointI);
        const scalar testDistance = mag(testProjPoint - pt);

        // Save if projected point is closest 
        if (testDistance < distance)
        {
            distance = testDistance;
            projPoint = testProjPoint;
            closestEdgeI = edgeI;
            closestEdgePointI = edgePointI;

            // Only set closest edge string index if the edge mesh size
            // matches string list size. Assumption here is that edge
            // mesh em is same (or at least very similar) to the
            // target edge mesh, so indices are correct.
            if (em.edges().size() == targetEdgeStrings.size())
            {
                closestEdgeStringI = targetEdgeStrings[edgeI];
            }
        }
    }

    if ((requiredStringI >= 0) and (closestEdgeStringI == UNDEF_LABEL))
    {
        FatalError << "Internal sanity check failed: Did not find any edges with string index " << requiredStringI << endl << abort(FatalError);
    }

    return 0;
}

// Classify boundary points, find corner points and target locations
// for the corners using the provided edge meshes

int classifyBoundaryPoints
(
    const fvMesh& mesh,
    const edgeMesh& initEdges,
    const edgeMesh& targetEdges,
    const labelList& layerPatchIds,
    const labelList& smoothingPatchIds,
    const boolList& isInternalPoint,
    boolList& isProcessorPoint,
    boolList& isConnectedToInternalPoint,
    boolList& isFeatureEdgePoint,
    boolList& isCornerPoint,
    vectorList& cornerPoints,
    boolList& isLayerSurfacePoint,
    boolList& isSmoothingSurfacePoint,
    boolList& isFrozenSurfacePoint,
    const labelList& targetEdgeStrings,
    const bool doBoundarySmoothing,
    const double distanceTolerance
)
{
    label nFeatureEdgePoints = 0;
    label nCornerPoints = 0;
    label nLayerSurfacePoints = 0;
    label nSmoothingSurfacePoints = 0;
    label nFrozenSurfacePoints = 0;

    boolList isVisitedPoint(mesh.nPoints(), false);

    forAll (mesh.boundary(), patchI)
    {
        const polyPatch& pp = mesh.boundaryMesh()[patchI];

        const label startI = mesh.boundary()[patchI].start();
        const label endI = startI + mesh.boundary()[patchI].Cf().size();

        for (label faceI = startI; faceI < endI; faceI++)
        {
            const face& f = mesh.faces()[faceI];
            forAll (f, facePointI)
            {
                const label pointI = mesh.faces()[faceI][facePointI];

                // Process each point only once
                if (isVisitedPoint[pointI])
                    continue;
                isVisitedPoint[pointI] = true;

                // Processor patch
                if (isA<processorPolyPatch>(pp))
                {
                    isProcessorPoint[pointI] = true;
                }

                // Skip rest of classifications if this is not a boundary point
                if (isInternalPoint[pointI])
                {
                    continue;
                }

                // Check if boundary point has connections to internal mesh point
                forAll (mesh.pointPoints()[pointI], pointPointI)
                {
                    const label i = mesh.pointPoints()[pointI][pointPointI];
                    if (isInternalPoint[i])
                    {
                        isConnectedToInternalPoint[pointI] = true;
                    }
                }

                // Classification of boundary smoothing points is done
                // only if information is available
                if ((initEdges.points().size() > 0) and (targetEdges.points().size() > 0))
                {
                    const point pt = mesh.points()[pointI];
                    point projPoint;
                    label dummy, dummy2;
                    label closestEdgePointI = UNDEF_LABEL;
                    findClosestEdgeInfo(pt, initEdges, -1, targetEdgeStrings, distanceTolerance, projPoint, dummy, dummy2, closestEdgePointI);

                    //if (mag(pt - vector(1, 1, 1)) < 1e-6)
                    //    Info << "pointI " << pointI << " at " << pt << " projPoint" << projPoint << " closestEdgePointI " << closestEdgePointI << endl;

                    // Corner points
                    if ((closestEdgePointI >= 0) and
                        (initEdges.pointEdges()[closestEdgePointI].size() != 2))
                    {
                        isCornerPoint[pointI] = true;
                        nCornerPoints++;

                        // Find and save the target coordinate for this corner
                        const label closestCornerPointI = findClosestEdgeMeshCornerPointIndex(pt, targetEdges);
                        cornerPoints[pointI] = targetEdges.points()[closestCornerPointI];
                        // Info << mesh.points()[pointI][0] << "," << mesh.points()[pointI][1] << "," << mesh.points()[pointI][2] << endl;
                    }

                    // Feature edges
                    else if (mag(pt - projPoint) < distanceTolerance)
                    {
                        isFeatureEdgePoint[pointI] = true;
                        nFeatureEdgePoints++;
                    }
                }

                // Layer treatment surface point
                const bool testIsLayerSurfacePoint = (findIndex(layerPatchIds, patchI) >= 0);
                if (testIsLayerSurfacePoint)
                {
                    isLayerSurfacePoint[pointI] = true;
                    nLayerSurfacePoints++;
                }

                // Smoothing surface points
                const bool testIsSmoothingSurfacePoint = (findIndex(smoothingPatchIds, patchI) >= 0);

                if ((doBoundarySmoothing) and (testIsSmoothingSurfacePoint))
                {
                    isSmoothingSurfacePoint[pointI] = true;
                    nSmoothingSurfacePoints++;
                    continue;
                }

                // Frozen boundary points are all the non-smoothed points
                else
                {
                    isFrozenSurfacePoint[pointI] = true;
                    nFrozenSurfacePoints++;
                    continue;
                }
            }
        }
    }

    // Summarize
    const label nSumCornerPoints  = returnReduce(nCornerPoints, sumOp<label>());
    const label nSumFeatureEdgePoints  = returnReduce(nFeatureEdgePoints, sumOp<label>());
    const label nSumLayerSurfacePoints  = returnReduce(nLayerSurfacePoints, sumOp<label>());
    const label nSumSmoothingSurfacePoints  = returnReduce(nSmoothingSurfacePoints, sumOp<label>());
    const label nSumFrozenSurfacePoints  = returnReduce(nFrozenSurfacePoints, sumOp<label>());

    Info << "Boundary point classification summary:" << endl;
    Info << "- Detected number of corner points: " << nSumCornerPoints << endl;
    Info << "- Detected number of feature edge points: " << nSumFeatureEdgePoints << endl;
    Info << "- Detected number of layer surface points: " << nSumLayerSurfacePoints << endl;
    Info << "- Detected number of smoothing surface points: " << nSumSmoothingSurfacePoints << endl;
    Info << "- Detected number of frozen surface points: " << nSumFrozenSurfacePoints << endl;
    Info << endl;

    return 0;
}

// Help function to get indices of continous edge mesh edges connected
// to a given edge.

int findContinuousEdgeMeshEdges
(
    const edgeMesh& em,
    const label edgeI,
    label& neighEdgeI1,
    label& neighEdgeI2
)
{
    neighEdgeI1 = UNDEF_LABEL;
    neighEdgeI2 = UNDEF_LABEL;

    // Process first end point
    const label pointI1 = em.edges()[edgeI][0];
    const label nEdges1 = em.pointEdges()[pointI1].size();

    if (nEdges1 == 2)
    {
        label edgeI1 = em.pointEdges()[pointI1][0];
        if (edgeI1 == edgeI)
        {
            edgeI1 = em.pointEdges()[pointI1][1];
        }
        neighEdgeI1 = edgeI1;
    }

    // Process second end point
    const label pointI2 = em.edges()[edgeI][1];
    const label nEdges2 = em.pointEdges()[pointI2].size();

    if (nEdges2 == 2)
    {
        label edgeI2 = em.pointEdges()[pointI2][0];
        if (edgeI2 == edgeI)
        {
            edgeI2 = em.pointEdges()[pointI2][1];
        }
        neighEdgeI2 = edgeI2;
    }

    return 0;
}

// Help function to assign a string index number to current edge mesh
// edge and recursively assing the same index to all neighbor edges in
// the string

int stringifyEdgeMeshEdges
(
    const edgeMesh& em,
    labelList& targetEdgeStrings,
    const label edgeI,
    const label neighEdgeI1,
    const label neighEdgeI2,
    label& nStrings
)
{
    // Find string indices of this edge and neighbor edges
    const label stringI0 = targetEdgeStrings[edgeI];

    label stringI1 = UNDEF_LABEL;
    if (neighEdgeI1 != UNDEF_LABEL)
    {
        stringI1 = targetEdgeStrings[neighEdgeI1];
    }

    label stringI2 = UNDEF_LABEL;
    if (neighEdgeI2 != UNDEF_LABEL)
    {
        stringI2 = targetEdgeStrings[neighEdgeI2];
    }

    const label maxStringI = max(max(stringI0, stringI1), stringI2);

    // Neighbors or self has no string index, generate a new value
    if (maxStringI == UNDEF_LABEL)
    {
        ++nStrings;
        targetEdgeStrings[edgeI] = nStrings;
    }

    // Neighbor or self already has a string value assigned
    else if (stringI0 == UNDEF_LABEL)
    {
        targetEdgeStrings[edgeI] = maxStringI;
    }

    // Recurse into neighbour edges if they don't already have a
    // string index and if the neighbor is not a corner
    if ((neighEdgeI1 != UNDEF_LABEL) and (stringI1 == UNDEF_LABEL) and (em.pointEdges()[neighEdgeI1].size() == 2))
    {
        label neighNeighEdgeI1(UNDEF_LABEL);
        label neighNeighEdgeI2(UNDEF_LABEL);
        findContinuousEdgeMeshEdges(em, neighEdgeI1, neighNeighEdgeI1, neighNeighEdgeI2);
        stringifyEdgeMeshEdges(em, targetEdgeStrings, neighEdgeI1, neighNeighEdgeI1, neighNeighEdgeI2, nStrings);
    }

    if ((neighEdgeI2 != UNDEF_LABEL) and (stringI2 == UNDEF_LABEL) and (em.pointEdges()[neighEdgeI2].size() == 2))
    {
        label neighNeighEdgeI1(UNDEF_LABEL);
        label neighNeighEdgeI2(UNDEF_LABEL);
        findContinuousEdgeMeshEdges(em, neighEdgeI2, neighNeighEdgeI1, neighNeighEdgeI2);
        stringifyEdgeMeshEdges(em, targetEdgeStrings, neighEdgeI2, neighNeighEdgeI1, neighNeighEdgeI2, nStrings);
    }

    return 0;
}


// Primary help function to start the identification and labeling of
// continuous edge strings in an edge mesh

int findEdgeMeshStrings
(
    labelList& targetEdgeStrings,
    const edgeMesh& em
)
{
    // Unique storage for next available string index number
    label nStrings(UNDEF_LABEL);

    // Initialize string list
    targetEdgeStrings.resize(em.edges().size());
    forAll (em.edges(), edgeI)
    {
        targetEdgeStrings[edgeI] = UNDEF_LABEL;
    }

    // Process all edge mesh edges
    forAll (em.edges(), edgeI)
    {
        // Skip edge if string index is already assigned
        if (targetEdgeStrings[edgeI] >= 0)
            continue;

        label neighEdgeI1(UNDEF_LABEL);
        label neighEdgeI2(UNDEF_LABEL);
        findContinuousEdgeMeshEdges(em, edgeI, neighEdgeI1, neighEdgeI2);
        stringifyEdgeMeshEdges(em, targetEdgeStrings, edgeI, neighEdgeI1, neighEdgeI2, nStrings);
    }

    return nStrings;
}

// Help function to find indices of non-corner and non-feature-edge
// boundary points next to given boundary point

labelList findNeighborSurfacePoints
(
    const fvMesh& mesh,
    const label pointI,
    const boolList& isInternalPoint,
    const boolList& isFeatureEdgePoint,
    const boolList& isCornerPoint
)
{
    labelList neighIs;

    forAll (mesh.pointPoints()[pointI], pointPointI)
    {
        const label neighI = mesh.pointPoints()[pointI][pointPointI];
        if (isInternalPoint[neighI])
            continue;
        if (isFeatureEdgePoint[neighI])
            continue;
        if (isCornerPoint[neighI])
            continue;
        neighIs.append(neighI);
    }

    return neighIs;
}


// Calculate new feature edge point positions and synchronizes among
// processors. Projects neighboring surface mesh points to closest
// feature edges and uses the median as a new feature edge point.

int calculateFeatureEdgeProjections
(
    const fvMesh& mesh,
    const edgeMesh& em,
    const boolList& isInternalPoint,
    const boolList& isFeatureEdgePoint,
    const boolList& isCornerPoint,
    const labelList& targetEdgeStrings,
    const labelList& pointStrings,
    const double distanceTolerance,
    vectorList& featureEdgeProjections,
    labelList& nFeatureEdgeProjections
)
{
    forAll(mesh.points(), pointI)
    {
        // This is done only for feature edge points
        if (! isFeatureEdgePoint[pointI])
            continue;

        // Find the valid surface neighbor points
        const labelList neighPointIs = findNeighborSurfacePoints(mesh, pointI, isInternalPoint, isFeatureEdgePoint, isCornerPoint);

        // Project neighbor points to target edges and add projected point
        forAll (neighPointIs, neighPointI)
        {
            const label pointStringI = pointStrings[pointI];
            const point neighPoint = mesh.points()[neighPointIs[neighPointI]];
            point projPoint;
            label dummy, dummy2, dummy3;
            findClosestEdgeInfo(neighPoint, em, pointStringI, targetEdgeStrings, distanceTolerance, projPoint, dummy, dummy2, dummy3);
            featureEdgeProjections[pointI] += projPoint;
            ++nFeatureEdgeProjections[pointI];
        }
    }

    // Synchronize among processors (using sum combination)
    syncTools::syncPointList
    (
        mesh,
        featureEdgeProjections,
        plusEqOp<vector>(),
        ZERO_VECTOR               // null value
    );

    syncTools::syncPointList
    (
        mesh,
        nFeatureEdgeProjections,
        plusEqOp<label>(),
        0                         // null value
    );

    return 0;
}


// Find intersection with triSurface faces in point normal direction

point findIntersection
(
    const indexedOctree<treeDataTriSurface>& tree,
    const label pointI,
    const point origPoint,
    const vector pointNormal,
    const double searchDistance
)
{
    if (pointNormal == ZERO_VECTOR)
    {
        FatalError
            << "pointNormal is zero for pointI " << pointI << " at " << origPoint << endl
            << abort(FatalError);
    }

    // Search towards normal direction
    vector hitPoint1(UNDEF_VECTOR);
    {
        const point endPoint = origPoint + pointNormal * searchDistance;
        const pointIndexHit hitInfo = tree.findLine(origPoint, endPoint);
        if (hitInfo.hit())
        {
            hitPoint1 = hitInfo.hitPoint();
        }
    }

    // Search towards opposite direction
    vector hitPoint2(UNDEF_VECTOR);
    {
        const point endPoint = origPoint - pointNormal * searchDistance;
        const pointIndexHit hitInfo = tree.findLine(origPoint, endPoint);
        if (hitInfo.hit())
        {
            hitPoint2 = hitInfo.hitPoint();
        }
    }

    // Return the closest hit point
    const double distance1 = mag(origPoint - hitPoint1);
    const double distance2 = mag(origPoint - hitPoint2);
    if (distance1 < distance2)
    {
        return hitPoint1;
    }
    else if (distance2 < distance1)
    {
        return hitPoint2;
    }

    // Otherwise search in between
    {
        const point endPoint1 = origPoint + pointNormal * searchDistance;
        const point endPoint2 = origPoint - pointNormal * searchDistance;
        const pointIndexHit hitInfo = tree.findLine(endPoint1, endPoint2);
        if (hitInfo.hit())
        {
            return hitInfo.hitPoint();
        }
    }

    return UNDEF_VECTOR;
}

// Help function to retrieve boundary face center for a given face
// index. Use of mesh.Cf()[X] for X larger than number of internal
// faces raises "index out of range" error in Debug build of OpenFOAM
// (although it does seem to work for Opt build of OpenFOAM).
// Access to boundary faces is done via patches.

const vector getBoundaryFaceCf
(
    const fvMesh& mesh,
    const label faceI
)
{
    forAll(mesh.boundaryMesh(), patchI)
    {
        const fvPatch& patch = mesh.boundary()[patchI];
        const label startI = patch.start();
        const label endI = startI + patch.size();

        if (faceI < endI)
        {
            const label patchFaceI = faceI - startI;
            const vector Cf = patch.Cf()[patchFaceI];
            return Cf;
        }
    }

    // Sanity check
    FatalError << "Could not find Cf for face index " << faceI << endl << abort(FatalError);

    return UNDEF_VECTOR;
}

// Help function to calculate centroid of boundary point boundary face
// centers

int calculateSurfaceCentroids
(
    const fvMesh& mesh,
    const boolList& isInternalPoint,
    vectorList& faceCentroids,
    labelList& nFaceCentroids
)
{
    forAll(faceCentroids, i)
    {
        faceCentroids[i] = Zero;
        nFaceCentroids[i] = 0;
    }

    const label firstBoundaryFaceI = mesh.boundary()[0].start();

    forAll(mesh.points(), pointI)
    {
        // Only consider boundary points
        if (isInternalPoint[pointI])
            continue;

        forAll(mesh.pointFaces()[pointI], pointFaceI)
        {
            const label faceI = mesh.pointFaces()[pointI][pointFaceI];

            // Ignore other than boundary faces
            if (faceI < firstBoundaryFaceI)
                continue;

            faceCentroids[pointI] += getBoundaryFaceCf(mesh, faceI);
            ++nFaceCentroids[pointI];
        }

        if (nFaceCentroids[pointI] == 0)
        {
            FatalError << "did not find faceNeighbour for point " << pointI << endl << abort(FatalError);
        }
    }

    // Synchronize among processors (using sum combination)
    syncTools::syncPointList
    (
        mesh,
        faceCentroids,
        plusEqOp<vector>(),
        ZERO_VECTOR               // null value
    );

    syncTools::syncPointList
    (
        mesh,
        nFaceCentroids,
        plusEqOp<label>(),
        0                         // null value
    );

    return 0;
}

// Projection of boundary points to feature edges and boundary surfaces

int projectBoundaryPointsToEdgesAndSurfaces
(
    const fvMesh& mesh,
    pointField& newPoints,
    const pointField& pointNormals,
    const boolList& isInternalPoint,
    const boolList& isSmoothingSurfacePoint,
    const boolList& isFeatureEdgePoint,
    const boolList& isCornerPoint,
    const vectorList& cornerPoints,
    const edgeMesh& targetEdges,
    const triSurface& surf,
    const indexedOctree<treeDataTriSurface>& tree,
    const double meshMinEdgeLength,
    const labelList& targetEdgeStrings,
    const labelList& pointStrings,
    const boolList& isSharpEdgePoint,
    const double distanceTolerance,
    boolList& isFrozenPoint
)
{
    // Calculate new feature edge point positions a priori
    vectorList featureEdgeProjections(mesh.nPoints(), ZERO_VECTOR);
    labelList nFeatureEdgeProjections(mesh.nPoints(), 0);
    calculateFeatureEdgeProjections(mesh, targetEdges, isInternalPoint, isFeatureEdgePoint, isCornerPoint, targetEdgeStrings, pointStrings, distanceTolerance, featureEdgeProjections, nFeatureEdgeProjections);

    // Calculate boundary surface face centroids a priori
    vectorList faceCentroids(mesh.nPoints(), ZERO_VECTOR);
    labelList nFaceCentroids(mesh.nPoints(), 0);
    calculateSurfaceCentroids(mesh, isInternalPoint, faceCentroids, nFaceCentroids);

    // Blending fraction of face centroidal
    // TODO: Needs more testing to understand how to make this stable.
    const double faceCentroidBlendingFraction = 0.0;

    forAll(mesh.points(), pointI)
    {
        if (isInternalPoint[pointI])
            continue;

        // Project to closest corner points
        // --------------------------------
        if (isCornerPoint[pointI])
        {
            newPoints[pointI] = cornerPoints[pointI];
            continue;
        }

        // Project to closest feature edges
        // --------------------------------
        if (isFeatureEdgePoint[pointI])
        {
            // Project neighboring surface mesh points to closest
            // feature edges and use the median as a new feature edge point

            const vector newPoint = featureEdgeProjections[pointI] / double(nFeatureEdgeProjections[pointI]);
            newPoints[pointI] = newPoint;
            continue;
        }

        // Freeze very sharp edge points which are not feature edge points
        if (isSharpEdgePoint[pointI])
        {
            isFrozenPoint[pointI] = true;
        }

        // Project to closest tri face in normal or opposite direction
        // -----------------------------------------------------------
        else if (isSmoothingSurfacePoint[pointI])
        {
            const point pointNormal = pointNormals[pointI];
            double searchDistance = distanceTolerance;

            // Blend face centroid with cell centroid
            const point newPoint = faceCentroidBlendingFraction * (faceCentroids[pointI] / double(nFaceCentroids[pointI])) + (1 - faceCentroidBlendingFraction) * newPoints[pointI];

            // Project new point to surface. Intersection search is
            // done with increasing search distance for accuracy.
            point surfPoint = UNDEF_VECTOR;
            for (label i = 0; i < 4; ++i)
            {
                searchDistance *= (1.0 / REL_TOL);
                surfPoint = findIntersection(tree, pointI, newPoint, pointNormal, searchDistance);
                if (surfPoint != UNDEF_VECTOR)
                {
                    newPoints[pointI] = surfPoint;
                    break;
                }
            }

            if (surfPoint == UNDEF_VECTOR)
            {
                FatalError
                    << "Did not find surface intersection for pointI " << pointI
                    << " at " << newPoint << " pointNormal " << pointNormal
                    << " searchDistance " << searchDistance << endl
                    << abort(FatalError);
            }
        }
    }

    return 0;
}

// ************************************************************************* //
