/*---------------------------------------------------------------------------*\
Library
    Boundary Point Smoothing

Description
    Projection of boundary points to feature edges and boundary surfaces
\*---------------------------------------------------------------------------*/

#include "fvMesh.H"

// Macros for value definitions
#define UNDEF_LABEL -1
#define UNDEF_VECTOR vector(GREAT, GREAT, GREAT)
#define ZERO_VECTOR vector(0, 0, 0)
#define ZERO 0.0

// Minimum distance for detection of points by coordinates.
// The value is fairly large to allow single precision point data from
// third party applications.
#define DISTANCE_TOLERANCE 1e-6

using namespace Foam;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Calculate and return the edge-edge angle (in radians,
// 0 < angle < pi) of two edges which share a common point
// at coordinate cCoords. The two edge end point coordinates are
// p1Coords and p2Coords.

// TODO: Find out how to import edgeEdgeAngle from helpFunctions.H and
// get rid of this copy

double edgeEdgeAngle2
(
    const vector cCoords,
    const vector p1Coords,
    const vector p2Coords
)
{
    vector vec1 = (p1Coords - cCoords);
    vector vec2 = (p2Coords - cCoords);
    vec1 /= mag(vec1);
    vec2 /= mag(vec2);

    const double cosA = vec1 & vec2;

    // Ensure cos angle is in sane range before calling arc cos
    const double MAX = 0.99999;
    const double cosAlpha = std::max(-MAX, std::min(MAX, cosA));
    const double angle = std::acos(cosAlpha);

    return angle;
}



// Find index of closest point in an edge mesh. Search may be limited
// to corner points only and/or to required string index number. Saves
// the closest edge mesh index to closestI and the edge string index
// of that point to closestStringI.

int findClosestEdgeMeshPointIndex
(
    const point pt,
    const edgeMesh& em,
    const bool mustBeCorner,
    const bool mustNotBeCorner,
    const bool matchTargetEdgeStrings,
    const label requiredStringI,
    const labelList& targetEdgeStrings,
    label& closestI,
    label& closestStringI
)
{
    // TODO: Replace linear search with an octree search for efficiency

    scalar distance = GREAT;
    closestI = UNDEF_LABEL;
    closestStringI = UNDEF_LABEL;

    forAll(em.points(), pointI)
    {
        // Disregard non-corner points if it's required
        if ((mustBeCorner) and (em.pointEdges()[pointI].size() == 2))
        {
            continue;
        }

        // Disregard corner point if it's required
        if ((mustNotBeCorner) and (em.pointEdges()[pointI].size() != 2))
        {
            continue;
        }

        // Disregard other edge strings if matching string index is
        // required
        if ((matchTargetEdgeStrings) and (targetEdgeStrings[pointI] != requiredStringI))
        {
            continue;
        }

        const point p = em.points()[pointI];
        const scalar testDistance = mag(pt - p);
        if (testDistance < distance)
        {
            distance = testDistance;
            closestI = pointI;
            closestStringI = targetEdgeStrings[pointI];
        }
    }

    if (closestI == UNDEF_LABEL)
    {
        FatalError << "Internal sanity check failed: Did not find any closest edge mesh points near " << pt << endl << abort(FatalError);
    }

    if ((mustNotBeCorner) and (closestStringI == UNDEF_LABEL))
    {
        FatalError << "Internal sanity check failed: Did not find string index for edge mesh points near " << pt << endl << abort(FatalError);
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
    const labelList& targetEdgeStrings
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
                    label closestInitPointI = UNDEF_LABEL;
                    label dummy;
                    findClosestEdgeMeshPointIndex(pt, initEdges, false, false, false, UNDEF_LABEL, targetEdgeStrings, closestInitPointI, dummy);
                    const point closestInitPoint = initEdges.points()[closestInitPointI];

                    // Corner points
                    if ((initEdges.pointEdges()[closestInitPointI].size() != 2)
                        and (mag(pt - closestInitPoint) < DISTANCE_TOLERANCE))
                    {
                        isCornerPoint[pointI] = true;
                        label closestCornerPointI = UNDEF_LABEL;
                        findClosestEdgeMeshPointIndex(closestInitPoint, targetEdges, true, false, false, UNDEF_LABEL, targetEdgeStrings, closestCornerPointI, dummy);
                        cornerPoints[pointI] = targetEdges.points()[closestCornerPointI];
                        nCornerPoints++;
                    }

                    // Feature edges
                    else if (mag(pt - closestInitPoint) < DISTANCE_TOLERANCE)
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

                if (testIsSmoothingSurfacePoint)
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

// Help function to find unique (no corner points) neighbour points of
// an edge mesh point

int findNeighborEdgeMeshPoints
(
    const edgeMesh& em,
    const label pointI,
    label& neighI1,
    label& neighI2
)
{
    neighI1 = UNDEF_LABEL;
    neighI2 = UNDEF_LABEL;

    const label nEdges = em.pointEdges()[pointI].size();

    // Do nothing if this point has more than two connected
    // edges (e.g. corner point)
    if (nEdges > 2)
    {
        return 0;
    }

    // Traverse first edge
    if (nEdges > 0)
    {
        const label edgeI1 = em.pointEdges()[pointI][0];
        label endPointI1 = em.edges()[edgeI1][0];
        if (endPointI1 == pointI)
        {
            endPointI1 = em.edges()[edgeI1][1];
        }
        neighI1 = endPointI1;
    }

    // Traverse second edge
    if (nEdges > 1)
    {
        const label edgeI2 = em.pointEdges()[pointI][1];
        label endPointI2 = em.edges()[edgeI2][0];
        if (endPointI2 == pointI)
        {
            endPointI2 = em.edges()[edgeI2][1];
        }
        neighI2 = endPointI2;
    }

    return 0;
}

// Help function to assign a string index number to current edge mesh
// point and recurse over all neighbor points in the string

int stringifyEdgeMeshPoints
(
    const edgeMesh& em,
    labelList& targetEdgeStrings,
    const label pointI,
    const label neighI1,
    const label neighI2,
    label& nStrings
)
{
    // Find string indices of this point and neighbour points
    const label stringI0 = targetEdgeStrings[pointI];

    label stringI1 = UNDEF_LABEL;
    if (neighI1 != UNDEF_LABEL)
    {
        stringI1 = targetEdgeStrings[neighI1];
    }

    label stringI2 = UNDEF_LABEL;
    if (neighI2 != UNDEF_LABEL)
    {
        stringI2 = targetEdgeStrings[neighI2];
    }

    const label maxStringI = max(max(stringI0, stringI1), stringI2);

    // New string, get new string index value
    if (maxStringI == UNDEF_LABEL)
    {
        ++nStrings;
        targetEdgeStrings[pointI] = nStrings;
    }
    // Neighbor or self already has a string value assigned
    else if (stringI0 == UNDEF_LABEL)
    {
        targetEdgeStrings[pointI] = maxStringI;
    }

    // Recurse into neighbour points if they don't already have a
    // string index and if the neighbor is not a corner
    if ((neighI1 != UNDEF_LABEL) and (stringI1 == UNDEF_LABEL) and (em.pointEdges()[neighI1].size() == 2))
    {
        label neighNeighI1(UNDEF_LABEL);
        label neighNeighI2(UNDEF_LABEL);
        findNeighborEdgeMeshPoints(em, neighI1, neighNeighI1, neighNeighI2);
        stringifyEdgeMeshPoints(em, targetEdgeStrings, neighI1, neighNeighI1, neighNeighI2, nStrings);
    }

    if ((neighI2 != UNDEF_LABEL) and (stringI2 == UNDEF_LABEL) and (em.pointEdges()[neighI2].size() == 2))
    {
        label neighNeighI1(UNDEF_LABEL);
        label neighNeighI2(UNDEF_LABEL);
        findNeighborEdgeMeshPoints(em, neighI2, neighNeighI1, neighNeighI2);
        stringifyEdgeMeshPoints(em, targetEdgeStrings, neighI2, neighNeighI1, neighNeighI2, nStrings);
    }

    return 0;
}


// Primary help function to start identification and labeling of
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
    targetEdgeStrings.resize(em.points().size());
    forAll (em.points(), pointI)
    {
        targetEdgeStrings[pointI] = UNDEF_LABEL;
    }

    // Process all edge mesh points
    forAll (em.points(), pointI)
    {
        // Skip point if string index is already assigned
        if (targetEdgeStrings[pointI] >= 0)
            continue;

        // Skip corner points
        if (em.pointEdges()[pointI].size() != 2)
            continue;

        label neighI1(UNDEF_LABEL);
        label neighI2(UNDEF_LABEL);
        findNeighborEdgeMeshPoints(em, pointI, neighI1, neighI2);
        stringifyEdgeMeshPoints(em, targetEdgeStrings, pointI, neighI1, neighI2, nStrings);
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

// Help function to project point coordinates to closest of the two
// edges connected to the closest point index in the edge mesh

point projectPointToClosestEdge
(
    const point pt,
    const edgeMesh& em,
    const label closestEdgePointI
)
{
    const point closestPoint = em.points()[closestEdgePointI];
    const point c2pt = pt - closestPoint;

    scalar minDistance = GREAT;
    vector deltaVec(ZERO_VECTOR);
    scalar finalDotProd(Zero);

    forAll(em.pointEdges()[closestEdgePointI], pointEdgeI)
    {
        const label edgeI = em.pointEdges()[closestEdgePointI][pointEdgeI];
        label endPointI = em.edges()[edgeI][0];
        if (endPointI == closestEdgePointI)
        {
            endPointI = em.edges()[edgeI][1];
        }

        const point endPoint = em.points()[endPointI];
        const point edgeVec = (endPoint - closestPoint) / mag(endPoint - closestPoint);
        const double dotProd = c2pt & edgeVec;
        const point projPoint = dotProd * edgeVec;
        const double distance = mag(projPoint - c2pt);

        if ((dotProd > 0) and (distance < minDistance))
        {
            minDistance = distance;
            deltaVec = projPoint;
            finalDotProd = dotProd;
        }

        if (finalDotProd > 1.0)
        {
            Info << "Warning: Detected extrapolation in feature edge projection of point " << pt << " by factor " << finalDotProd << endl;
        }

        if (finalDotProd > 10.0)
        {
            FatalError << "Error: Detected too large extrapolation in feature edge projection of point " << pt << " by factor " << finalDotProd << ". One reason might be that the target feature edges are too far away from initial feature edges, so that the mapping from intial to target (by edge point proximity) fails. Second reason might be that the smoothing step length is too large, and some feature edges get jumped over in smoothing." << endl << abort(FatalError);
        }
    }

    // Return the projection with a smallest positive vector scaling
    // factor (the point projection hits the edge), or the closest
    // point if neither edge results in a good projection
    if (minDistance < GREAT)
    {
        return closestPoint + deltaVec;
    }

    return closestPoint;
}

// Calculate new feature edge point positions and synchronize.
// Projects neighboring surface mesh points to closest
// feature edges and use the median as a new feature edge point.

int calculateFeatureEdgeProjections
(
    const fvMesh& mesh,
    const edgeMesh& em,
    const labelList& closestEdgePointIs,
    const boolList& isInternalPoint,
    const boolList& isFeatureEdgePoint,
    const boolList& isCornerPoint,
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
        const labelList neighIs = findNeighborSurfacePoints(mesh, pointI, isInternalPoint, isFeatureEdgePoint, isCornerPoint);

        // Project and accumulate data
        forAll (neighIs, neighI)
        {
            const label closestEdgePointI = closestEdgePointIs[pointI];
            const point neighPoint = mesh.points()[neighIs[neighI]];
            const point projPoint = projectPointToClosestEdge(neighPoint, em, closestEdgePointI);
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

// Help function to check if given edge mesh point indices are
// neighbors

bool isEdgePointNeighbor
(
    const edgeMesh& em,
    const label pointI1,
    const label pointI2
)
{
    forAll(em.pointEdges()[pointI1], pointEdgeI)
    {
        const label edgeI = em.pointEdges()[pointI1][pointEdgeI];
        const label endPointI1 = em.edges()[edgeI][0];
        const label endPointI2 = em.edges()[edgeI][1];

        if ((endPointI1 == pointI1) and (endPointI2 == pointI2))
        {
            return true;
        }

        if ((endPointI1 == pointI2) and (endPointI2 == pointI1))
        {
            return true;
        }
    }

    return false;
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

    FatalError
        << "Did not find surface intersection for pointI " << pointI
        << " at " << origPoint << endl
        << abort(FatalError);

    return UNDEF_VECTOR;
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
    labelList& closestEdgePointIs,
    const triSurface& surf,
    const indexedOctree<treeDataTriSurface>& tree,
    const double meshMaxEdgeLength,
    const labelList& targetEdgeStrings,
    const labelList& pointStrings
)
{
    // Calculate new feature edge point positions a priori
    vectorList featureEdgeProjections(mesh.nPoints(), ZERO_VECTOR);
    labelList nFeatureEdgeProjections(mesh.nPoints(), 0);
    calculateFeatureEdgeProjections(mesh, targetEdges, closestEdgePointIs, isInternalPoint, isFeatureEdgePoint, isCornerPoint, featureEdgeProjections, nFeatureEdgeProjections);

    forAll(mesh.points(), pointI)
    {
        if (isInternalPoint[pointI])
            continue;

        const vector origPoint = newPoints[pointI];

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

            // Update closest edge mesh point index if needed
            const label closestEdgePointI = closestEdgePointIs[pointI];
            label newClosestEdgePointI;
            label newStringI = UNDEF_LABEL;
            findClosestEdgeMeshPointIndex(newPoint, targetEdges, false, true, true, pointStrings[pointI], targetEdgeStrings, newClosestEdgePointI, newStringI);

            if (newClosestEdgePointI != closestEdgePointI)
            {
                closestEdgePointIs[pointI] = newClosestEdgePointI;
                // Info << "newClosestEdgePoint for " << pointI << " is " << newClosestEdgePointI << endl;
            }

            if (pointStrings[pointI] != newStringI)
            {
                FatalError << "Internal sanity check failed: String number changed for pointI " << pointI << " from " << pointStrings[pointI] << " to " << newStringI << endl << abort(FatalError);
            }

            continue;
        }

        // Project to closest tri face in normal or opposite direction
        // -----------------------------------------------------------
        if (isSmoothingSurfacePoint[pointI])
        {
            const point pointNormal = pointNormals[pointI];
            const double searchDistance = 10 * meshMaxEdgeLength;
            newPoints[pointI] = findIntersection(tree, pointI, origPoint, pointNormal, searchDistance);
        }
    }

    return 0;
}

// ************************************************************************* //
