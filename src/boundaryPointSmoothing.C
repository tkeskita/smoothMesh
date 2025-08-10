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
// to corner points only.

label findClosestEdgeMeshPointIndex
(
    const point pt,
    const edgeMesh& em,
    const bool mustBeCorner,
    const bool mustNotBeCorner
)
{
    // TODO: Replace linear search with an octree search for efficiency

    scalar distance = GREAT;
    label closestI = UNDEF_LABEL;

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

        const point p = em.points()[pointI];
        const scalar testDistance = mag(pt - p);
        if (testDistance < distance)
        {
            distance = testDistance;
            closestI = pointI;
        }
    }

    return closestI;
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
    boolList& isFeatureEdgePoint,
    boolList& isCornerPoint,
    vectorList& cornerPoints,
    boolList& isLayerSurfacePoint,
    boolList& isSmoothingSurfacePoint,
    boolList& isFrozenSurfacePoint
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

                if ((initEdges.points().size() > 0) and (targetEdges.points().size() > 0))
                {
                    const point pt = mesh.points()[pointI];
                    const label closestInitPointI = findClosestEdgeMeshPointIndex(pt, initEdges, false, false);
                    const point closestInitPoint = initEdges.points()[closestInitPointI];

                    // Corner points
                    if ((initEdges.pointEdges()[closestInitPointI].size() != 2)
                        and (mag(pt - closestInitPoint) < DISTANCE_TOLERANCE))
                    {
                        isCornerPoint[pointI] = true;
                        const label closestCornerPointI = findClosestEdgeMeshPointIndex(closestInitPoint, targetEdges, true, false);
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
                if (findIndex(layerPatchIds, patchI) >= 0)
                {
                    isLayerSurfacePoint[pointI] = true;
                    nLayerSurfacePoints++;
                }

                // Smoothing surface points
                if (findIndex(smoothingPatchIds, patchI) >= 0)
                {
                    isSmoothingSurfacePoint[pointI] = true;
                    nSmoothingSurfacePoints++;
                    continue;
                }
                // Frozen surface points
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

    Info << "Detected number of corner points: " << nSumCornerPoints << endl;
    Info << "Detected number of feature edge points: " << nSumFeatureEdgePoints << endl;
    Info << "Detected number of layer surface points: " << nSumLayerSurfacePoints << endl;
    Info << "Detected number of smoothing surface points: " << nSumSmoothingSurfacePoints << endl;
    Info << "Detected number of frozen surface points: " << nSumFrozenSurfacePoints << endl;
    Info << endl;

    return 0;
}


// Project point to closest edge in edge mesh

point projectPointToClosestEdge
(
    const point pt,
    const edgeMesh& em,
    const label meshPointI
)
{
    const label pointI = findClosestEdgeMeshPointIndex(pt, em, false, false);
    const point closestPoint = em.points()[pointI];
    const point c2pt = pt - closestPoint;

    scalar minDistance = GREAT;
    vector deltaVec(ZERO_VECTOR);

    forAll(em.pointEdges()[pointI], i)
    {
        const label edgeI = em.pointEdges()[pointI][i];
        label endPointI = em.edges()[edgeI][0];
        if (endPointI == pointI)
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
    const boolList& isFeatureEdgePoint,
    const boolList& isCornerPoint,
    const vectorList& cornerPoints,
    const edgeMesh& targetEdges,
    const triSurface& surf,
    const indexedOctree<treeDataTriSurface>& tree,
    const double meshMaxEdgeLength
)
{
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
            newPoints[pointI] = projectPointToClosestEdge(origPoint, targetEdges, pointI);
            continue;
        }

        // Project to closest tri face in normal or opposite direction
        // -----------------------------------------------------------
        const point pointNormal = pointNormals[pointI];
        const double searchDistance = 10 * meshMaxEdgeLength;

        newPoints[pointI] = findIntersection(tree, pointI, origPoint, pointNormal, searchDistance);
    }

    return 0;
}

// ************************************************************************* //
