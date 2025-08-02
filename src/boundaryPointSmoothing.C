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

// Find corner points and target locations for the corners using the
// provided edge meshes

int findCornersAndFeatureEdges
(
    const fvMesh& mesh,
    const edgeMesh& initEdges,
    const edgeMesh& targetEdges,
    const boolList& isInternalPoint,
    boolList& isFeatureEdgePoint,
    boolList& isCornerPoint,
    vectorList& cornerPoints
)
{
    label nCornerPoints = 0;
    label nFeatureEdgePoints = 0;

    forAll(mesh.points(), pointI)
    {
        if (isInternalPoint[pointI])
            continue;

        const point pt = mesh.points()[pointI];
        const label closestInitPointI = findClosestEdgeMeshPointIndex(pt, initEdges, false, false);
        const point closestInitPoint = initEdges.points()[closestInitPointI];

        // Identify and handle corner points
        if (initEdges.pointEdges()[closestInitPointI].size() != 2)
        {
            isCornerPoint[pointI] = true;
            const label closestCornerPointI = findClosestEdgeMeshPointIndex(closestInitPoint, targetEdges, true, false);
            cornerPoints[pointI] = targetEdges.points()[closestCornerPointI];
            nCornerPoints++;
            continue;
        }

        // Identify feature edges
        if (mag(pt - closestInitPoint) < DISTANCE_TOLERANCE)
        {
            isFeatureEdgePoint[pointI] = true;
            nFeatureEdgePoints++;
        }
    }

    const label nSumCornerPoints  = returnReduce(nCornerPoints, sumOp<label>());
    const label nSumFeatureEdgePoints  = returnReduce(nFeatureEdgePoints, sumOp<label>());
    Info << "Detected number of corner points: " << nSumCornerPoints << endl;
    Info << "Detected number of feature edge points: " << nSumFeatureEdgePoints << endl;
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

// Projection of boundary points to feature edges and boundary surfaces

int projectBoundaryPointsToEdgesAndSurfaces
(
    const fvMesh& mesh,
    pointField& newPoints,
    const boolList& isInternalPoint,
    const boolList& isFeatureEdgePoint,
    const boolList& isCornerPoint,
    const vectorList& cornerPoints,
    const edgeMesh& targetEdges,
    const triSurface& surf,
    const indexedOctree<treeDataTriSurface>& tree,
    const double meshMaxEdgeLength,
    const double boundaryMaxPointBlendingFraction
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

        // Project to closest tri face
        // ---------------------------
        // Find closest point and closest face info
        const pointIndexHit hitInfo = tree.findNearest(origPoint, 8.0*sqr(meshMaxEdgeLength));
        if (! hitInfo.hit())
        {
            FatalError
                << "findNearest for surface failed for pointI " << pointI
                << " located at " << origPoint << ". Search distance was "
                << meshMaxEdgeLength << "." << endl
                << exit(FatalError);
        }

        const point& closestPoint = hitInfo.hitPoint();
        const vector newCoords =
            (1 - boundaryMaxPointBlendingFraction) * origPoint
            + boundaryMaxPointBlendingFraction * closestPoint;
        newPoints[pointI] = newCoords;
    }

    return 0;
}

// ************************************************************************* //
