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

// TODO: Find out how to get rid of the need to copy this function here

double edgeEdgeAngleLocalCopy
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


// Find index of point in triSurface mesh which is part of a given
// face index faceI and is located nearest to coordinates coords

label findVertexPointI
(
    const triSurface& surf,
    const label faceI,
    const point& coords
)
{
    scalar distance = VGREAT;
    label vertexPointI = -1;

    forAll(surf.localFaces()[faceI], pointI)
    {
        const label localPointI = surf.localFaces()[faceI][pointI];
        const point& v = surf.localPoints()[localPointI];
        const scalar testDistance = mag(coords - v);
        if (testDistance < distance)
        {
            vertexPointI = localPointI;
            distance = testDistance;
        }
    }

    return vertexPointI;
}


// Calculate the angle between edges of a tri face meeting at a point

double calcTriSurfaceEdgeAngle
(
    const triSurface& surf,
    const label faceI,
    const label vertexI
)
{
    vector endPoint1(UNDEF_VECTOR);
    vector endPoint2(UNDEF_VECTOR);

    // Deduce end points
    forAll(surf.localFaces()[faceI], i)
    {
        const label pointI = surf.localFaces()[faceI][i];
        if (pointI==vertexI)
        {
            continue;
        }
        else if (endPoint1 == UNDEF_VECTOR)
        {
            endPoint1 = surf.localPoints()[pointI];
        }
        else
        {
            endPoint2 = surf.localPoints()[pointI];
        }
    }

    const point startPoint = surf.localPoints()[vertexI];

    return edgeEdgeAngleLocalCopy(startPoint, endPoint1, endPoint2);
}


// Project the given point coordinate to the normal planes of given
// faces of the triSurface and save projected coordinate to projPoints.
// Additionally calculated angle of each edges of each face meeting at
// vertexPointI.

int projectPointToFaces
(
    const triSurface& surf,
    const point& coords,
    const label vertexPointI,
    const labelList& faceIs,
    vectorList& projPoints,
    scalarList& edgeAngles
)
{
    forAll (faceIs, i)
    {
        // Project point
        const label faceI = faceIs[i];
        const point fCoords = surf.faceCentres()[faceI];
        const vector fNormal =
            surf.faceNormals()[faceI] / mag(surf.faceNormals()[faceI]);
        const vector cf = coords - fCoords;
        const double dotProd = cf & fNormal;
        const vector resCoords = coords - dotProd * fNormal;
        projPoints[i] = resCoords;
        edgeAngles[i] = calcTriSurfaceEdgeAngle(surf, faceI, vertexPointI);
    }

    return 0;
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
            continue;
        }

        // Identify feature edges
        if (mag(pt - closestInitPoint) < DISTANCE_TOLERANCE)
        {
            isFeatureEdgePoint[pointI] = true;
        }
    }

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

    // if (meshPointI == 11)
    //   Info << "  pointEdges " << em.pointEdges()[pointI].size() << endl;
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

        // if (meshPointI == 11)
        //    Info << "    pt " << pt << " endPoint " << endPoint << " edgeVec " << edgeVec << " dotProd " << dotProd << " projPoint " << projPoint << " distance " << distance << endl;

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
        // if (meshPointI == 11)
        //     Info << "    deltaVec " << deltaVec << " from closestpoint " << closestPoint << endl;
        return closestPoint + deltaVec;
    }

    // if (meshPointI == 11)
    //    Info << "        using closest point" << endl;
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
        // Info << "Projecting pointI " << pointI << " at " << origPoint << endl;

        // if (pointI == 11)
        //      Info << pointI << " at " << origPoint << " iscorner " << isCornerPoint[pointI] << " isFeature " << isFeatureEdgePoint[pointI] << endl;

        // Project to closest corner points
        // --------------------------------
        if (isCornerPoint[pointI])
        {
            newPoints[pointI] = cornerPoints[pointI];
            // Info << "corner point " << pointI << " to " << newPoints[pointI] << endl;
            continue;
        }

        // Project to closest feature edges
        // --------------------------------
        if (isFeatureEdgePoint[pointI])
        {
            newPoints[pointI] = projectPointToClosestEdge(origPoint, targetEdges, pointI);
            // Info << "edge point " << pointI << " to " << newPoints[pointI] << endl;
            // if (pointI == 11)
            //     Info << "   projected final point is " << newPoints[pointI] << endl;
            continue;
        }

        // Project to closest tri face
        // ---------------------------
        // Find closest point and closest face info
        const pointIndexHit hitInfo = tree.findNearest(origPoint, sqr(meshMaxEdgeLength));
        if (! hitInfo.hit())
        {
            FatalError
                << "No nearest points for pointI " << pointI << endl
                << " within distance " << meshMaxEdgeLength
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
