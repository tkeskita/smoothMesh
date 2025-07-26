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

// Projection of boundary points to feature edges and boundary surfaces

int projectBoundaryPointsToEdgesAndSurfaces
(
    const fvMesh& mesh,
    pointField& newPoints,
    const boolList& isInternalPoint,
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
        const label faceI = hitInfo.index();

        // Search closest mesh point index
        const label vertexPointI = findVertexPointI(surf, faceI, closestPoint);
        const point vertexPoint = surf.localPoints()[vertexPointI];

        // Get faces surrounding the vertex
        const labelList closestFaces = surf.pointFaces()[vertexPointI];

        // Project original point to each face normal plane and
        // calculate edge angles at the surface point
        vectorList projPoints(closestFaces.size(), ZERO_VECTOR);
        scalarList edgeAngles(closestFaces.size(), ZERO);
        projectPointToFaces(surf, origPoint, vertexPointI, closestFaces, projPoints, edgeAngles);

        // Primary data includes only points which are closer to face
        // center compared to vertex point
        point primaryPoint(ZERO_VECTOR);
        scalar sumPrimaryDist(ZERO);

        // Secondary (fall back) data includes all projected points
        point secondaryPoint(ZERO_VECTOR);
        scalar sumSecondaryDist(ZERO);

        point newPoint(ZERO_VECTOR);
        bool isProjectedPointOnFace(false);
        forAll (closestFaces, i)
        {
            const point faceCenter = surf.faceCentres()[closestFaces[i]];
            const point projPoint = projPoints[i];

            const scalar vertexPointToFaceCenter = mag(vertexPoint - faceCenter);
            const scalar projPointToFaceCenter = mag(projPoint - faceCenter);
            const scalar origPointToPlaneDist = mag(origPoint - projPoint);
            // const scalar invDist = 1.0 / pow4(origPointToPlaneDist);
            const scalar invDist = 1.0 / sqr(origPointToPlaneDist);
            // const scalar invDist = 1.0 / origPointToPlaneDist;
            // const scalar invDist = 1.0;

            const scalar angleFactor = edgeAngles[i] / sum(edgeAngles);

            // If point is on a plane, projection is already complete
            if (origPointToPlaneDist < VSMALL)
            {
                isProjectedPointOnFace = true;
                primaryPoint = projPoint;
                break;
            }

            // Add projection to secondary data
            secondaryPoint += projPoint * invDist * angleFactor;
            sumSecondaryDist += invDist * angleFactor;

            // Add projection to primary data, but only if the
            // projected point is closer to face centers than the
            // vertex point.
            if (projPointToFaceCenter < vertexPointToFaceCenter)
            {
                primaryPoint += projPoint * invDist * angleFactor;
                sumPrimaryDist += invDist * angleFactor;
            }
        }

        // New point coordinate
        if (isProjectedPointOnFace)
        {
            newPoint = primaryPoint;
        }
        else if (sumPrimaryDist > VSMALL)
        {
            newPoint = primaryPoint / sumPrimaryDist;
        }
        else
        {
            newPoint = secondaryPoint / sumSecondaryDist;
        }

        // Project the new point on the closest face normal plane
        // along the vector from original point to new point.
        const point fc = surf.faceCentres()[faceI];
        const vector fNormal =
            surf.faceNormals()[faceI] / mag(surf.faceNormals()[faceI]);
        const vector o2fc = fc - origPoint;
        const double fcProj = o2fc & fNormal;
        const vector o2n = newPoint - origPoint;
        const double nProj = o2n & fNormal;
        const vector projCoords = origPoint + fcProj/nProj * o2n;
        newPoint = projCoords;

        // Blend and save
        const vector newCoords =
            (1 - boundaryMaxPointBlendingFraction) * origPoint
            + boundaryMaxPointBlendingFraction * newPoint;
        newPoints[pointI] = newCoords;
    }

    return 0;
}

// ************************************************************************* //
