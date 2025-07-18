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

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Find index of point in triSurface mesh which is part of a given
// face index faceI and is located nearest to coordinates coords

label findClosestPointI
(
    const triSurface& surf,
    const label faceI,
    const point& coords
)
{
    scalar distance = VGREAT;
    label closestPointI = -1;

    forAll(surf.localFaces()[faceI], pointI)
    {
        const label localPointI = surf.localFaces()[faceI][pointI];
        const point& v = surf.localPoints()[localPointI];
        const scalar testDistance = mag(coords - v);
        if (testDistance < distance)
        {
            closestPointI = localPointI;
            distance = testDistance;
        }
    }

    return closestPointI;
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

        // const point& pt = mesh.points()[pointI];
        const vector cCoords = newPoints[pointI];
        pointIndexHit hitInfo = tree.findNearest(cCoords, sqr(meshMaxEdgeLength));

        if (! hitInfo.hit())
        {
            FatalErrorInFunction
                << "No nearest points for pointI " << pointI << endl
                << " within distance " << meshMaxEdgeLength
                << exit(FatalError);
        }

        // Closest point and face info
        const point& closestPoint = hitInfo.hitPoint();
        const label faceI = hitInfo.index();
        const point fCoords = surf.faceCentres()[faceI];
        const vector faceNormal =
            surf.faceNormals()[faceI] / mag(surf.faceNormals()[faceI]);

        // Simple and naive: Project point to the normal plane of the
        // closest face
        const vector cf = cCoords - fCoords;
        const double dotProd = cf & faceNormal;
        const vector pCoords = cCoords - dotProd * faceNormal;
        // Info << "  - projected point coords " << pCoords << endl;
        // Info << "    closest point coords   " << closestPoint << endl;

        // Search closest point index
        const label closestPointI = findClosestPointI(surf, faceI, closestPoint);
        Info << "Closest pointI is " << closestPointI << " at " << surf.localPoints()[closestPointI] << endl;

        // Info << "Hit point for pointI " << pointI << " at " << pt << " is " << closestPoint << " with face index " << hitInfo.index() << " and face center " << faceCenter << endl;

        // Blend and save
        const vector newCoords =
            (1 - boundaryMaxPointBlendingFraction) * cCoords
            + boundaryMaxPointBlendingFraction * pCoords;
        newPoints[pointI] = newCoords;
    }

    return 0;
}

// ************************************************************************* //
