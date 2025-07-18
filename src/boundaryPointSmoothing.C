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

        // Closest face info
        // const point& closestPoint = hitInfo.hitPoint();
        const label faceI = hitInfo.index();
        const point fCoords = surf.faceCentres()[faceI];
        const vector faceNormal =
            surf.faceNormals()[faceI] / mag(surf.faceNormals()[faceI]);

        // Simple and naive: Project point to the normal plane of the
        // closest face
        const vector cf = cCoords - fCoords;
        const double dotProd = cf & faceNormal;
        const vector pCoords = cCoords - dotProd * faceNormal;

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
