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
    const triSurfaceSearch& searchSurfaces,
    const double boundaryMaxPointBlendingFraction
)
{
    forAll(mesh.points(), pointI)
    {
        // TBA
        // const vector newCoords = cCoords - boundaryMaxPointBlendingFraction * pVec;
        // Update point coordinates
        //newPoints[pointI] = newCoords;
        continue;
    }

    return 0;
}

// ************************************************************************* //
