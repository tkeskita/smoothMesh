/*---------------------------------------------------------------------------*\
Library
    Layer Treatment

Description
    Special treatment of prismatic boundary layers, with aim to
    control the length of prismatic side edges on boundary cell layers
    (cell thickness).
\*---------------------------------------------------------------------------*/

#include "fvMesh.H"
#include "smoothMeshCommon.H"

using namespace Foam;

// Maximum number of prism islands allowed in one processor
const label maxIds = 10000;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Help function to find an unprocessed prismatic boundary point

label findNewPrismBoundaryPointI
(
    const fvMesh& mesh,
    boolList& isVisitedPoint,
    const boolList& isPrismaticPoint,
    const boolList& isLayerSurfacePoint
)
{
    forAll(mesh.points(), pointI)
    {
        if (! isPrismaticPoint[pointI])
            continue;
        if (! isLayerSurfacePoint[pointI])
            continue;
        if (isVisitedPoint[pointI])
            continue;

        // Mark this point as processed
        isVisitedPoint[pointI] = true;

        // Check that all neighbor boundary points are prismatic
        // before accepting
        forAll (mesh.pointPoints()[pointI], pointPointI)
        {
            const label neighI = mesh.pointPoints()[pointI][pointPointI];
            if ((isLayerSurfacePoint[neighI]) and (! isPrismaticPoint[neighI]))
                continue;
        }

        // Accept point
        return pointI;
    }

    // All points processed
    return UNDEF_LABEL;
}


// Help function to set a given island id for a point index. If no
// island slots are available, reset all slots to ignored value and
// return false to indicate reset.

bool addIslandIdForPoint
(
    const label pointI,
    const label id,
    labelList& prismIslands1,
    labelList& prismIslands2,
    labelList& prismIslands3
)
{
    if (prismIslands1[pointI] == UNDEF_LABEL)
    {
        prismIslands1[pointI] = id;
        return true;
    }
    else if (prismIslands2[pointI] == UNDEF_LABEL)
    {
        prismIslands2[pointI] = id;
        return true;
    }
    else if (prismIslands3[pointI] == UNDEF_LABEL)
    {
        prismIslands3[pointI] = id;
        return true;
    }

    // All slots are already in use, reset slots to ignored
    prismIslands1[pointI] = IGNORED_LABEL;
    prismIslands2[pointI] = IGNORED_LABEL;
    prismIslands3[pointI] = IGNORED_LABEL;
    return false;
}


// Help function which propagates the island id of given point index
// to neighboring prismatic points and to island edge points on the
// boundary

label propagatePrismIndex
(
    const fvMesh& mesh,
    const label pointI,
    boolList& isVisitedPoint,
    const boolList& isPrismaticPoint,
    const boolList& isLayerSurfacePoint,
    labelList& prismIslands1,
    labelList& prismIslands2,
    labelList& prismIslands3
)
{
    label n = 1;
    label nTot = 1;
    const label id = prismIslands1[pointI];
    if (id < 0)
    {
        FatalError << "Illegal propagation id " << id << endl << abort(FatalError);
    }

    // Propagate island id in a loop until propagation stops
    while (n > 0)
    {
        n = 0;

        forAll(mesh.points(), pointI)
        {
            if (! isPrismaticPoint[pointI])
                continue;
            if (! isLayerSurfacePoint[pointI])
                continue;
            if (prismIslands1[pointI] != UNDEF_LABEL)
                continue;

            forAll (mesh.pointPoints()[pointI], pointPointI)
            {
                const label neighI = mesh.pointPoints()[pointI][pointPointI];
                if (! isLayerSurfacePoint[neighI])
                    continue;
                if (isVisitedPoint[pointI])
                    continue;

                if (prismIslands1[neighI] == id)
                {
                    prismIslands1[pointI] = id;
                    isVisitedPoint[pointI] = true;
                    ++n;
                    ++nTot;
                    Info << mesh.points()[pointI] << endl;
                }
            }
        }
    }

    return nTot;
}


// Function to identify and number the prismatic islands on the
// boundaries

int identifyPrismaticBoundaryIslands
(
    const fvMesh& mesh,
    const boolList& isPrismaticPoint,
    const boolList& isLayerSurfacePoint,
    labelList& prismIslands1,
    labelList& prismIslands2,
    labelList& prismIslands3
)
{
    // Processor number (MPI rank)
    const label myProcNo = Pstream::myProcNo();

    // Next available island ID
    label id = myProcNo * maxIds;

    // Storage of processed points
    boolList isVisitedPoint(mesh.nPoints(), false);

    labelList dummy;

    // Main identification loop
    while (true)
    {
        // Find an unprocessed prism boundary point
        const label pointI = findNewPrismBoundaryPointI(mesh, isVisitedPoint, isPrismaticPoint, isLayerSurfacePoint);

        // No more unprocessed prism points, stop
        if (pointI == UNDEF_LABEL)
        {
            break;
        }

        // Add and propagate the island id to all unprocessed prism
        // boundary points and island edge points on this island
        addIslandIdForPoint(pointI, id, prismIslands1, dummy, dummy);
        Info << "Starting pointI " << pointI << " at " << mesh.points()[pointI] << endl;

        const label n = propagatePrismIndex(mesh, pointI, isVisitedPoint, isPrismaticPoint, isLayerSurfacePoint, prismIslands1, prismIslands2, prismIslands3);

        FatalError << "WIP stop. Island " << id << " has " << n << " prism points" << endl << abort(FatalError);
        
        if (n < 1)
        {
            break;
        }

        // Reserve next id
        ++id;

        if (id >= (myProcNo + 1) * maxIds)
        {
            FatalError << "Exceeded maximum number of islands " << maxIds << endl << abort(FatalError);
        }
    }

    // WIP TBA: Synchronize island IDs among processors

    return 0;
}
