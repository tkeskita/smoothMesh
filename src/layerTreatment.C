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


// Help function to test whether a point belongs to a given island

bool isPointInIsland
(
    const label pointI,
    const label islandI,
    const labelList& prismIslands1,
    const labelList& prismIslands2,
    const labelList& prismIslands3
)
{
    return ((prismIslands1[pointI] == islandI) or
        (prismIslands2[pointI] == islandI) or
        (prismIslands3[pointI] == islandI)
    );
}


// Help function to set a given island info for a point index. If no
// island slots are available, reset all slots to ignored value and
// return zero to indicate that point has been reset.

label addIslandInfoForPoint
(
    const label pointI,
    const label pointNormalSourceI,
    const label islandI,
    const vectorList& pointNormals,
    labelList& prismIslands1,
    labelList& prismIslands2,
    labelList& prismIslands3,
    labelList& pointHops1,
    labelList& pointHops2,
    labelList& pointHops3,
    labelList& pointNormalSource1,
    labelList& pointNormalSource2,
    labelList& pointNormalSource3,
    vectorList& pointNormals1,
    vectorList& pointNormals2,
    vectorList& pointNormals3
)
{
    if (isPointInIsland(pointI, islandI, prismIslands1, prismIslands2, prismIslands3))
    {
        FatalError << "Point " << pointI << " is already in island " << islandI << endl << abort(FatalError);
    }

    if (prismIslands1[pointI] == UNDEF_LABEL)
    {
        prismIslands1[pointI] = islandI;
        pointNormalSource1[pointI] = pointNormalSourceI;
        pointNormals1[pointI] = pointNormals[pointNormalSourceI];
        pointHops1[pointI] = 0;
        return 1;
    }
    else if (prismIslands2[pointI] == UNDEF_LABEL)
    {
        prismIslands2[pointI] = islandI;
        pointNormalSource2[pointI] = pointNormalSourceI;
        pointNormals2[pointI] = pointNormals[pointNormalSourceI];
        pointHops2[pointI] = 0;
        return 2;
    }
    else if (prismIslands3[pointI] == UNDEF_LABEL)
    {
        prismIslands3[pointI] = islandI;
        pointNormalSource3[pointI] = pointNormalSourceI;
        pointNormals3[pointI] = pointNormals[pointNormalSourceI];
        pointHops3[pointI] = 0;
        return 3;
    }

    // All slots are already in use, reset slots to ignored
    prismIslands1[pointI] = IGNORED_LABEL;
    prismIslands2[pointI] = IGNORED_LABEL;
    prismIslands3[pointI] = IGNORED_LABEL;
    pointNormalSource1[pointI] = UNDEF_LABEL;
    pointNormalSource2[pointI] = UNDEF_LABEL;
    pointNormalSource3[pointI] = UNDEF_LABEL;
    pointNormals1[pointI] = Zero;
    pointNormals2[pointI] = Zero;
    pointNormals3[pointI] = Zero;
    pointHops1[pointI] = UNDEF_LABEL;
    pointHops2[pointI] = UNDEF_LABEL;
    pointHops3[pointI] = UNDEF_LABEL;
    return 0;
}


// Help function which propagates the island id of given point index
// to neighboring prismatic points on the boundary

label propagateIslandInfoOnBoundary
(
    const fvMesh& mesh,
    const label startPointI,
    boolList& isVisitedPoint,
    const boolList& isPrismaticPoint,
    const boolList& isLayerSurfacePoint,
    const vectorList& pointNormals,
    labelList& prismIslands1,
    labelList& prismIslands2,
    labelList& prismIslands3,
    labelList& pointHops1,
    labelList& pointHops2,
    labelList& pointHops3,
    labelList& pointNormalSource1,
    labelList& pointNormalSource2,
    labelList& pointNormalSource3,
    vectorList& pointNormals1,
    vectorList& pointNormals2,
    vectorList& pointNormals3
)
{
    label n = 1;
    label nTot = 1;
    const label islandI = prismIslands1[startPointI];
    if (islandI < 0)
    {
        FatalError << "Illegal propagation islandI " << islandI << endl << abort(FatalError);
    }

    // Propagate island id among prismatic points in a loop until
    // propagation stops
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

                if (prismIslands1[neighI] == islandI)
                {
                    isVisitedPoint[pointI] = true;
                    addIslandInfoForPoint(pointI, pointI, islandI, pointNormals, prismIslands1, prismIslands2, prismIslands3, pointHops1, pointHops2, pointHops3, pointNormalSource1, pointNormalSource2, pointNormalSource3, pointNormals1, pointNormals2, pointNormals3);
                    ++n;
                    ++nTot;
                    // Info << mesh.points()[pointI] << endl;
                }
            }
        }
    }

    return nTot;
}


// Help function to add non-prismatic edge points to the island

label addEdgePointsToIsland
(
    const fvMesh& mesh,
    const label islandI,
    const boolList& isPrismaticPoint,
    const boolList& isLayerSurfacePoint,
    const vectorList& pointNormals,
    labelList& prismIslands1,
    labelList& prismIslands2,
    labelList& prismIslands3,
    labelList& pointHops1,
    labelList& pointHops2,
    labelList& pointHops3,
    labelList& pointNormalSource1,
    labelList& pointNormalSource2,
    labelList& pointNormalSource3,
    vectorList& pointNormals1,
    vectorList& pointNormals2,
    vectorList& pointNormals3
)
{
    labelList edgePoints;

    // Identify edge points
    forAll(mesh.boundary(), patchI)
    {
        const label startI = mesh.boundary()[patchI].start();
        const label endI = startI + mesh.boundary()[patchI].Cf().size();

        for (label faceI = startI; faceI < endI; faceI++)
        {
            const face& f = mesh.faces()[faceI];

            // Is this face part of island?
            bool isIslandFace = false;
            forAll (f, facePointI)
            {
                const label pointI = mesh.faces()[faceI][facePointI];
                if (isPointInIsland(pointI, islandI, prismIslands1, prismIslands2, prismIslands3))
                {
                    isIslandFace = true;
                    break;
                }
            }

            // Add non-island face points to edge points
            if (isIslandFace)
            {
                forAll (f, facePointI)
                {
                    const label pointI = mesh.faces()[faceI][facePointI];
                    if (! isPointInIsland(pointI, islandI, prismIslands1, prismIslands2, prismIslands3))
                    {
                        if (findIndex(edgePoints, pointI) == -1)
                        {
                            edgePoints.append(pointI);
                        }
                    }
                }
            }
        }
    }

    // Propagate point normal source to edge points
    label n = 1;
    label nTot = 0;

    while (n > 0)
    {
        n = 0;

        forAll(edgePoints, edgePointI)
        {
            const label pointI = edgePoints[edgePointI];

            if (! isLayerSurfacePoint[pointI])
                continue;
            if (isPointInIsland(pointI, islandI, prismIslands1, prismIslands2, prismIslands3))
                continue;

            forAll (mesh.pointPoints()[pointI], pointPointI)
            {
                const label neighI = mesh.pointPoints()[pointI][pointPointI];
                if (! isLayerSurfacePoint[neighI])
                    continue;

                if (isPointInIsland(neighI, islandI, prismIslands1, prismIslands2, prismIslands3))
                {
                    addIslandInfoForPoint(pointI, pointNormalSource1[neighI], islandI, pointNormals, prismIslands1, prismIslands2, prismIslands3, pointHops1, pointHops2, pointHops3, pointNormalSource1, pointNormalSource2, pointNormalSource3, pointNormals1, pointNormals2, pointNormals3);
                    ++n;
                    ++nTot;
                    // Pout << "  Added edge point " << pointI << " to island " << islandI << " res " << res << endl;
                    break;
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
    const vectorList& pointNormals,
    labelList& prismIslands1,
    labelList& prismIslands2,
    labelList& prismIslands3,
    labelList& pointHops1,
    labelList& pointHops2,
    labelList& pointHops3,
    labelList& pointNormalSource1,
    labelList& pointNormalSource2,
    labelList& pointNormalSource3,
    vectorList& pointNormals1,
    vectorList& pointNormals2,
    vectorList& pointNormals3
)
{
    // Processor number (MPI rank)
    const label myProcNo = Pstream::myProcNo();

    // Next available island ID
    label islandI = myProcNo * maxIds;

    // Storage of processed points
    boolList isVisitedPoint(mesh.nPoints(), false);

    // Main identification loop
    while (true)
    {
        // Find an unprocessed prism boundary point
        const label startPointI = findNewPrismBoundaryPointI(mesh, isVisitedPoint, isPrismaticPoint, isLayerSurfacePoint);

        // No more unprocessed prism points, stop
        if (startPointI == UNDEF_LABEL)
        {
            break;
        }

        // Add and propagate the island id to all unprocessed prism
        // boundary points points on this island
        addIslandInfoForPoint(startPointI, startPointI, islandI, pointNormals, prismIslands1, prismIslands2, prismIslands3, pointHops1, pointHops2, pointHops3, pointNormalSource1, pointNormalSource2, pointNormalSource3, pointNormals1, pointNormals2, pointNormals3);
        // Info << "Starting pointI " << pointI << " at " << mesh.points()[pointI] << endl;

        const label n = propagateIslandInfoOnBoundary(mesh, startPointI, isVisitedPoint, isPrismaticPoint, isLayerSurfacePoint, pointNormals, prismIslands1, prismIslands2, prismIslands3, pointHops1, pointHops2, pointHops3, pointNormalSource1, pointNormalSource2, pointNormalSource3, pointNormals1, pointNormals2, pointNormals3);

        // Add the non-prismatic edge points to this island
        const label ne = addEdgePointsToIsland(mesh, islandI, isPrismaticPoint, isLayerSurfacePoint, pointNormals, prismIslands1, prismIslands2, prismIslands3, pointHops1, pointHops2, pointHops3, pointNormalSource1, pointNormalSource2, pointNormalSource3, pointNormals1, pointNormals2, pointNormals3);

        Pout << "Island " << islandI << " has " << n << " prism points" << " and " << ne << " edge points" << endl;

        if (n < 1)
        {
            break;
        }

        // Reserve next id
        ++islandI;

        if (islandI >= (myProcNo + 1) * maxIds)
        {
            FatalError << "Exceeded maximum number of islands " << maxIds << endl << abort(FatalError);
        }
    }

    // Synchronize island IDs among processors



    return 0;
}
