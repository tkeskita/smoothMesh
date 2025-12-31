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

// Maximum number of prism islands allowed in one processor. Since
// minimum number of processors is one, this number equals the total
// allowed maximum number of islands in the mesh.
const label maxIslands = 10000;

// Passive index start number
const label passiveIndexStart = 1000000000;

// Note: Maximum allowed number of processors is passiveIndexStart/maxIslands

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
        bool allNeighsArePrismatic = true;
        forAll (mesh.pointPoints()[pointI], pointPointI)
        {
            const label neighI = mesh.pointPoints()[pointI][pointPointI];
            if ((isLayerSurfacePoint[neighI]) and (! isPrismaticPoint[neighI]))
            {
                allNeighsArePrismatic = false;
            }
        }

        // Accept point
        if (allNeighsArePrismatic)
        {
            return pointI;
        }
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
    return
    (
        (prismIslands1[pointI] == islandI) or
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
    const label nLayer,
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
        pointHops1[pointI] = nLayer;
        return 1;
    }
    else if (prismIslands2[pointI] == UNDEF_LABEL)
    {
        prismIslands2[pointI] = islandI;
        pointNormalSource2[pointI] = pointNormalSourceI;
        pointNormals2[pointI] = pointNormals[pointNormalSourceI];
        pointHops2[pointI] = nLayer;
        return 2;
    }
    else if (prismIslands3[pointI] == UNDEF_LABEL)
    {
        prismIslands3[pointI] = islandI;
        pointNormalSource3[pointI] = pointNormalSourceI;
        pointNormals3[pointI] = pointNormals[pointNormalSourceI];
        pointHops3[pointI] = nLayer;
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
    vectorList& pointNormals3,
    const boolList& isProcessorPoint,
    label& nProcPrisms
)
{
    label n = 1;
    label nTot = 1;
    nProcPrisms = 0;
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
        Pout << "  Starting propagation round, nTot=" << nTot << endl;

        forAll(mesh.points(), pointI)
        {
            if (! isPrismaticPoint[pointI])
                continue;
            if (! isLayerSurfacePoint[pointI])
                continue;
            if (prismIslands1[pointI] != UNDEF_LABEL)
                continue;
            // Pout << "    Considering pointI " << pointI << endl;

            forAll (mesh.pointPoints()[pointI], pointPointI)
            {
                const label neighI = mesh.pointPoints()[pointI][pointPointI];
                if (! isLayerSurfacePoint[neighI])
                    continue;
                // Pout << "      Considering neighI " << neighI << endl;

                if (prismIslands1[neighI] == islandI)
                {
                    isVisitedPoint[pointI] = true;
                    addIslandInfoForPoint(pointI, pointI, islandI, 0, pointNormals, prismIslands1, prismIslands2, prismIslands3, pointHops1, pointHops2, pointHops3, pointNormalSource1, pointNormalSource2, pointNormalSource3, pointNormals1, pointNormals2, pointNormals3);
                    ++n;
                    ++nTot;
                    if (isProcessorPoint[pointI])
                    {
                        ++nProcPrisms;
                    }
                    // Pout << "    Added pointI " << pointI << " at " << mesh.points()[pointI] << endl;
                    break;
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
            // Skip faces on processor patches
            const polyPatch& pp = mesh.boundaryMesh()[patchI];
            if (isA<processorPolyPatch>(pp))
                continue;

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
                    addIslandInfoForPoint(pointI, pointNormalSource1[neighI], islandI, 0, pointNormals, prismIslands1, prismIslands2, prismIslands3, pointHops1, pointHops2, pointHops3, pointNormalSource1, pointNormalSource2, pointNormalSource3, pointNormals1, pointNormals2, pointNormals3);
                    ++n;
                    ++nTot;
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
    vectorList& pointNormals3,
    const boolList& isProcessorPoint,
    labelList& islandIs
)
{
    // Processor number (MPI rank)
    const label myProcNo = Pstream::myProcNo();

    // Check for exceeding maximum processes
    const double maxProcs = 1.0 * passiveIndexStart / maxIslands;
    if (myProcNo >= maxProcs)
    {
        FatalError << "Maximum supported number of processes is " << maxProcs << endl << abort(FatalError);
    }

    // Next available island ID
    label islandI = myProcNo * maxIslands;

    // Storage of processed points
    boolList isVisitedPoint(mesh.nPoints(), false);

    // Main identification loop
    label nProcPrisms = 0;
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
        addIslandInfoForPoint(startPointI, startPointI, islandI, 0, pointNormals, prismIslands1, prismIslands2, prismIslands3, pointHops1, pointHops2, pointHops3, pointNormalSource1, pointNormalSource2, pointNormalSource3, pointNormals1, pointNormals2, pointNormals3);
        Pout << "Starting island " << islandI << " pointI " << startPointI << " at " << mesh.points()[startPointI] << endl;

        const label n = propagateIslandInfoOnBoundary(mesh, startPointI, isVisitedPoint, isPrismaticPoint, isLayerSurfacePoint, pointNormals, prismIslands1, prismIslands2, prismIslands3, pointHops1, pointHops2, pointHops3, pointNormalSource1, pointNormalSource2, pointNormalSource3, pointNormals1, pointNormals2, pointNormals3, isProcessorPoint, nProcPrisms);

        if (n < 1)
        {
            break;
        }

        // Reserve next id
        ++islandI;

        if (islandI >= (myProcNo + 1) * maxIslands)
        {
            FatalError << "Exceeded maximum number of islands " << maxIslands << endl << abort(FatalError);
        }
    }

    // Synchronize island IDs among processors

    labelList prismIslands1Sync(mesh.nPoints(), UNDEF_LABEL);
    label oldI = UNDEF_LABEL;
    label newI = UNDEF_LABEL;
    label nRenumbered = 1;
    label i = 0;

    while (nRenumbered > 0)
    {
        nRenumbered = 0;
        oldI = UNDEF_LABEL;
        newI = UNDEF_LABEL;
        ++i;
        Info << "Renumbering round " << i << endl;

        forAll(mesh.points(), pointI)
        {
            prismIslands1Sync[pointI] = prismIslands1[pointI];
        }

        syncTools::syncPointList
            (
                mesh,
                prismIslands1Sync,
                maxEqOp<label>(),
                UNDEF_LABEL               // null value
            );

        // Find old and new index numbers from synced result
        forAll(mesh.points(), pointI)
        {
            if ((prismIslands1[pointI]) < 0)
                continue;
            if (prismIslands1Sync[pointI] > prismIslands1[pointI])
            {
                oldI = prismIslands1[pointI];
                newI = prismIslands1Sync[pointI];
            }
        }

        // Change old index number to new one
        if (newI > 0)
        {
            forAll(mesh.points(), pointI)
            {
                if (prismIslands1[pointI] == oldI)
                {
                    prismIslands1[pointI] = newI;
                    ++nRenumbered;
                }
            }

            Pout << "  Renumbered island " << oldI << " to " << newI << endl;
        }

        // Sync work switch among processors
        const label nRenumberedMax = returnReduce(nRenumbered, maxOp<label>());
        nRenumbered = nRenumberedMax;
    }

    // Find out final island id numbers
    forAll (mesh.points(), pointI)
    {
        if (prismIslands1[pointI] < 0)
            continue;

        if (findIndex(islandIs, prismIslands1[pointI]) == -1)
        {
            islandIs.append(prismIslands1[pointI]);
        }
    }

    // Add the non-prismatic edge points to final islands
    for (const label islandI : islandIs)
    {
        // const label ne =
        addEdgePointsToIsland(mesh, islandI, isPrismaticPoint, isLayerSurfacePoint, pointNormals, prismIslands1, prismIslands2, prismIslands3, pointHops1, pointHops2, pointHops3, pointNormalSource1, pointNormalSource2, pointNormalSource3, pointNormals1, pointNormals2, pointNormals3);

        // Pout << "Island " << islandI << " has " << ne << " edge points" << endl;
    }

    return 0;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Help function to convert active island id number to passive or vice versa

label invertIslandI
(
    const label islandI
)
{
    if (islandI < 0)
    {
        FatalError << "island " << islandI << " is illegal for conversion" << endl;
    }

    if (islandI < passiveIndexStart)
    {
        return passiveIndexStart + islandI;
    }
    else
    {
        return islandI - passiveIndexStart;
    }
}

// Help function to calculate number of free (undefined) islands,
// number of all front islands points, and passive front points
// neighboring each point

int countFrontPoints
(
    const fvMesh& mesh,
    const label islandI,
    labelList& freePoints,
    labelList& frontPoints,
    labelList& passivePoints,
    const labelList& prismIslands1,
    const labelList& prismIslands2,
    const labelList& prismIslands3
)
{
    forAll (mesh.points(), pointI)
    {
        forAll (mesh.pointPoints()[pointI], pointPointI)
        {
            const label neighI = mesh.pointPoints()[pointI][pointPointI];

            // Free points
            if
            (
                (! isPointInIsland(neighI, islandI, prismIslands1, prismIslands2, prismIslands3)) and
                (! isPointInIsland(neighI, invertIslandI(islandI), prismIslands1, prismIslands2, prismIslands3))
            )
            {
                freePoints[pointI] += 1;
                if (pointI==32)
                    Pout << "32 free at " << neighI << endl;
            }

            // All front points
            if
            (
                (isPointInIsland(neighI, islandI, prismIslands1, prismIslands2, prismIslands3)) or
                (isPointInIsland(neighI, invertIslandI(islandI), prismIslands1, prismIslands2, prismIslands3))
            )
            {
                frontPoints[pointI] += 1;
                if (pointI==32)
                    Pout << "32 front at " << neighI
                         << " islandI " <<  isPointInIsland(neighI, islandI, prismIslands1, prismIslands2, prismIslands3)
                         << " invIslandI " <<  isPointInIsland(neighI, invertIslandI(islandI), prismIslands1, prismIslands2, prismIslands3)
                         << " prismIslands1 " << prismIslands1[32]
                         << " prismIslands2 " << prismIslands2[32]
                         << " prismIslands3 " << prismIslands3[32]
                         << endl;
            }

            // Passive front points
            if (isPointInIsland(neighI, invertIslandI(islandI), prismIslands1, prismIslands2, prismIslands3))
            {
                passivePoints[pointI] += 1;
            }
        }
    }

    // Synchronize free points
    syncTools::syncPointList
    (
        mesh,
        freePoints,
        plusEqOp<label>(),
        0                         // null value
    );

    // Synchronize front points
    syncTools::syncPointList
    (
        mesh,
        frontPoints,
        plusEqOp<label>(),
        0                         // null value
    );

    // Synchronize passive points
    syncTools::syncPointList
    (
        mesh,
        passivePoints,
        plusEqOp<label>(),
        0                         // null value
    );

    Pout << "32 free " << freePoints[32] << " front " << frontPoints[32] << " passive " << passivePoints[32] << endl;

    return 0;
}


// Help function to find current front and candidate points for front propagation

int findPropagationFrontPointIs
(
    const fvMesh& mesh,
    const label islandI,
    labelList& frontPointIs,
    labelList& candidatePointIs,
    boolList& isCandidatePrismatic,
    boolList& isPropagationModeActive,
    const labelList& prismIslands1,
    const labelList& prismIslands2,
    const labelList& prismIslands3
)
{
    // label debugV = 0;

    // Count global number of free (undefined) islands, number of
    // front islands and passive island points connected to each point
    labelList freePoints(mesh.nPoints(), Zero);
    labelList frontPoints(mesh.nPoints(), Zero);
    labelList passivePoints(mesh.nPoints(), Zero);
    countFrontPoints(mesh, islandI, freePoints, frontPoints, passivePoints, prismIslands1, prismIslands2, prismIslands3);

    forAll (mesh.points(), pointI)
    {
        // Info << "islandI " << islandI << " pointI " << pointI << endl;

        // Consider only active front points next to
        if (! isPointInIsland(pointI, islandI, prismIslands1, prismIslands2, prismIslands3))
            continue;

        // Find and add pairs of front-to-free points into lists
        forAll (mesh.pointPoints()[pointI], pointPointI)
        {
            const label neighI = mesh.pointPoints()[pointI][pointPointI];

            // Skip neighbor if it is not free
            if
            (
                (isPointInIsland(neighI, islandI, prismIslands1, prismIslands2, prismIslands3)) or
                (isPointInIsland(neighI, invertIslandI(islandI), prismIslands1, prismIslands2, prismIslands3))
            )
            {
                continue;
            }

            if (frontPoints[neighI] > 0)
            {
                // Front point
                frontPointIs.append(pointI);

                // Candidate free point
                candidatePointIs.append(neighI);

                // Is connection prismatic? Yes if front point has
                // only one free point connection (it's the
                // candidate), and the candidate has only one front
                // point connection (it's the front point)
                if ((freePoints[pointI] == 1) and (frontPoints[neighI] == 1))
                {
                    isCandidatePrismatic.append(true);
                }
                else
                {
                    isCandidatePrismatic.append(false);
                }

                // Is propagation for this connection active or
                // passive? Mode is passive if connection is not
                // prismatic, or if the front point is connected to
                // passive front point. Otherwise propagation is
                // active.
                const label lastI = isCandidatePrismatic.size() - 1;
                if
                (
                    (! isCandidatePrismatic[lastI]) or
                    (passivePoints[pointI] > 0)
                )
                {
                    isPropagationModeActive.append(false);
                }
                else
                {
                    isPropagationModeActive.append(true);

                    // // OBJ format debug printout for viewing prismatic edges
                    // if (islandI == 0)
                    // {
                    //     Pout << "v " << mesh.points()[pointI][0] << " " << mesh.points()[pointI][1] << " " << mesh.points()[pointI][2] << endl;
                    //     Pout << "v " << mesh.points()[neighI][0] << " " << mesh.points()[neighI][1] << " " << mesh.points()[neighI][2] << endl;
                    //     Pout << "l " << debugV + 1 << " " << debugV + 2 << endl;
                    //     debugV += 2;
                    // }
                }
            }
        }
    }

    // FatalError << "DEBUG STOP" << endl << abort(FatalError);

    return 0;
}


// Function to run one island info propagation iteration

int propagateIslandFronts
(
    const fvMesh& mesh,
    const label nLayer,
    const labelList& islandIs,
    const vectorList& pointNormals,
    const labelList& nProcessorsOnPoint,
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
    vectorList& pointNormals3,
    labelList& innerPrismPointLabels1,
    labelList& innerPrismPointLabels2,
    labelList& innerPrismPointLabels3,
    labelList& outerPrismPointLabels1,
    labelList& outerPrismPointLabels2,
    labelList& outerPrismPointLabels3
)
{
    scalar nAddedPrisms = 0.0;
    labelList candidatePointIs;
    labelList frontPointIs;
    boolList isCandidatePrismatic;
    boolList isPropagationModeActive;

    for (const label islandI : islandIs)
    {
        // Find next front points (free unassigned point indices next
        // to active island points), and their neighboring active
        // current front point indices
        candidatePointIs.clear();
        frontPointIs.clear();
        isCandidatePrismatic.clear();
        isPropagationModeActive.clear();
        findPropagationFrontPointIs(mesh, islandI, frontPointIs, candidatePointIs, isCandidatePrismatic, isPropagationModeActive, prismIslands1, prismIslands2, prismIslands3);

        // Info << "size of candidatePointIs " << candidatePointIs.size() << endl;
        // Info << "size of frontPointIs " << frontPointIs.size() << endl;
        // Info << "size of isCandidatePrismatic " << isCandidatePrismatic.size() << endl;

        // Process the candidate points
        forAll (candidatePointIs, pointI)
        {
            // Propagate active island index to candidate point only
            // if candidate is prismatic and propagation type is active
            if ((isCandidatePrismatic[pointI]) and (isPropagationModeActive[pointI]))
            {
                addIslandInfoForPoint(candidatePointIs[pointI], frontPointIs[pointI], islandI, nLayer, pointNormals, prismIslands1, prismIslands2, prismIslands3, pointHops1, pointHops2, pointHops3, pointNormalSource1, pointNormalSource2, pointNormalSource3, pointNormals1, pointNormals2, pointNormals3);
                // WIP addPrismaticMappingsForPoint(mesh, candidatePointIs[pointI], frontPointIs[pointI], islandI, prismIslands1, prismIslands2, prismIslands3, innerPrismPointLabels1, innerPrismPointLabels2, innerPrismPointLabels3, outerPrismPointLabels1, outerPrismPointLabels2, outerPrismPointLabels3);
                nAddedPrisms += 1.0 / nProcessorsOnPoint[candidatePointIs[pointI]];
            }

            // Otherwise propagate the passive island index number
            else
            {
                addIslandInfoForPoint(candidatePointIs[pointI], frontPointIs[pointI], invertIslandI(islandI), nLayer, pointNormals, prismIslands1, prismIslands2, prismIslands3, pointHops1, pointHops2, pointHops3, pointNormalSource1, pointNormalSource2, pointNormalSource3, pointNormals1, pointNormals2, pointNormals3);
            }
        }

        // Synchronize propagation variables among processors

        syncTools::syncPointList
        (
            mesh,
            prismIslands1,
            maxEqOp<label>(),
            UNDEF_LABEL           // null value
        );

        syncTools::syncPointList
        (
            mesh,
            prismIslands2,
            maxEqOp<label>(),
            UNDEF_LABEL           // null value
        );

        syncTools::syncPointList
        (
            mesh,
            prismIslands3,
            maxEqOp<label>(),
            UNDEF_LABEL           // null value
        );

        const scalar nSumAddedPrisms = returnReduce(nAddedPrisms, sumOp<scalar>());
        Info << "Layer " << nLayer << " island " << islandI << " added prisms " << nSumAddedPrisms << endl;
    }

    return 0;
}
