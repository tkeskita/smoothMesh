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
const label MAX_ISLANDS = 10000;

// Passive island index start number
const label PASSIVE_START_LABEL = 1000000000;

// Ignored point island index number. Set to largest used value so
// that parallel synchronization with maxEqOp favors it.
const label IGNORED_LABEL = 2 * PASSIVE_START_LABEL;

// Note: Maximum allowed number of processors is PASSIVE_START_LABEL/MAX_ISLANDS

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
    const fvMesh& mesh,  // TODO: Can remove this after debugging
    const label pointI,
    const label pointNormalSourceI,
    const label islandI,
    const label nLayer,
    labelList& prismIslands1,
    labelList& prismIslands2,
    labelList& prismIslands3,
    labelList& nIslandSlotsUsed,
    labelList& pointHops1,
    labelList& pointHops2,
    labelList& pointHops3,
    labelList& pointNormalSource1,
    labelList& pointNormalSource2,
    labelList& pointNormalSource3
)
{
    if (isPointInIsland(pointI, islandI, prismIslands1, prismIslands2, prismIslands3))
    {
        FatalError << "Point " << pointI << " at " << mesh.points()[pointI] << " is already in island " << islandI << endl << abort(FatalError);
    }

    if (isPointInIsland(pointI, IGNORED_LABEL, prismIslands1, prismIslands2, prismIslands3))
    {
        FatalError << "Point " << pointI << " at " << mesh.points()[pointI] << " is set to ignored state" << endl << abort(FatalError);
    }

    if (prismIslands1[pointI] == UNDEF_LABEL)
    {
        prismIslands1[pointI] = islandI;
        pointNormalSource1[pointI] = pointNormalSourceI;
        pointHops1[pointI] = nLayer;
        nIslandSlotsUsed[pointI] = 1;
        return 1;
    }
    else if (prismIslands2[pointI] == UNDEF_LABEL)
    {
        prismIslands2[pointI] = islandI;
        pointNormalSource2[pointI] = pointNormalSourceI;
        pointHops2[pointI] = nLayer;
        nIslandSlotsUsed[pointI] = 2;
        return 2;
    }
    else if (prismIslands3[pointI] == UNDEF_LABEL)
    {
        prismIslands3[pointI] = islandI;
        pointNormalSource3[pointI] = pointNormalSourceI;
        pointHops3[pointI] = nLayer;
        nIslandSlotsUsed[pointI] = 3;
        return 3;
    }

    // Fail if storage slot runs out
    FatalError << "Maximum limit of three prism islands exceeded for pointI " << pointI << " at " << mesh.points()[pointI] << ". Can't add islandI " << islandI << ". Already assigned islandsIs are " << prismIslands1[pointI] << ", " << prismIslands2[pointI] << ", " << prismIslands3[pointI] << "." << abort(FatalError);

    // Fail if storage slot runs out. TODO: This needs to be improved
    // in such a way that layer propagation is simply not done for any
    // island for that point if the storage slots run out on that
    // round. However, deducing that condition before changing state
    // variables requires thinking through how to program it for
    // parallel processing. Fail for now.

    // FatalError << "Maximum limit of three prism islands exceeded for pointI " << pointI << " at " << mesh.points()[pointI] << ". Can't add islandI " << islandI << ". You need to either lower maxLayers or remove patch(es) from layerPatches, until storage overflow handling is improved." << endl << abort(FatalError);

    // All slots are already in use

    // // Pout << "Reset island pointI " << pointI << " at " << mesh.points()[pointI] << endl;
    // prismIslands1[pointI] = IGNORED_LABEL;
    // prismIslands2[pointI] = IGNORED_LABEL;
    // prismIslands3[pointI] = IGNORED_LABEL;
    // pointNormalSource1[pointI] = IGNORED_LABEL;
    // pointNormalSource2[pointI] = IGNORED_LABEL;
    // pointNormalSource3[pointI] = IGNORED_LABEL;
    // pointNormals1[pointI] = Zero;
    // pointNormals2[pointI] = Zero;
    // pointNormals3[pointI] = Zero;
    // pointHops1[pointI] = IGNORED_LABEL;
    // pointHops2[pointI] = IGNORED_LABEL;
    // pointHops3[pointI] = IGNORED_LABEL;

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
    labelList& nIslandSlotsUsed,
    labelList& pointHops1,
    labelList& pointHops2,
    labelList& pointHops3,
    labelList& pointNormalSource1,
    labelList& pointNormalSource2,
    labelList& pointNormalSource3,
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
        // Pout << "  Starting propagation round for islandI " << islandI << ", nTot=" << nTot << endl;

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

                    if (isPointInIsland(pointI, islandI, prismIslands1, prismIslands2, prismIslands3))
                        continue;

                    addIslandInfoForPoint(mesh, pointI, pointI, islandI, 0, prismIslands1, prismIslands2, prismIslands3, nIslandSlotsUsed, pointHops1, pointHops2, pointHops3, pointNormalSource1, pointNormalSource2, pointNormalSource3);
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


// Help function to add face index for calculation of normal direction

int addFaceNormalInfoForPoint
(
    const fvMesh& mesh,
    const label pointI,
    const label faceNormalSourceI,
    const label islandI,
    labelList& prismIslands1,
    labelList& prismIslands2,
    labelList& prismIslands3,
    labelListList& faceNormalSource1,
    labelListList& faceNormalSource2,
    labelListList& faceNormalSource3
)
{
    if (prismIslands1[pointI] == islandI)
    {
        faceNormalSource1[pointI].append(faceNormalSourceI);
    }
    else if (prismIslands2[pointI] == islandI)
    {
        faceNormalSource2[pointI].append(faceNormalSourceI);
    }
    else if (prismIslands3[pointI] == islandI)
    {
        faceNormalSource3[pointI].append(faceNormalSourceI);
    }
    else
    {
        FatalError << "addFaceNormalInfoForPoint did not find islandI " << islandI << " for pointI " << pointI << endl << abort(FatalError);
    }

    return 0;
}


// Help function to add non-prismatic boundary edge points surrounding
// the island to the island

label addEdgePointsToIsland
(
    const fvMesh& mesh,
    const label islandI,
    const boolList& isPrismaticPoint,
    const boolList& isLayerSurfacePoint,
    const vectorList& pointNormals,
    boolList& isIslandEdgePoint,
    labelList& prismIslands1,
    labelList& prismIslands2,
    labelList& prismIslands3,
    labelList& nIslandSlotsUsed,
    labelList& pointHops1,
    labelList& pointHops2,
    labelList& pointHops3,
    labelList& pointNormalSource1,
    labelList& pointNormalSource2,
    labelList& pointNormalSource3,
    labelListList& faceNormalSource1,
    labelListList& faceNormalSource2,
    labelListList& faceNormalSource3
)
{
    label nTot = 0;
    labelList edgeFaces;

    // Traverse all boundary faces to find island edge faces

    forAll(mesh.boundary(), patchI)
    {
        // Skip faces on processor patches
        const polyPatch& pp = mesh.boundaryMesh()[patchI];
        if (isA<processorPolyPatch>(pp))
            continue;

        const label startI = mesh.boundary()[patchI].start();
        const label endI = startI + mesh.boundary()[patchI].Cf().size();

        for (label faceI = startI; faceI < endI; faceI++)
        {
            const face& f = mesh.faces()[faceI];

            label nIslandPoints = 0;
            label nFreePoints = 0;

            // Count island and free points on the face

            forAll (f, facePointI)
            {
                const label pointI = mesh.faces()[faceI][facePointI];
                if (! isLayerSurfacePoint[pointI])
                    continue;

                if (isPointInIsland(pointI, islandI, prismIslands1, prismIslands2, prismIslands3))
                {
                    ++nIslandPoints;
                }
                else
                {
                    ++nFreePoints;
                }
            }

            // If both free and island points exist, then face needs
            // to be processed

            if ((nIslandPoints > 0) and (nFreePoints > 0))
            {
                edgeFaces.append(faceI);
            }
        }
    }

    // Process edge faces

    forAll (edgeFaces, edgeFaceI)
    {
        const label faceI = edgeFaces[edgeFaceI];
        const face& f = mesh.faces()[faceI];

        forAll (f, facePointI)
        {
            const label pointI = f[facePointI];

            if (! isPointInIsland(pointI, islandI, prismIslands1, prismIslands2, prismIslands3))
            {
                addIslandInfoForPoint(mesh, pointI, UNDEF_LABEL, islandI, 0, prismIslands1, prismIslands2, prismIslands3, nIslandSlotsUsed, pointHops1, pointHops2, pointHops3, pointNormalSource1, pointNormalSource2, pointNormalSource3);

                // Mark as edge point
                isIslandEdgePoint[pointI] = true;

                // Add face index to face normal sources
                addFaceNormalInfoForPoint(mesh, pointI, faceI, islandI, prismIslands1, prismIslands2, prismIslands3, faceNormalSource1, faceNormalSource2, faceNormalSource3);

                ++nTot;
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
    boolList& isIslandEdgePoint,
    labelList& prismIslands1,
    labelList& prismIslands2,
    labelList& prismIslands3,
    labelList& nIslandSlotsUsed,
    labelList& pointHops1,
    labelList& pointHops2,
    labelList& pointHops3,
    labelList& pointNormalSource1,
    labelList& pointNormalSource2,
    labelList& pointNormalSource3,
    labelListList& faceNormalSource1,
    labelListList& faceNormalSource2,
    labelListList& faceNormalSource3,
    const boolList& isProcessorPoint,
    labelList& islandIs
)
{
    // Processor number (MPI rank)
    const label myProcNo = Pstream::myProcNo();

    // Check for exceeding maximum processes
    const double maxProcs = 1.0 * PASSIVE_START_LABEL / MAX_ISLANDS;
    if (myProcNo >= maxProcs)
    {
        FatalError << "Maximum supported number of processes is " << maxProcs << endl << abort(FatalError);
    }

    // Next available island ID
    label islandI = myProcNo * MAX_ISLANDS;

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
        addIslandInfoForPoint(mesh, startPointI, startPointI, islandI, 0, prismIslands1, prismIslands2, prismIslands3, nIslandSlotsUsed, pointHops1, pointHops2, pointHops3, pointNormalSource1, pointNormalSource2, pointNormalSource3);
        // Pout << "Starting island " << islandI << " pointI " << startPointI << " at " << mesh.points()[startPointI] << endl;

        const label n = propagateIslandInfoOnBoundary(mesh, startPointI, isVisitedPoint, isPrismaticPoint, isLayerSurfacePoint, pointNormals, prismIslands1, prismIslands2, prismIslands3, nIslandSlotsUsed, pointHops1, pointHops2, pointHops3, pointNormalSource1, pointNormalSource2, pointNormalSource3, isProcessorPoint, nProcPrisms);
        // Pout << "IslandI " << islandI << " startPointI " << startPointI << " at " << mesh.points()[startPointI] << " nIslandPoints " << n << endl;

        if (n < 1)
        {
            break;
        }

        // Reserve next id
        ++islandI;

        if (islandI >= (myProcNo + 1) * MAX_ISLANDS)
        {
            FatalError << "Exceeded maximum number of islands " << MAX_ISLANDS << endl << abort(FatalError);
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
        Info << "  Island label synchronization iteration " << i << endl;

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

        // Find old and new index numbers from synced result.
        // Only prismIslands1 needs to be checked.
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
            // Pout << "  Renumbered island " << oldI << " to " << newI << endl;
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
        const label ne =
        addEdgePointsToIsland(mesh, islandI, isPrismaticPoint, isLayerSurfacePoint, pointNormals, isIslandEdgePoint, prismIslands1, prismIslands2, prismIslands3, nIslandSlotsUsed, pointHops1, pointHops2, pointHops3, pointNormalSource1, pointNormalSource2, pointNormalSource3, faceNormalSource1, faceNormalSource2, faceNormalSource3);
    }

    return 0;
}


// Help function to synchronize and sort island ids among processors

int mergeAndSortIslandIs
(
    labelList& islandIs
)
{
    // Sort local islands to build stack in correct order
    labelList sortedIs;
    // Foam::sortedOrder(islandIs, sortedIs, maxOp<label>()); // Does not work right, bug in OF?
    Foam::sortedOrder(islandIs, sortedIs);

    // Populate a stack with local island ids
    std::stack<label> islandStack;
    islandStack.push(UNDEF_LABEL);  // terminator label
    for (label i = sortedIs.size() - 1; i >= 0; --i)
    {
        islandStack.push(islandIs[sortedIs[i]]);
    }

    islandIs.clear();

    label firstIsland = 0;

    while (firstIsland != UNDEF_LABEL)
    {
        // Get smallest id among processors
        firstIsland = islandStack.top();
        const label syncedFirst = returnReduce(firstIsland, minOp<label>());

        // Stop if done
        if (syncedFirst == UNDEF_LABEL)
        {
            break;
        }

        // Add to island list and clean up if needed
        islandIs.append(syncedFirst);

        if (firstIsland == syncedFirst)
        {
            islandStack.pop();
        }
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

    if (islandI < PASSIVE_START_LABEL)
    {
        return PASSIVE_START_LABEL + islandI;
    }
    else
    {
        return islandI - PASSIVE_START_LABEL;
    }
}

// Help function to calculate number of free (undefined) islands,
// number of all front islands points, and passive front points
// neighboring each point. If island slots are full for a point, it is
// considered as passive

int countFrontPoints
(
    const fvMesh& mesh,
    const label islandI,
    const labelList& nProcessorsOnPoint,
    labelList& nFreePoints,
    labelList& nFrontPoints,
    labelList& nPassivePoints,
    const labelList& prismIslands1,
    const labelList& prismIslands2,
    const labelList& prismIslands3
)
{
    vectorList freePoints(mesh.nPoints(), UNDEF_VECTOR);
    vectorList frontPoints(mesh.nPoints(), UNDEF_VECTOR);
    vectorList passivePoints(mesh.nPoints(), UNDEF_VECTOR);

    vectorList freePointsSync(mesh.nPoints(), UNDEF_VECTOR);
    vectorList frontPointsSync(mesh.nPoints(), UNDEF_VECTOR);
    vectorList passivePointsSync(mesh.nPoints(), UNDEF_VECTOR);

    forAll (mesh.points(), pointI)
    {
        forAll (mesh.pointPoints()[pointI], pointPointI)
        {
            const label neighI = mesh.pointPoints()[pointI][pointPointI];
            const vector neighP = mesh.points()[neighI];

            // Free points
            if
            (
                (! isPointInIsland(neighI, islandI, prismIslands1, prismIslands2, prismIslands3)) and
                (! isPointInIsland(neighI, invertIslandI(islandI), prismIslands1, prismIslands2, prismIslands3))
            )
            {
                freePoints[pointI] = neighP;
                freePointsSync[pointI] = neighP;
                nFreePoints[pointI] += 1;
            }

            // All front points
            if
            (
                (isPointInIsland(neighI, islandI, prismIslands1, prismIslands2, prismIslands3)) or
                (isPointInIsland(neighI, invertIslandI(islandI), prismIslands1, prismIslands2, prismIslands3))
            )
            {
                frontPoints[pointI] = neighP;
                frontPointsSync[pointI] = neighP;
                nFrontPoints[pointI] += 1;
            }

            // Passive front points
            if
            (
                (isPointInIsland(neighI, invertIslandI(islandI), prismIslands1, prismIslands2, prismIslands3))
            )
            {
                passivePoints[pointI] = neighP;
                passivePointsSync[pointI] = neighP;
                nPassivePoints[pointI] += 1;
            }
        }
    }

    // Synchronize point coordinates

    syncTools::syncPointList
    (
        mesh,
        freePointsSync,
        eqOp<vector>(),
        UNDEF_VECTOR              // null value
    );

    syncTools::syncPointList
    (
        mesh,
        frontPointsSync,
        eqOp<vector>(),
        UNDEF_VECTOR              // null value
    );

    syncTools::syncPointList
    (
        mesh,
        passivePointsSync,
        eqOp<vector>(),
        UNDEF_VECTOR              // null value
    );

    // Check for differences in point coordinates, which indicates
    // that different processors found a different point

    forAll (mesh.points(), pointI)
    {
        if ((freePointsSync[pointI] != UNDEF_VECTOR) and
            (freePoints[pointI] != freePointsSync[pointI]))
        {
            nFreePoints[pointI] += 1;
        }

        if ((frontPointsSync[pointI] != UNDEF_VECTOR) and
            (frontPoints[pointI] != frontPointsSync[pointI]))
        {
            nFrontPoints[pointI] += 1;
        }

        if ((passivePointsSync[pointI] != UNDEF_VECTOR) and
            (passivePoints[pointI] != passivePointsSync[pointI]))
        {
            nPassivePoints[pointI] += 1;
        }
    }

    // Finally get maximum count as the result. Note that this is not
    // a true count for values above one, because not every point
    // coordinate was compared. However this is fine, since only the
    // counts of zero, one and "more than one" are interesting.

    syncTools::syncPointList
    (
        mesh,
        nFreePoints,
        maxEqOp<label>(),
        0                         // null value
    );

    syncTools::syncPointList
    (
        mesh,
        nFrontPoints,
        maxEqOp<label>(),
        0                         // null value
    );

    syncTools::syncPointList
    (
        mesh,
        nPassivePoints,
        maxEqOp<label>(),
        0                         // null value
    );

    return 0;
}


// Help function to find current front and candidate points for front propagation

int findPropagationFrontPointIs
(
    const fvMesh& mesh,
    const label nLayer,
    const label islandI,
    const labelList& nProcessorsOnPoint,
    labelList& frontPointIs,
    labelList& candidatePointIs,
    boolList& isCandidatePrismatic,
    boolList& isPropagationModeActive,
    labelList& frontIslandIs,
    const labelList& prismIslands1,
    const labelList& prismIslands2,
    const labelList& prismIslands3
)
{
    // OBJ output for debugging
    const bool writeObj = false;
    const label debugIsland = 3;
    label debugV = 0;
    std::ofstream myfile;
    if ((writeObj) and (islandI == debugIsland))
    {
        myfile.open("debugPropagationFrontPoints_" + std::to_string(nLayer) + "_proc" + std::to_string(Pstream::myProcNo()) + ".obj");
        myfile << "o obj\n";
    }


    // Count global number of free (undefined) islands, number of
    // front islands and passive island points connected to each point
    labelList nFreePoints(mesh.nPoints(), Zero);
    labelList nFrontPoints(mesh.nPoints(), Zero);
    labelList nPassivePoints(mesh.nPoints(), Zero);
    countFrontPoints(mesh, islandI, nProcessorsOnPoint, nFreePoints, nFrontPoints, nPassivePoints, prismIslands1, prismIslands2, prismIslands3);

    forAll (mesh.points(), pointI)
    {
        // Info << "islandI " << islandI << " pointI " << pointI << endl;
        // Consider only active front points
        if (! isPointInIsland(pointI, islandI, prismIslands1, prismIslands2, prismIslands3))
            continue;

        // Find and add pairs of front-to-free points into lists
        forAll (mesh.pointPoints()[pointI], pointPointI)
        {
            const label neighI = mesh.pointPoints()[pointI][pointPointI];
            // Info << "debug islandI " << islandI << " pointI " << pointI << " neighI " << neighI << endl;

            // Skip neighbor if it's ignored // TODO: Remove?
            if (prismIslands1[neighI] == IGNORED_LABEL)
                continue;

            // Skip neighbor if it is not free
            if
            (
                (isPointInIsland(neighI, islandI, prismIslands1, prismIslands2, prismIslands3)) or
                (isPointInIsland(neighI, invertIslandI(islandI), prismIslands1, prismIslands2, prismIslands3))
            )
            {
                continue;
            }

            if (nFrontPoints[neighI] > 0)
            {
                // Front point
                frontPointIs.append(pointI);

                // Candidate free point
                candidatePointIs.append(neighI);

                // Island number
                frontIslandIs.append(islandI);

                // Is connection prismatic? Yes if front point has
                // only one free point connection (it's the
                // candidate), and the candidate has only one front
                // point connection (it's the front point)
                if ((nFreePoints[pointI] == 1) and (nFrontPoints[neighI] == 1))
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
                    (nPassivePoints[pointI] > 0)
                )
                {
                    isPropagationModeActive.append(false);
                }
                else
                {
                    isPropagationModeActive.append(true);

                    // OBJ format debug printout for viewing prismatic edges
                    if ((writeObj) and (islandI == debugIsland))
                    {
                        myfile << "v " << mesh.points()[pointI][0] << " " << mesh.points()[pointI][1] << " " << mesh.points()[pointI][2] << "\n";
                        myfile << "v " << mesh.points()[neighI][0] << " " << mesh.points()[neighI][1] << " " << mesh.points()[neighI][2] << "\n";
                        myfile << "l " << debugV + 1 << " " << debugV + 2 << "\n";
                        debugV += 2;
                    }
                }
            }
        }
    }

    if ((writeObj) and (islandI == debugIsland))
    {
        myfile.close();
    }
    // FatalError << "DEBUG STOP" << endl << abort(FatalError);

    return 0;
}


// Help function to set a label value in list

void setLabel
(
    labelList& aList,
    const label index,
    const label value,
    const string debugI
)
{
    if (aList[index] == value)
    {
        return;
    }

    if (aList[index] != UNDEF_LABEL)
    {
        FatalError << "List " << debugI << " index " << index << " value is not UNDEF_LABEL but " << aList[index] << ", would have set value to " << value << endl << abort(FatalError);
    }

    aList[index] = value;
}


// Help function to add prismatic mappings for a point pair

int addPrismaticMappingsForPoint
(
    const fvMesh& mesh,
    const label islandI,
    const label frontI,
    const label candidateI,
    const labelList& prismIslands1,
    const labelList& prismIslands2,
    const labelList& prismIslands3,
    labelList& innerPrismPointLabels1,
    labelList& innerPrismPointLabels2,
    labelList& innerPrismPointLabels3,
    labelList& outerPrismPointLabels1,
    labelList& outerPrismPointLabels2,
    labelList& outerPrismPointLabels3
)
{
    // Pout << "islandI " << islandI << " frontI " << frontI << " candidateI " << candidateI << endl;

    // Set outer label
    if (prismIslands1[candidateI] == islandI)
    {
        setLabel(outerPrismPointLabels1, candidateI, frontI, "outerPrismPointLabels1");
    }
    else if (prismIslands2[candidateI] == islandI)
    {
        setLabel(outerPrismPointLabels2, candidateI, frontI, "outerPrismPointLabels2");
    }
    else if (prismIslands3[candidateI] == islandI)
    {
        setLabel(outerPrismPointLabels3, candidateI, frontI, "outerPrismPointLabels3");
    }
    else
    {
        FatalError << "frontI " << frontI << " at " << mesh.points()[frontI] << " is not part of islandI " << islandI << endl << abort(FatalError);
    }

    // Set inner label
    if (prismIslands1[frontI] == islandI)
    {
        setLabel(innerPrismPointLabels1, frontI, candidateI, "innerPrismPointLabels1");
    }
    else if (prismIslands2[frontI] == islandI)
    {
        setLabel(innerPrismPointLabels2, frontI, candidateI, "innerPrismPointLabels2");
    }
    else if (prismIslands3[frontI] == islandI)
    {
        setLabel(innerPrismPointLabels3, frontI, candidateI, "innerPrismPointLabels3");
    }
    else
    {
        FatalError << "candidateI " << candidateI << " at " << mesh.points()[candidateI] << " is not part of islandI " << islandI << endl << abort(FatalError);
    }

    return 0;
}

// Help function to count needed number of slots for points

int countRequiredNumberOfSlots
(
    const fvMesh& mesh,
    const label islandI,
    const labelList& candidatePointIs,
    const labelList& frontIslandIs,
    labelList& nIslandSlotsUsedTest,
    const boolList isDisabledIs
)
{
    forAll(candidatePointIs, pointI)
    {
        if (isDisabledIs[pointI])
            continue;

        if (frontIslandIs[pointI] == islandI)
        {
            const label candidatePointI = candidatePointIs[pointI];
            nIslandSlotsUsedTest[candidatePointI] += 1;
        }
    }

    syncTools::syncPointList
    (
        mesh,
        nIslandSlotsUsedTest,
        maxEqOp<label>(),
        0                     // null value
    );

    return 0;
}

// Help function to disable points if there will be opposing fronts
// meeting, either from previous rounds or this round

int disablePointsForOppositeFronts
(
    const fvMesh& mesh,
    boolList& isDisabledIs,
    const labelList& candidatePointIs,
    const labelList& frontPointIs,
    const labelList& outerPrismPointLabels1,
    const labelList& outerPrismPointLabels2,
    const labelList& outerPrismPointLabels3
)
{
    forAll (frontPointIs, pointI)
    {
        const label frontPointI = frontPointIs[pointI];
        const label candidatePointI = candidatePointIs[pointI];

        // Check for existing front
        if ((outerPrismPointLabels1[frontPointI] == candidatePointI) or
            (outerPrismPointLabels2[frontPointI] == candidatePointI) or
            (outerPrismPointLabels3[frontPointI] == candidatePointI))
        {
            isDisabledIs[pointI] = true;
        }

        forAll (candidatePointIs, pointJ)
        {
            // Check for new front point pairs
            const label frontPointJ = frontPointIs[pointJ];
            const label candidatePointJ = candidatePointIs[pointJ];

            if ((frontPointI == candidatePointJ) and
               (candidatePointI == frontPointJ))
            {
                isDisabledIs[pointI] = true;
                isDisabledIs[pointJ] = true;
            }
        }
    }

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
    labelList& nIslandSlotsUsed,
    labelList& pointHops1,
    labelList& pointHops2,
    labelList& pointHops3,
    labelList& pointNormalSource1,
    labelList& pointNormalSource2,
    labelList& pointNormalSource3,
    labelList& innerPrismPointLabels1,
    labelList& innerPrismPointLabels2,
    labelList& innerPrismPointLabels3,
    labelList& outerPrismPointLabels1,
    labelList& outerPrismPointLabels2,
    labelList& outerPrismPointLabels3
)
{
    label nAddedPrisms = 0;
    labelList candidatePointIs;
    labelList frontPointIs;
    boolList isCandidatePrismatic;
    boolList isPropagationModeActive;
    boolList isVisitedPoint(mesh.nPoints(), false);
    labelList frontIslandIs;

    // Debug option to set true for printing edges as a STL file
    // Best visualized as wireframe in Paraview
    const bool exportEdgesAsStl = true;
    std::ofstream myfile;
    if (exportEdgesAsStl)
    {
        myfile.open("debugPrismEdgesAsStl_" + std::to_string(nLayer) + "_proc" + std::to_string(Pstream::myProcNo()) + ".stl");
        myfile << "solid edgesAsStl\n";
    }

    // Generate lists of prismatic edge information for all islands.
    // Find next front points (free unassigned point indices next
    // to active island points), and their neighboring active
    // current front point indices

    for (const label islandI : islandIs)
    {
        findPropagationFrontPointIs(mesh, nLayer, islandI, nProcessorsOnPoint, frontPointIs, candidatePointIs, isCandidatePrismatic, isPropagationModeActive, frontIslandIs, prismIslands1, prismIslands2, prismIslands3);
    }

    // Disable points if there will be opposing fronts meeting, either
    // from previous rounds or this round
    boolList isDisabledIs(candidatePointIs.size(), false);
    disablePointsForOppositeFronts(mesh, isDisabledIs, candidatePointIs, frontPointIs, outerPrismPointLabels1, outerPrismPointLabels2, outerPrismPointLabels3);

    // Disable points if there will be lack of island slots
    labelList nIslandSlotsUsedTest(nIslandSlotsUsed);
    for (const label islandI : islandIs)
    {
        countRequiredNumberOfSlots(mesh, islandI, candidatePointIs, frontIslandIs, nIslandSlotsUsedTest, isDisabledIs);
    }

    boolList isDisabledPropagationPoint(mesh.nPoints(), false);
    forAll(mesh.points(), pointI)
    {
        if (nIslandSlotsUsedTest[pointI] > 3)
        {
            isDisabledPropagationPoint[pointI] = true;
        }
    }

    // Carry out propagation, island by island to keep island info synced

    for (const label islandI : islandIs)
    {
        forAll (candidatePointIs, pointI)
        {
            const label frontPointI = frontPointIs[pointI];
            const label candidatePointI = candidatePointIs[pointI];
            const label frontIslandI = frontIslandIs[pointI];

            if (frontIslandI != islandI)
                continue;
            //if (isDisabledPropagationPoint[frontPointI])
            //    continue;
            if (isDisabledPropagationPoint[candidatePointI])
                continue;
            if (isDisabledIs[pointI])
                continue;

            // Propagate active island index to candidate point only
            // if candidate is prismatic and propagation type is active
            if ((isCandidatePrismatic[pointI]) and (isPropagationModeActive[pointI]))
            {
                addIslandInfoForPoint(mesh, candidatePointI, frontPointI, frontIslandI, nLayer, prismIslands1, prismIslands2, prismIslands3, nIslandSlotsUsed, pointHops1, pointHops2, pointHops3, pointNormalSource1, pointNormalSource2, pointNormalSource3);

                addPrismaticMappingsForPoint(mesh, frontIslandI, frontPointI, candidatePointI, prismIslands1, prismIslands2, prismIslands3, innerPrismPointLabels1, innerPrismPointLabels2, innerPrismPointLabels3, outerPrismPointLabels1, outerPrismPointLabels2, outerPrismPointLabels3);
                ++nAddedPrisms;
            }

            // Otherwise propagate the passive island index number
            // unless it's already propagated
            else
            {
                if (isPointInIsland(candidatePointI, invertIslandI(frontIslandI), prismIslands1, prismIslands2, prismIslands3))
                    continue;

                addIslandInfoForPoint(mesh, candidatePointI, frontPointI, invertIslandI(frontIslandI), nLayer, prismIslands1, prismIslands2, prismIslands3, nIslandSlotsUsed, pointHops1, pointHops2, pointHops3, pointNormalSource1, pointNormalSource2, pointNormalSource3);
            }

            // Debugging: export edge as a sliver triangle in STL ascii format,
            // with fake normal direction.
            if (exportEdgesAsStl)
            {
                const label i1 = candidatePointI;
                const label i2 = frontPointI;
                if (i1 == i2)
                    FatalError << "candidate point " << i1 << " is same as front point" << endl << abort(FatalError);
                myfile << "facet normal 0 0 0" << "\n"
                       << " outer loop" << "\n"
                       << "  vertex "
                       << mesh.points()[i1][0] << " "
                       << mesh.points()[i1][1] << " "
                       << mesh.points()[i1][2] << "\n"
                       << "  vertex "
                       << mesh.points()[i2][0] << " "
                       << mesh.points()[i2][1] << " "
                       << mesh.points()[i2][2] << "\n"
                       << "  vertex "
                       << mesh.points()[i1][0] * (1.0 + ABS_TOL) + ABS_TOL << " "
                       << mesh.points()[i1][1] * (1.0 + ABS_TOL) + ABS_TOL << " "
                       << mesh.points()[i1][2] * (1.0 + ABS_TOL) + ABS_TOL << "\n"
                       << " endloop" << "\n"
                       << "endfacet" << "\n";
            }
        }

        // Synchronize propagation variables among processors

        syncTools::syncPointList
        (
            mesh,
            nIslandSlotsUsed,
            maxEqOp<label>(),
            UNDEF_LABEL           // null value
        );

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

        //

        syncTools::syncPointList
        (
            mesh,
            pointHops1,
            maxEqOp<label>(),
            UNDEF_LABEL           // null value
        );

        syncTools::syncPointList
        (
            mesh,
            pointHops2,
            maxEqOp<label>(),
            UNDEF_LABEL           // null value
        );

        syncTools::syncPointList
        (
            mesh,
            pointHops3,
            maxEqOp<label>(),
            UNDEF_LABEL           // null value
        );

        // Info << "Layer " << nLayer << " island " << islandI << " added prisms " << nSumAddedPrisms << endl;
    }

    if (exportEdgesAsStl)
    {
        myfile << "endsolid\n";
        myfile.close();
    }

    const label nSumAddedPrisms = returnReduce(nAddedPrisms, sumOp<label>());
    return nSumAddedPrisms;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Help function to update neighboring point vector property to
// points. Used in the propagation of e.g. locations of either inner
// or outer points, depending on which label and vector lists are
// given as arguments

int updatePointVectorValues
(
    const fvMesh& mesh,
    const vectorList& vectorValues,
    const labelList& pointLabels1,
    const labelList& pointLabels2,
    const labelList& pointLabels3,
    vectorList& points1,
    vectorList& points2,
    vectorList& points3
)
{
    // Set current values
    forAll(vectorValues, pointI)
    {
        if (pointLabels1[pointI] == UNDEF_LABEL)
        {
            points1[pointI] = UNDEF_VECTOR;
        }
        else
        {
            points1[pointI] = vectorValues[pointLabels1[pointI]];
        }

        if (pointLabels2[pointI] == UNDEF_LABEL)
        {
            points2[pointI] = UNDEF_VECTOR;
        }
        else
        {
            points2[pointI] = vectorValues[pointLabels2[pointI]];
        }

        if (pointLabels3[pointI] == UNDEF_LABEL)
        {
            points3[pointI] = UNDEF_VECTOR;
        }
        else
        {
            points3[pointI] = vectorValues[pointLabels3[pointI]];
        }
    }

    // Synchronize
    syncTools::syncPointList
    (
        mesh,
        points1,
        minMagSqrEqOp<vector>(),
        UNDEF_VECTOR           // null value
    );

    syncTools::syncPointList
    (
        mesh,
        points2,
        minMagSqrEqOp<vector>(),
        UNDEF_VECTOR           // null value
    );

    syncTools::syncPointList
    (
        mesh,
        points3,
        minMagSqrEqOp<vector>(),
        UNDEF_VECTOR           // null value
    );

    return 0;
}


// Help function to set potentially updated boundary point normals
// from pointNormals to island slots. Point normals can change if
// boundary patches deform during smoothing.

int updateBoundaryPointNormals
(
    const fvMesh& mesh,
    const boolList& isIslandEdgePoint,
    const vectorList& pointNormals,
    const labelList& pointHops1,
    const labelList& pointHops2,
    const labelList& pointHops3,
    const labelListList& faceNormalSource1,
    const labelListList& faceNormalSource2,
    const labelListList& faceNormalSource3,
    vectorList& pointNormals1,
    vectorList& pointNormals2,
    vectorList& pointNormals3
)
{
    forAll(mesh.points(), pointI)
    {
        if (pointHops1[pointI] == 0)
        {
            if (isIslandEdgePoint[pointI])
            {
                pointNormals1[pointI] = Zero;
                forAll(faceNormalSource1[pointI], i)
                {
                    pointNormals1[pointI] += getBoundaryNf(mesh, faceNormalSource1[pointI][i]);
                }
            }
            else
            {
                pointNormals1[pointI] = pointNormals[pointI];
            }
        }

        if (pointHops2[pointI] == 0)
        {
            if (isIslandEdgePoint[pointI])
            {
                pointNormals2[pointI] = Zero;
                forAll(faceNormalSource2[pointI], i)
                {
                    pointNormals2[pointI] += getBoundaryNf(mesh, faceNormalSource2[pointI][i]);
                }
            }
            else
            {
                pointNormals2[pointI] = pointNormals[pointI];
            }
        }

        if (pointHops3[pointI] == 0)
        {
            if (isIslandEdgePoint[pointI])
            {
                pointNormals3[pointI] = Zero;
                forAll(faceNormalSource3[pointI], i)
                {
                    pointNormals3[pointI] += getBoundaryNf(mesh, faceNormalSource3[pointI][i]);
                }
            }
            else
            {
                pointNormals3[pointI] = pointNormals[pointI];
            }
        }
    }

    // Synchronize

    syncTools::syncPointList
    (
        mesh,
        pointNormals1,
        plusEqOp<vector>(),
        UNDEF_VECTOR           // null value
    );

    syncTools::syncPointList
    (
        mesh,
        pointNormals2,
        plusEqOp<vector>(),
        UNDEF_VECTOR           // null value
    );

    syncTools::syncPointList
    (
        mesh,
        pointNormals3,
        plusEqOp<vector>(),
        UNDEF_VECTOR           // null value
    );

    // Normalize

    forAll(mesh.points(), pointI)
    {
        if (mag(pointNormals1[pointI]) > 0.0)
            pointNormals1[pointI] /= mag(pointNormals1[pointI]);

        if (mag(pointNormals2[pointI]) > 0.0)
            pointNormals2[pointI] /= mag(pointNormals2[pointI]);

        if (mag(pointNormals3[pointI]) > 0.0)
            pointNormals3[pointI] /= mag(pointNormals3[pointI]);
    }

    return 0;
}


// Help function to return vector corresponding to a given island at a
// point

vector getSameIslandVectorFromPoint
(
    const label pointI,
    const label islandI,
    const labelList& prismIslands1,
    const labelList& prismIslands2,
    const labelList& prismIslands3,
    const vectorList& points1,
    const vectorList& points2,
    const vectorList& points3
)
{
    if (prismIslands1[pointI] == islandI)
    {
        return points1[pointI];
    }
    else if (prismIslands2[pointI] == islandI)
    {
        return points2[pointI];
    }
    else if (prismIslands3[pointI] == islandI)
    {
        return points3[pointI];
    }
    else
    {
        FatalError << "getSameIslandVectorFromPoint failed to find islandI " << islandI << " for pointI " << pointI << endl << abort(FatalError);
    }

    return UNDEF_VECTOR;
}


// Help function to propagate neighboring point vector property to
// points. Propagation is done for points in given hop order.

int propagatePointVectorValues
(
    const fvMesh& mesh,
    const labelList& hopOrder,
    const labelList& prismIslands1,
    const labelList& prismIslands2,
    const labelList& prismIslands3,
    const labelList& pointHops1,
    const labelList& pointHops2,
    const labelList& pointHops3,
    const labelList& pointLabels1,
    const labelList& pointLabels2,
    const labelList& pointLabels3,
    vectorList& points1,
    vectorList& points2,
    vectorList& points3
)
{
    for (const label nHops : hopOrder)
    {
        // Set current values
        forAll(mesh.points(), pointI)
        {
            if (pointHops1[pointI] == nHops)
            {
                if (pointLabels1[pointI] == UNDEF_LABEL)
                {
                    points1[pointI] = UNDEF_VECTOR;
                }
                else
                {
                    const label islandI = prismIslands1[pointI];
                    const label sourceI = pointLabels1[pointI];
                    points1[pointI] = getSameIslandVectorFromPoint(sourceI, islandI, prismIslands1, prismIslands2, prismIslands3, points1, points2, points3);
                }
            }

            if (pointHops2[pointI] == nHops)
            {
                if (pointLabels2[pointI] == UNDEF_LABEL)
                {
                    points2[pointI] = UNDEF_VECTOR;
                }
                else
                {
                    const label islandI = prismIslands2[pointI];
                    const label sourceI = pointLabels2[pointI];
                    points2[pointI] = getSameIslandVectorFromPoint(sourceI, islandI, prismIslands1, prismIslands2, prismIslands3, points1, points2, points3);
                }
            }

            if (pointHops3[pointI] == nHops)
            {
                if (pointLabels3[pointI] == UNDEF_LABEL)
                {
                    points3[pointI] = UNDEF_VECTOR;
                }
                else
                {
                    const label islandI = prismIslands3[pointI];
                    const label sourceI = pointLabels3[pointI];
                    points3[pointI] = getSameIslandVectorFromPoint(sourceI, islandI, prismIslands1, prismIslands2, prismIslands3, points1, points2, points3);
                }
            }
        }
    }

    // Synchronize
    syncTools::syncPointList
    (
        mesh,
        points1,
        minMagSqrEqOp<vector>(),
        UNDEF_VECTOR           // null value
    );

    syncTools::syncPointList
    (
        mesh,
        points2,
        minMagSqrEqOp<vector>(),
        UNDEF_VECTOR           // null value
    );

    syncTools::syncPointList
    (
        mesh,
        points3,
        minMagSqrEqOp<vector>(),
        UNDEF_VECTOR           // null value
    );

    return 0;
}


// Main help function to update boundary point normals and propagate point
// normals

int updateAndPropagatePointNormals
(
    const fvMesh& mesh,
    const label maxLayers,
    const boolList& isIslandEdgePoint,
    const vectorList& pointNormals,
    const labelList& prismIslands1,
    const labelList& prismIslands2,
    const labelList& prismIslands3,
    const labelList& pointHops1,
    const labelList& pointHops2,
    const labelList& pointHops3,
    const labelList& pointNormalSource1,
    const labelList& pointNormalSource2,
    const labelList& pointNormalSource3,
    const labelListList& faceNormalSource1,
    const labelListList& faceNormalSource2,
    const labelListList& faceNormalSource3,
    vectorList& pointNormals1,
    vectorList& pointNormals2,
    vectorList& pointNormals3
)
{
    // Update new point normals for boundary points
    updateBoundaryPointNormals(mesh, isIslandEdgePoint, pointNormals, pointHops1, pointHops2, pointHops3, faceNormalSource1, faceNormalSource2, faceNormalSource3, pointNormals1, pointNormals2, pointNormals3);

    // Propagate point normals towards inner mesh
    labelList hopOrder;
    for (label i = maxLayers; i > 0; i--)
    {
        hopOrder.append(i);
    }

    propagatePointVectorValues(mesh, hopOrder, prismIslands1, prismIslands2, prismIslands3, pointHops1, pointHops2, pointHops3, pointNormalSource1, pointNormalSource2, pointNormalSource3, pointNormals1, pointNormals2, pointNormals3);

    return 0;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// Help function to calculate a line-plane intersecion point. The line
// is defined with first point to second point, and the plane is a
// normal plane, located at a given distance from the second point
// towards the normal direction. Using equation for "d" for the
// algebraic form from
// https://en.wikipedia.org/wiki/Line-plane_intersection

point calcPointSlideOnLine
(
    const point p1,
    const point p2,
    const double distance,
    const vector normalVec
)
{
    if (mag(p2 - p1) < ABS_TOL)
    {
        FatalError << "Points " << p1 << " and " << p2 << " are too close to each other" << endl << abort(FatalError);
    }

    const vector p0 = p2 + distance * normalVec;
    const double c1 = (p0 - p1) & normalVec;
    const vector lVec = (p2 - p1) / mag(p2 - p1);
    const double c2 = lVec & normalVec;
    if (fabs(c2) < ABS_TOL)
    {
        FatalError << "Line vector for points " << p1 << " and " << p2 << ", and normal vector " << normalVec << " are almost the same, c2 is " << c2 << endl << abort(FatalError);
    }

    const vector target = p1 + (c1 / c2) * lVec;
    return target;
}

// Help function to calculate a layer point location for one island
// slot

vector calcLayerPointMove
(
    const fvMesh& mesh,
    const label pointI,
    const labelList& pointHops,
    const vectorList& outerPrismPoints,
    const vectorList& pointNormals,
    const double layerEdgeLength,
    const double layerExpansionRatio,
    const label debugI
)
{
    if (pointNormals[pointI] == UNDEF_VECTOR)
    {
        FatalError
            << "Sanity broken, pointNormals" << debugI << " is zero for pointI"
            << pointI << " at " << mesh.points()[pointI] << endl << abort(FatalError);
    }

    if (pointHops[pointI] < 1)
    {
        FatalError
            << "Sanity broken, pointHops" << debugI << " is " << pointHops[pointI]
            << " for pointI " << pointI << " at " << mesh.points()[pointI]
            << endl << abort(FatalError);
    }

    // Target length of edge towards boundary
    const label nHops = pointHops[pointI] - 1;
    const double layerThickness = layerEdgeLength * pow(layerExpansionRatio, nHops);

    const point movingPoint = mesh.points()[pointI];
    const point refPoint = outerPrismPoints[pointI];
    const vector endNormal = pointNormals[pointI];

    if (fabs(mag(endNormal) - 1.0) > ABS_TOL)
    {
        FatalError << "pointI " << pointI << " normal vector " << endNormal << " length is not unity" << endl << abort(FatalError);
    }

    const vector layerTargetPoint = calcPointSlideOnLine(movingPoint, refPoint, layerThickness, endNormal);

    return layerTargetPoint;
}


// Function to calculate new coordinates according to layer treatment

int blendWithLayerPoints
(
    const fvMesh& mesh,
    pointField& newPoints,
    const boolList& isInternalPoint,
    const labelList& pointHops1,
    const labelList& pointHops2,
    const labelList& pointHops3,
    const vectorList& outerPrismPoints1,
    const vectorList& outerPrismPoints2,
    const vectorList& outerPrismPoints3,
    const vectorList& pointNormals1,
    const vectorList& pointNormals2,
    const vectorList& pointNormals3,
    const double layerMaxBlendingFraction,
    const double layerEdgeLength,
    const double layerExpansionRatio,
    const double minLayers,
    const double maxLayers
)
{
    const bool writeCsv = true;
    std::ofstream myfile;
    if (writeCsv)
    {
        myfile.open("debugLayerPoints.csv");
        myfile << "x,y,z\n";
    }

    forAll(mesh.points(), pointI)
    {
        label n = 0;
        vector newCoords = ZERO_VECTOR;
        labelList nHops;

        // Calculate point movements for each slot and add to newCoords

        if (outerPrismPoints1[pointI] != UNDEF_VECTOR)
        {
            const vector newP = calcLayerPointMove(mesh, pointI, pointHops1, outerPrismPoints1, pointNormals1, layerEdgeLength, layerExpansionRatio, 1);
            newCoords += newP;
            ++n;
            nHops.append(pointHops1[pointI]);
        }

        if (outerPrismPoints2[pointI] != UNDEF_VECTOR)
        {
            const vector newP = calcLayerPointMove(mesh, pointI, pointHops2, outerPrismPoints2, pointNormals2, layerEdgeLength, layerExpansionRatio, 2);
            newCoords += newP;
            ++n;
            nHops.append(pointHops2[pointI]);
        }

        if (outerPrismPoints3[pointI] != UNDEF_VECTOR)
        {
            const vector newP = calcLayerPointMove(mesh, pointI, pointHops3, outerPrismPoints3, pointNormals3, layerEdgeLength, layerExpansionRatio, 3);
            newCoords += newP;
            ++n;
            nHops.append(pointHops3[pointI]);
        }

        if (n > 0)
        {
            // New point coordinates from layer treatment
            const point layerPoint = newCoords / n;

            if (writeCsv)
                myfile << layerPoint[0] << ","
                     << layerPoint[1] << ","
                     << layerPoint[2] << "\n";

            // Artificially slow down smoothing of boundary points, to
            // get internal points to smoothen faster. This reduces
            // squishing of internal cells, but is not good enough.
            // if (! isInternalPoint[pointI])
            // {
            //     const double fac = 0.5;  // blending factor
            //     layerPoint = (1.0 - fac) * newPoints[pointI] + fac * layerPoint;
            // }

            // Target blending fraction
            const double slope = -layerMaxBlendingFraction / (maxLayers + 1.0 - minLayers);
            const double y0 = -slope * (maxLayers + 1.0);
            const double y = y0 + slope * min(nHops); // CHECKME: is nHops 0 OK?
            const double blendFrac = max(0.0, min(y, layerMaxBlendingFraction));

            const point newPoint = newPoints[pointI];
            const vector blendedPoint = blendFrac * layerPoint +
                (1.0 - blendFrac) * newPoint;

            // Update point coordinates
            newPoints[pointI] = blendedPoint;
        }
    }

    if (writeCsv)
    {
        myfile.close();
    }

    return 0;
}
