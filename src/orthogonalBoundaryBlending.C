/*---------------------------------------------------------------------------*\
Library
    Orthogonal Boundary Blending

Description
    Special treatment of prismatic boundary layers, with aim to
    increase orthogonality and control the thickness of side edges
    (prismatic edges) on boundary cell layers.
\*---------------------------------------------------------------------------*/

#include "fvMesh.H"

// Macros for value definitions
#define UNDEF_LABEL -1
#define UNDEF_VECTOR vector(GREAT, GREAT, GREAT)
#define ZERO_VECTOR vector(0, 0, 0)

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Calculate the minimum number of edge hops required to reach
// any boundary point for all mesh points

int calculatePointHopsToBoundary
(
    const fvMesh& mesh,
    labelList& pointHopsToBoundary,
    const label boundaryMaxLayers
)
{
    // Set boundary patch points to zero hops

    forAll(mesh.boundary(), patchI)
    {
        const polyPatch& pp = mesh.boundaryMesh()[patchI];

        // Skip processor and empty patches
        if (isA<processorPolyPatch>(pp))
            continue;
        if (isA<emptyPolyPatch>(pp))
            continue;

        const label startI = mesh.boundary()[patchI].start();
        const label endI = startI + mesh.boundary()[patchI].Cf().size();

        for (label faceI = startI; faceI < endI; faceI++)
        {
            const face& f = mesh.faces()[faceI];
            forAll (f, facePointI)
            {
                const label pointI = mesh.faces()[faceI][facePointI];
                pointHopsToBoundary[pointI] = 0;
            }
        }
    }

    // Storage for new hop counts
    labelList newHopCounts(mesh.nPoints(), -1);

    // Propagate hop information to internal mesh points
    for (label iter = 0; iter < boundaryMaxLayers; ++iter)
    {
        forAll(mesh.points(), pointI)
        {
            // Skip the point if a hop value exists already
            if (pointHopsToBoundary[pointI] >= 0)
                continue;

            // Find the maximum count of neighbour hops
            label maxHops = -1;

            forAll(mesh.pointPoints(pointI), pointPpI)
            {
                const label neighI = mesh.pointPoints(pointI)[pointPpI];
                if (pointHopsToBoundary[neighI] > maxHops)
                    maxHops = pointHopsToBoundary[neighI];
            }

            // If maximum is > 0, then assign maximum + 1 to current point
            if (maxHops >= 0)
            {
                newHopCounts[pointI] = maxHops + 1;
            }
        }

        // Merge new values to hop list
        forAll(mesh.points(), pointI)
        {
            if (newHopCounts[pointI] > pointHopsToBoundary[pointI])
                pointHopsToBoundary[pointI] = newHopCounts[pointI];
        }

        // Synchronize hop list among processors (using max combination)
        syncTools::syncPointList
        (
            mesh,
            pointHopsToBoundary,
            maxEqOp<label>(),
            UNDEF_LABEL               // null value
        );
    }

    return 0;
}

// Calculate point normals of boundary points starting from
// polyMesh. Store point normals to pointNormals field. Point normals
// are not calculated for processor and empty patch points nor for
// internal mesh points.

int calculateBoundaryPointNormals
(
    const fvMesh& mesh,
    pointField& pointNormals
)
{
    forAll(mesh.boundaryMesh(), patchI)
    {
        const polyPatch& pp = mesh.boundaryMesh()[patchI];

        // Skip processor and empty patches
        if (isA<processorPolyPatch>(pp))
            continue;
        if (isA<emptyPolyPatch>(pp))
            continue;

        const label startI = mesh.boundary()[patchI].start();
        const label endI = startI + mesh.boundary()[patchI].Cf().size();

        // Add inversed unit normal vectors of patch faces to all
        // pointNormals of the face points
        for (label faceI = startI; faceI < endI; faceI++)
        {
            // Sf is unit normal vector multiplied by surface area, so
            // need to normalise it before use
            vector Sf = mesh.Sf()[faceI];
            Sf.normalise();

            const face& f = mesh.faces()[faceI];
            forAll (f, facePointI)
            {
                const label pointI = f[facePointI];
                pointNormals[pointI] -= Sf;
            }
        }
    }

    // Synchronize among processors (using sum combination)
    syncTools::syncPointList
    (
        mesh,
        pointNormals,
        plusEqOp<vector>(),
        UNDEF_VECTOR               // null value
    );

    // Normalize the vectors
    forAll(pointNormals, pointI)
    {
        if (pointNormals[pointI] != ZERO_VECTOR)
            pointNormals[pointI].normalise();
    }

    return 0;
}

// A sentinel function to check if a point is a boundary point and not
// part of a patch eligible for boundary layer treatment

bool boundaryPatchCheck
(
    const fvMesh& mesh,
    const label pointI,
    const labelList& patchIds,
    const label nHops
)
{
    // Pass if point is not at boundary
    if (nHops > 1)
        return true;

    // Pass if point is found in an allowed patch
    for (const label patchI : patchIds)
    {
        const polyPatch& pp = mesh.boundaryMesh()[patchI];

        // Skip processor and empty patches
        if (isA<processorPolyPatch>(pp))
            continue;
        if (isA<emptyPolyPatch>(pp))
            continue;

        const label startI = mesh.boundary()[patchI].start();
        const label endI = startI + mesh.boundary()[patchI].Cf().size();

        for (label faceI = startI; faceI < endI; faceI++)
        {
            const face& f = mesh.faces()[faceI];
            forAll (f, facePointI)
            {
                const label testPointI = mesh.faces()[faceI][facePointI];
                if (testPointI == pointI)
                    return true;
            }
        }
    }

    // Fail otherwise
    return false;
}

// Propagate the point normal vectors from boundary points to internal
// points, to be used for orthogonal alignment of internal edges.
// This is done only for the internal mesh points which have a
// unique shortest edge hop route to one and only one boundary point
// by an edge.

int propagateNeighInfo
(
    const fvMesh& mesh,
    const labelList& patchIds,
    bitSet& isNeighInProc,
    labelList& pointToNeighPointMap,
    pointField& pointNormals,
    const labelList& pointHopsToBoundary,
    const label boundaryMaxLayers
)
{
    // Neighbour search is done iteratively to propagate information from
    // boundary towards internal mesh points

    for (label iter = 1; iter < boundaryMaxLayers + 1; ++iter)
    {
        forAll(mesh.points(), pointI)
        {
            // Number of hops this point has to boundary
            const label nHops = pointHopsToBoundary[pointI];

            // Process only the points with a correct hop number
            if (nHops != iter)
                continue;

            // Number of neighbour points with a lower hop count
            label nNeighHops = 0;

            // Index to neighbour point
            label neighPointI = UNDEF_LABEL;

            // Go through all neighbours to find ones with a lower hop count
            forAll(mesh.pointPoints(pointI), pointPpI)
            {
                const label neighI = mesh.pointPoints(pointI)[pointPpI];
                if (pointHopsToBoundary[neighI] == (nHops - 1))
                {
                    ++nNeighHops;
                    neighPointI = neighI;
                }
            }
            
            // If exactly one neighbour point was found with a lower
            // hop number, then there is a mapping to boundary
            if (nNeighHops == 1)
            {
                // If point is a boundary point, then do nothing if
                // the point does not belong to a patch eligible for
                // boundary layer treatment
                if (! boundaryPatchCheck(mesh, neighPointI, patchIds, nHops))
                    continue;

                // Mark that the neighbour point is inside this
                // processor domain
                isNeighInProc.set(pointI);

                // Add index of neighbour to map
                pointToNeighPointMap[pointI] = neighPointI;

                // Copy the point normal from neighbour to this point
                pointNormals[pointI] = pointNormals[neighPointI];
            }
        }

        // Synchronize point normals among processors (using
        // maximum magnitude combine). isNeighInProc and
        // pointToNeighPointMap are local info only, so they are
        // not synced.
        syncTools::syncPointList
        (
             mesh,
             pointNormals,
             maxMagSqrEqOp<vector>(),
             UNDEF_VECTOR               // null value
        );
    }

    return 0;
}

// Update neighbour coordinates among processors

int updateNeighCoords
(
    const fvMesh& mesh,
    bitSet& isNeighInProc,
    labelList& pointToNeighPointMap,
    pointField& neighCoords
)
{
    forAll(mesh.points(), pointI)
    {
        // Reset points whose neighbour is not in processor domain
        if (! isNeighInProc[pointI])
        {
            neighCoords[pointI] = UNDEF_VECTOR;
            continue;
        }

        const label neighI = pointToNeighPointMap[pointI];
        if (neighI < 0)
            FatalError << "Sanity broken, neighI does not exist for pointI "
                       << pointI << endl << abort(FatalError);

        // Save coordinates of the neighbour point
        neighCoords[pointI] = mesh.points()[neighI];
    }

    // Synchronize neighbour coordinates among processors (using
    // min magnitude combine)
    syncTools::syncPointList
    (
         mesh,
         neighCoords,
         minMagSqrEqOp<vector>(),
         UNDEF_VECTOR               // null value
    );

    return 0;
}

// Calculate the point coordinates for the orthogonally optimal point
// location and blend with the given new point coordinates from
// centroidal smoothing

int blendWithOrthogonalPoints
(
    const polyMesh& mesh,
    pointField& newPoints,
    const bitSet& isInternalPoint,
    const labelList& pointHopsToBoundary,
    const pointField& pointNormals,
    const pointField& neighCoords,
    const double boundaryMaxBlendingFraction,
    const double boundaryEdgeLength,
    const double boundaryExpansionRatio,
    const double boundaryMinLayers,
    const double boundaryMaxLayers
)
{
    forAll(mesh.points(), pointI)
    {
        // Skip points without required information
        if (pointNormals[pointI] == ZERO_VECTOR)
            continue;
        if (! isInternalPoint[pointI])
            continue;

        const label nHops = pointHopsToBoundary[pointI];
        if (nHops < 1)
            FatalError << "Sanity broken, nHops<1 for pointI "
                       << pointI << endl << abort(FatalError);

        const vector pointNormal = pointNormals[pointI];
        if (pointNormal == ZERO_VECTOR)
            FatalError << "Sanity broken, pointNormal is zero for pointI "
                       << pointI << endl << abort(FatalError);

        const vector neighCoord = neighCoords[pointI];
        if (neighCoord == UNDEF_VECTOR)
            FatalError << "Sanity broken, neighCoord is zero for pointI "
                       << pointI << endl << abort(FatalError);

        // Target length of edge towards boundary
        const label maxHops = min((nHops - 1), boundaryMaxLayers);
        const double length = boundaryEdgeLength * pow(boundaryExpansionRatio, maxHops);

        // Target blending fraction
        const double slope = -boundaryMaxBlendingFraction / (boundaryMaxLayers - boundaryMinLayers);
        const double y0 = -slope * boundaryMaxLayers;
        const double y = y0 + slope * nHops;
        const double blendFrac = max(0.0, min(y, boundaryMaxBlendingFraction));

        const vector newPoint = newPoints[pointI];
        const vector orthoPoint = neighCoord + length * pointNormal;
        const vector blendedPoint = blendFrac * orthoPoint +
            (1.0 - blendFrac) * newPoint;

        // Update point coordinates
        newPoints[pointI] = blendedPoint;
    }

    return 0;
}

// ************************************************************************* //
