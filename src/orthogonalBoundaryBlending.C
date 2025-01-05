/*---------------------------------------------------------------------------*\
Library
    Orthogonal Boundary Blending

Description
    Special treatment of boundary cell layer, with aim to increase
    orthogonality of internal side faces on the cell layer.

    Idea is to find internal mesh points which have only one edge
    connection to the boundary points. The ideal direction for those
    edges is opposite of the boundary point vertex normal direction.
    The edge length is calculated by trigonometric projection from
    the current point location to the (inverted) vertex normal vector.
    The resulting point is called Orthogonal point coordinate.

    To avoid intersections caused by non-orthogonal meshes, the new
    point location is a weighted average of the point coordinate from
    main smoothing method and the Orthogonal point coordinate using
    a blending weight factor (orthogonalBlendingFraction).
\*---------------------------------------------------------------------------*/

#include "fvMesh.H"

// Value for initializing lists with an "undefined" value
#define UNDEF_LABEL -1
#define UNDEF_VECTOR vector(GREAT, GREAT, GREAT)
#define ZERO_VECTOR vector(0, 0, 0)

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Calculate the minimum number of edge hops required to reach
// boundary point for all mesh points

int calculatePointHopsToBoundary
(
    const fvMesh& mesh,
    labelList& pointHopsToBoundary
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

    // Propagate distance to boundary to internal mesh points until
    // no more points are processed
    label nProcessedPoints = 1;
    while (nProcessedPoints > 0)
    {
        nProcessedPoints = 0;

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
                ++nProcessedPoints;
                newHopCounts[pointI] = maxHops + 1;
                // Info << "Set pointI " << pointI << " to " << maxHops + 1 << endl;
            }
        }

        // Merge new values with old
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
// polyMesh. Stores point normals to pointNormals field. Point
// normals are not calculated for processor and empty patch points nor
// for internal mesh points.

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
            // const vector Sf = mesh.boundary()[patchI].Sf()[faceI];
            vector Sf = mesh.Sf()[faceI];
            Sf.normalise();
            const face& f = mesh.faces()[faceI];
            forAll (f, facePointI)
            {
                const label pointI = f[facePointI];
                pointNormals[pointI] -= Sf;
                // Info << "Setting pointI " << pointI << " normal " << Sf << endl;
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
        // Info << "pointI " << pointI << " normal " << pointNormals[pointI] << endl;
        if (pointNormals[pointI] != ZERO_VECTOR)
            pointNormals[pointI].normalise();
    }

    return 0;
}


// Calculate label map from internal point to the boundary point, and
// another map from internal point to the neighbour point which is
// along the path to the boundary point. Maps are created only for
// those internal mesh points which have a unique shortest edge hop
// route to one and only one boundary point by an edge.

int calculateBoundaryPointMap
(
    const fvMesh& mesh,
    labelList& pointToBoundaryPointMap,
    labelList& pointToNeighPointMap,
    const labelList& pointHopsToBoundary
)
{
    // Initialize boundary point map so that boundary points point to
    // themselves

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

        for (label faceI = startI; faceI < endI; faceI++)
        {
            const face& f = mesh.faces()[faceI];
            forAll (f, facePointI)
            {
                const label pointI = f[facePointI];
                pointToBoundaryPointMap[pointI] = pointI;
            }
        }
    }

    // Maps are created iteratively to propagate information from
    // boundary towards internal mesh points until all points which
    // can be processed have been processed

    label nProcessedPoints = 1;
    while (nProcessedPoints > 0)
    {
        nProcessedPoints = 0;

        forAll(mesh.points(), pointI)
        {
            // Skip if this point is already processed
            if (pointToBoundaryPointMap[pointI] != UNDEF_LABEL)
                continue;

            // Number of hops this point has to boundary
            const label nHops = pointHopsToBoundary[pointI];

            // Number of neighbour points with a lower hop count
            label nNeighHops = 0;

            // Index to boundary point
            label boundaryPointI = UNDEF_LABEL;

            // Index to neighbour point
            label neighPointI = UNDEF_LABEL;

            // Go through all neighbours with a lower hop count
            forAll(mesh.pointPoints(pointI), pointPpI)
            {
                const label neighI = mesh.pointPoints(pointI)[pointPpI];
                if (pointHopsToBoundary[neighI] == (nHops - 1))
                {
                    ++nNeighHops;
                    boundaryPointI = pointToBoundaryPointMap[neighI];
                    neighPointI = neighI;
                }
            }
            
            // If there is exactly one neighbour point with a lower hop number,
            // then get the mapping to boundary using that neighbour.
            if (nNeighHops == 1)
            {
                ++nProcessedPoints;
                pointToBoundaryPointMap[pointI] = boundaryPointI;
                pointToNeighPointMap[pointI] = neighPointI;
            }
        }
    }
    return 0;
}

// Calculate the point coordinates for the orthogonally optimal point
// location and blend with the given new point coordinates

int blendWithOrthogonalPoints
(
    const polyMesh& mesh,
    pointField& newPoints,
    const bitSet& isInternalPoint,
    const labelList& pointToBoundaryPointMap,
    const labelList& pointToNeighPointMap,
    const labelList& pointHopsToBoundary,
    const pointField& pointNormals,
    const double orthogonalBlendingFraction
)
{
    forAll(mesh.points(), pointI)
    {
        // Skip points without mapping information and boundary points
        if (pointToBoundaryPointMap[pointI] == UNDEF_LABEL)
            continue;
        if (pointHopsToBoundary[pointI] < 1)
            continue;
        if (! isInternalPoint[pointI])
            continue;

        const label boundaryPointI = pointToBoundaryPointMap[pointI];
        if (boundaryPointI < 0)
            FatalError << "Sanity broken, boundaryPointI not defined for pointI "
                       << pointI << endl << abort(FatalError);

        const label neighPointI = pointToNeighPointMap[pointI];
        if (neighPointI < 0)
            FatalError << "Sanity broken, neighPointI not defined for pointI "
                       << pointI << endl << abort(FatalError);

        const label nHops = pointHopsToBoundary[pointI];
        if (nHops < 1)
            FatalError << "Sanity broken, nHops<1 for pointI "
                       << pointI << endl << abort(FatalError);

        const vector pointNormal = pointNormals[boundaryPointI];
        if (pointNormal == ZERO_VECTOR)
            FatalError << "Sanity broken, pointNormal is zero for pointI "
                       << pointI << endl << abort(FatalError);

        // Target length of edge between neighbour and self
        const double length = 0.02 * pow(1.2, nHops);

        // Target blending fraction
        const double blendFrac = orthogonalBlendingFraction * pow(0.8, nHops);

        const vector newPoint = newPoints[pointI];
        const vector neighCoords = mesh.points()[neighPointI];
        const vector orthoPoint = neighCoords + length * pointNormal;
        const vector blendedPoint = blendFrac * orthoPoint +
            (1.0 - blendFrac) * newPoint;

        // Update point coordinates
        newPoints[pointI] = blendedPoint;
    }

    return 0;
}

// ************************************************************************* //
