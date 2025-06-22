/*---------------------------------------------------------------------------*\
Library
    Orthogonal Boundary Blending

Description
    Special treatment of prismatic boundary layers, with aim to
    increase orthogonality and control the thickness of side edges
    (prismatic edges) on boundary cell layers.
\*---------------------------------------------------------------------------*/

#include "fvMesh.H"
#include <fstream> // for debug printing only (exportEdgesAsStl)

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
// internal mesh points. Calculate also isFlatPatchPoint, which marks
// boundary points whose boundary surfaces all share approximately
// same normal direction.

int calculateBoundaryPointNormals
(
    const fvMesh& mesh,
    pointField& pointNormals,
    boolList& isFlatPatchPoint
)
{
    // Storage for number of boundary faces for points
    labelList nFaces(mesh.nPoints(), 0);

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
            const vector cSf = mesh.Sf()[faceI];
            vector Sf = cSf / mag(cSf);

            const face& f = mesh.faces()[faceI];
            forAll (f, facePointI)
            {
                const label pointI = f[facePointI];
                pointNormals[pointI] -= Sf;
                ++nFaces[pointI];
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

    syncTools::syncPointList
    (
        mesh,
        nFaces,
        plusEqOp<label>(),
        UNDEF_LABEL               // null value
    );

    // Calculate flatness of faces from the length of the accumulated
    // point normals. Note: This is lightweight to calculate, but
    // might not be good enough for star points where an unusual
    // amount of faces meet.

    forAll(pointNormals, pointI)
    {
        if (nFaces[pointI] < 1)
            continue;

        const double idealLength = double(nFaces[pointI]);
        const double magNorm = mag(pointNormals[pointI]);
        const double flatness = magNorm / idealLength;

        // Flatness limit
        if (flatness > 0.9999)
        {
            isFlatPatchPoint[pointI] = true;
            // Info << "pointI " << pointI << " flatness " << flatness << endl;
        }

        // Zero the normal vector for baffle edge points
        if (magNorm < 0.1)
        {
            pointNormals[pointI] = ZERO_VECTOR;
        }
    }

    // Normalise the point normal vectors

    forAll(pointNormals, pointI)
    {
        if (pointNormals[pointI] != ZERO_VECTOR)
        {
            pointNormals[pointI] /= mag(pointNormals[pointI]);
        }
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

// Generates mapping information (isNeighInProc and
// pointToOuterPointMap) for propagation of information to outer
// points. Propagates also the pointNormals vectors from boundary
// points to internal points, to be used for orthogonal alignment of
// internal edges. This is done only for the internal mesh points
// which have a unique shortest edge hop route to one and only one
// boundary point by an edge.

int propagateOuterNeighInfo
(
    const fvMesh& mesh,
    const labelList& patchIds,
    boolList& isNeighInProc,
    labelList& pointToOuterPointMap,
    pointField& pointNormals,
    const labelList& pointHopsToBoundary,
    const label boundaryMaxLayers
)
{
    // Debug option to set true for printing edges as a STL file
    // Best visualized as wireframe in Paraview
    bool exportEdgesAsStl = false;
    std::ofstream myfile;
    if (exportEdgesAsStl)
    {
        myfile.open ("debugEdgesAsStl.stl");
        myfile << "solid edgesAsStl\n";
    }

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
            // hop number, then there is a mapping towards boundary
            if (nNeighHops == 1)
            {
                // If point is a boundary point, then check if the
                // point belongs to a patch eligible for boundary
                // layer treatment. If not, do nothing.
                if (! boundaryPatchCheck(mesh, neighPointI, patchIds, nHops))
                    continue;

                // Mark that the neighbour point is inside this
                // processor domain
                isNeighInProc[pointI] = true;

                // Add index of neighbour to map
                pointToOuterPointMap[pointI] = neighPointI;

                // Copy the point normal from neighbour to this point
                pointNormals[pointI] = pointNormals[neighPointI];

                // Debugging: print edges as triangles in STL ascii format,
                // with fake normal direction
                if (exportEdgesAsStl)
                {
                    myfile << "facet normal 0 0 0" << "\n"
                           << " outer loop" << "\n"
                           << "  vertex "
                           << mesh.points()[pointI][0] << " "
                           << mesh.points()[pointI][1] << " "
                           << mesh.points()[pointI][2] << "\n"
                           << "  vertex "
                           << mesh.points()[neighPointI][0] << " "
                           << mesh.points()[neighPointI][1] << " "
                           << mesh.points()[neighPointI][2] << "\n"
                           << "  vertex "
                           << mesh.points()[pointI][0] + 1e-4 << " "
                           << mesh.points()[pointI][1] + 1e-4 << " "
                           << mesh.points()[pointI][2] + 1e-4 << "\n"
                           << " endloop" << "\n"
                           << "endfacet" << "\n";
                }
            }
        }

        // Synchronize point normals among processors (using
        // maximum magnitude combine). isNeighInProc and
        // pointToOuterPointMap are local info only, so they are
        // not synced.
        syncTools::syncPointList
        (
             mesh,
             pointNormals,
             maxMagSqrEqOp<vector>(),
             UNDEF_VECTOR               // null value
        );
    }

    if (exportEdgesAsStl)
    {
        myfile << "endsolid\n";
        myfile.close();
    }

    return 0;
}

// Generates the mapping required for propagation of inner point
// coordinates to the boundary points

int propagateInnerNeighInfo
(
    const fvMesh& mesh,
    const labelList& patchIds,
    boolList& isNeighInProc,
    labelList& pointToInnerPointMap,
    const labelList& pointHopsToBoundary
)
{
    // Neighbour search is done only for the boundary layer points, as
    // the resulting mapping is used only to smooth boundary points.

    forAll(mesh.points(), pointI)
    {
        // Number of hops this point has to boundary
        const label nHops = pointHopsToBoundary[pointI];

        // Process only the points with a correct hop number
        if (nHops != 0)
            continue;

        // Number of neighbour points with a lower hop count
        label nNeighHops = 0;

        // Index to neighbour point
        label neighPointI = UNDEF_LABEL;

        // Go through all neighbours to find ones with a higher hop count
        forAll(mesh.pointPoints(pointI), pointPpI)
        {
            const label neighI = mesh.pointPoints(pointI)[pointPpI];
            if (pointHopsToBoundary[neighI] == (nHops + 1))
            {
                ++nNeighHops;
                neighPointI = neighI;
            }
        }

        // If exactly one neighbour point was found with a lower
        // hop number, then there is a mapping towards boundary
        if (nNeighHops == 1)
        {

            // If point is a boundary point, then check if the
            // point belongs to a patch eligible for boundary
            // layer treatment. If not, do nothing.
            if (! boundaryPatchCheck(mesh, pointI, patchIds, nHops))
                continue;

            // Mark that the neighbour point is inside this
            // processor domain
            isNeighInProc[pointI] = true;

            // Add index of neighbour to map
            pointToInnerPointMap[pointI] = neighPointI;
        }
    }

    return 0;
}


// Update neighbour coordinates among processors

int updateNeighCoords
(
    const fvMesh& mesh,
    boolList& isNeighInProc,
    labelList& pointToPointMap,
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

        const label neighI = pointToPointMap[pointI];
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
    const boolList& isInternalPoint,
    const labelList& pointHopsToBoundary,
    const pointField& pointNormals,
    const pointField& outerNeighCoords,
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

        const vector outerNeighCoord = outerNeighCoords[pointI];
        if (outerNeighCoord == UNDEF_VECTOR)
            FatalError << "Sanity broken, outerNeighCoord is zero for pointI "
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
        const vector orthoPoint = outerNeighCoord + length * pointNormal;
        const vector blendedPoint = blendFrac * orthoPoint +
            (1.0 - blendFrac) * newPoint;

        // Update point coordinates
        newPoints[pointI] = blendedPoint;
    }

    return 0;
}

// Help function to check if given point index is among points of
// given patch ids

bool isPointOnPatches
(
    const fvMesh& mesh,
    const label pointI,
    const labelList& patchIds
)
{
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

    return false;
}


// Smoothing of boundary points on flat patches. Project the first
// layer prismatic points orthogonally towards boundary, to be used as
// new coordinates for the boundary points.

int projectBoundaryPoints
(
    const fvMesh& mesh,
    pointField& newPoints,
    const boolList& isFlatPatchPoint,
    const labelList& pointHopsToBoundary,
    const pointField& pointNormals,
    const pointField& innerNeighCoords,
    const labelList& patchIds,
    const double boundaryMaxPointBlendingFraction
)
{
    forAll(mesh.points(), pointI)
    {
        const label nHops = pointHopsToBoundary[pointI];
        const vector pointNormal = pointNormals[pointI];
        const vector innerNeighCoord = innerNeighCoords[pointI];

        // Info << "pointI " << pointI << " isFlat " << isFlatPatchPoint[pointI] << " nHops " << nHops << " pointNormal " << pointNormal << " innerNeighCoord " << innerNeighCoord << endl;

        // Skip points without required information
        if (! isFlatPatchPoint[pointI])
            continue;
        if (nHops != 0)
            continue;
        if (pointNormal == ZERO_VECTOR)
            continue;
        if (innerNeighCoord == UNDEF_VECTOR)
            continue;
        if (! isPointOnPatches(mesh, pointI, patchIds))
            continue;

        // Info << "Boundary smoothing for pointI " << pointI << endl;

        // Calculate new coordinates for this boundary point
        const vector cCoords = newPoints[pointI];
        const vector neighVec = cCoords - innerNeighCoord;
        const double dotProd = neighVec & pointNormal;
        const vector pVec = neighVec - dotProd * pointNormal;
        const vector newCoords = cCoords - boundaryMaxPointBlendingFraction * pVec;

        // Update point coordinates
        newPoints[pointI] = newCoords;
    }

    return 0;
}

// ************************************************************************* //
