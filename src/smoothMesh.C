/*---------------------------------------------------------------------------*\
Application
    smoothMesh

Description
    Smooth internal mesh points to improve mesh quality
\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "polyMesh.H"
#include "Time.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "ReadFields.H"
#include "syncTools.H"
#include "meshTools.H"
#include <bits/stdc++.h> // For std::stack
#include "vectorList.H"
#include "stringListOps.H"  // For stringListOps::findMatching()
#include "ListOps.H"
#include "wordReList.H"
#include "triSurface.H"
#include "triSurfaceSearch.H"
#include "edgeMesh.H"
#include <fstream> // for fileExists

// Macros for value definitions
#define UNDEF_LABEL -1
#define UNDEF_VECTOR vector(GREAT, GREAT, GREAT)
#define ZERO_VECTOR vector(0, 0, 0)

// Boolean for developer mode. Set to false to use bleeding edge
// work in progress features.
// #define USE_STABLE_FEATURES_ONLY true

// #include <typeinfo>
// Typeinfo is needed only for getting types while debugging, for example:
// Info << "Type is " << typeid(x).name() << endl;
// Use terminal command like this to demangle the mangled name:
// c++filt -t N4Foam4faceE

#include "orthogonalBoundaryBlending.C"
#include "boundaryPointSmoothing.C"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Help function to find internal points of the argument mesh.
// Updates argument boolList accordingly: true for internal points
// (including processor points) and false for boundary points.

int findInternalMeshPoints
(
    const fvMesh& mesh,
    boolList& isInternalPoint
)
{
    // Start from all points in
    isInternalPoint = true;

    // Remove points on boundary patches, except not
    // processor patches nor empty patches
    const faceList& faces = mesh.faces();
    forAll(mesh.boundary(), patchI)
    {
        const polyPatch& pp = mesh.boundaryMesh()[patchI];

        // Skip processor and empty patches
        if (isA<processorPolyPatch>(pp))
            continue;

        // No smoothing for non-3D meshes
        if (isA<emptyPolyPatch>(pp))
        {
            FatalError
                << "Smoothing of non-3D meshes (meshes with type empty patches) is not supported"
                << endl << abort(FatalError);
        }

        const label startI = mesh.boundary()[patchI].start();
        const label endI = startI + mesh.boundary()[patchI].Cf().size();

        for (label faceI = startI; faceI < endI; faceI++)
        {
            const face& f = mesh.faces()[faceI];
            forAll (f, pointI)
            {
                const label i = faces[faceI][pointI];
                isInternalPoint[i] = false;
            }
        }
    }

    label nPoints = 0;
    forAll(isInternalPoint, pointI)
    {
        if (isInternalPoint[pointI])
            ++nPoints;
    }

    const label sumNPoints = returnReduce(nPoints, sumOp<label>());
    return sumNPoints;
}


// Function for centroidal smoothing of internal mesh points.

Foam::tmp<Foam::pointField> centroidalSmoothing
(
    const fvMesh& mesh,
    const boolList isInternalPoint,
    const bool doBoundarySmoothing
)
{
    // Centroidal smoothing algorithm

    // Calculate number and sum of surrounding cell center
    // coordinates

    vectorField cellPoints(mesh.points().size(), Zero);

    labelList nPoints(mesh.points().size(), Zero);

    forAll(mesh.points(), pointI)
    {
        // Carry out centroidal smoothing only for internal points
        // if boundary point smoothing is disabled
        if ((! doBoundarySmoothing) and (! isInternalPoint[pointI]))
        {
            continue;
        }

        const labelList& pCells = mesh.pointCells(pointI);

        // Number of points
        nPoints[pointI] = pCells.size();

        // Sum of cell center coordinates
        for (const label celli : pCells)
        {
            cellPoints[pointI] += mesh.cellCentres()[celli];
        }
    }

    // Synchronize among processors
    syncTools::syncPointList
    (
        mesh,
        cellPoints,
        plusEqOp<vector>(),
        UNDEF_VECTOR               // null value
    );

    syncTools::syncPointList
    (
        mesh,
        nPoints,
        plusEqOp<label>(),
        label(0)                   // null value
    );

    // Calculate new point locations
    // - init with fallback value (the current mesh points)
    tmp<pointField> tnewPoints(new pointField(mesh.points()));
    auto& newPoints = tnewPoints.ref();

    forAll(newPoints, pointI)
    {
        // centroidal point
        if (nPoints[pointI])
        {
            newPoints[pointI] =
                cellPoints[pointI] / nPoints[pointI];
        }
    }

    return tnewPoints;
}


// Help function to return distance to current mesh point with index
// pointI from argument coords.

double getPointDistance
(
    const vector coords1,
    const vector coords2
)
{
    const vector v = coords2 - coords1;
    return mag(v);
}


/////////////////////////////////////////////
// Help functions for aspectRatioSmoothing //
/////////////////////////////////////////////

// Generate a lists of neighbour points from pointPoints where the
// point pair share a cell. Used in aspectRatioSmoothing.

int generatePointNeighPoints
(
    const fvMesh& mesh,
    labelListList& pointNeighPoints
)
{
    forAll (mesh.points(), pointI)
    {
        forAll (mesh.pointCells(pointI), pointCellI)
        {
            const label cellI = mesh.pointCells(pointI)[pointCellI];
            forAll (mesh.cellPoints(cellI), pointCPI)
            {
                const label pointPointI = mesh.cellPoints(cellI)[pointCPI];
                if (pointI == pointPointI)
                    continue;

                // Add if not already added
                if (findIndex(pointNeighPoints[pointI], pointPointI) == -1)
                {
                    pointNeighPoints[pointI].append(pointPointI);
                }
            }
        }
    }

    return 0;
}

// Help function which compares values of vectors element by element
// and returns true if first argument is smaller than second argument.

bool isSmallerByVectorElements
(
    const vector vec1,
    const vector vec2
)
{
    if (vec1.size() != vec2.size())
        FatalError << vec1 << " size does not equal size of " << vec2 << endl << abort(FatalError);

    forAll(vec1, i)
    {
        if (vec1[i] < vec2[i])
            return true;
        else if (vec1[i] > vec2[i])
            return false;
    }
    return false;
}


// Help function to determine if point1 is closer to origin than
// point2. This function takes into account the special case where
// distance is equal but point coordinates are different.

bool isCloserPoint
(
    const vector point1,
    const vector point2
)
{
    if (point1 == point2)
    {
        return false;
    }

    // point1 is clearly closer than point2
    const double deltaDistance = mag(point1) - mag(point2);
    if (deltaDistance < VSMALL)
    {
        return true;
    }

    // point1 and point2 are at same distance, prefer the point with
    // smaller x, y, or z coordinate
    else if ((abs(deltaDistance) < VSMALL) and (isSmallerByVectorElements(point1, point2)))
    {
        return true;
    }

    return false;
}


// Help function for aspectRatioSmoothing to find out closest
// three points (connected by edges) to each point

int findClosestPoints
(
    const fvMesh& mesh,
    const boolList& isInternalPoint,
    vectorList& closestPoints1,
    vectorList& closestPoints2,
    vectorList& closestPoints3,
    boolList& hasCommonCell,
    const labelListList& pointNeighPoints
)
{
    // Initialize with local information
    forAll (mesh.points(), pointI)
    {
        if (! isInternalPoint[pointI])
        {
            closestPoints1[pointI] = ZERO_VECTOR;
            closestPoints2[pointI] = ZERO_VECTOR;
            closestPoints3[pointI] = ZERO_VECTOR;
            hasCommonCell[pointI] = false;
            continue;
        }

        const vector cCoords = mesh.points()[pointI];
        const labelList& pointPoints = mesh.pointPoints()[pointI];
        const label nPoints = pointPoints.size();
        List<scalar> edgeLengths(nPoints, 0.0);

        // Calculate edge lengths and put into list
        forAll(pointPoints, neighPointI)
        {
            const label pointI = pointPoints[neighPointI];
            const double edgeLength = getPointDistance(mesh.points()[pointI], cCoords);
            edgeLengths[neighPointI] = edgeLength;
        }

        // Get order of sorted edgeLength
        // openfoam.com :
        // const labelList sLabels(Foam::sortedOrder(edgeLengths));
        //
        // openfoam.com and openfoam.org :
        labelList sLabels;
        Foam::sortedOrder(edgeLengths, sLabels);

        // Save local values
        closestPoints1[pointI] = mesh.points()[pointPoints[sLabels[0]]] - cCoords;
        closestPoints2[pointI] = mesh.points()[pointPoints[sLabels[1]]] - cCoords;
        closestPoints3[pointI] = mesh.points()[pointPoints[sLabels[2]]] - cCoords;

        // Check if closest two points share a cell
        if (findIndex(pointNeighPoints[pointPoints[sLabels[0]]], pointPoints[sLabels[1]]) >= 0)
            hasCommonCell[pointI] = true;
        else
            hasCommonCell[pointI] = false;
    }

    // For each position, deduce the globally closest points among processors

    // Position 1
    // ----------

    // Storage for sharing relative point locations among processors
    vectorField syncPoints(mesh.points().size(), Zero);

    forAll (mesh.points(), pointI)
    {
        syncPoints[pointI] = closestPoints1[pointI];
    }

    syncTools::syncPointList
    (
        mesh,
        syncPoints,
        minMagSqrEqOp<vector>(),
        UNDEF_VECTOR               // null value
    );

    forAll (mesh.points(), pointI)
    {
        if (isCloserPoint(syncPoints[pointI], closestPoints1[pointI]))
        {
            closestPoints3[pointI] = closestPoints2[pointI];
            closestPoints2[pointI] = closestPoints1[pointI];
            closestPoints1[pointI] = syncPoints[pointI];
            hasCommonCell[pointI] = false;
        }
    }

    // Position 2
    // ----------

    forAll (mesh.points(), pointI)
    {
        syncPoints[pointI] = closestPoints2[pointI];
    }

    syncTools::syncPointList
    (
        mesh,
        syncPoints,
        minMagSqrEqOp<vector>(),
        UNDEF_VECTOR               // null value
    );

    forAll (mesh.points(), pointI)
    {
        if (isCloserPoint(syncPoints[pointI], closestPoints2[pointI]))
        {
            closestPoints3[pointI] = closestPoints2[pointI];
            closestPoints2[pointI] = syncPoints[pointI];
            hasCommonCell[pointI] = false;
        }
    }

    // Position 3
    // ----------

    forAll (mesh.points(), pointI)
    {
        syncPoints[pointI] = closestPoints3[pointI];
    }

    syncTools::syncPointList
    (
        mesh,
        syncPoints,
        minMagSqrEqOp<vector>(),
        UNDEF_VECTOR               // null value
    );

    forAll (mesh.points(), pointI)
    {
        if (isCloserPoint(syncPoints[pointI], closestPoints3[pointI]))
        {
            closestPoints3[pointI] = syncPoints[pointI];
        }
    }

    // Synchronize cell sharing boolean
    syncTools::syncPointList
    (
         mesh,
         hasCommonCell,
         orEqOp<bool>(),
         false                      // null value
    );

    return 0;
}

// Help function for aspectRatioSmoothing to identify point locations
// of two short edge points which don't share cells.  Calculates also
// a blending fraction used for weighting the coordinates with
// centroidal smoothing coordinates. Blending fraction is set to zero
// if points share a cell, to prevent smoothing in that case.

double calcARSmoothingRatio
(
    const vector& closestPoint1,
    const vector& closestPoint2,
    const vector& closestPoint3,
    const bool hasCommonCell
)
{
    // Do nothing if the points are part of a common cell
    if (hasCommonCell)
    {
        return 0.0;
    }

    if ((closestPoint1 == ZERO_VECTOR) or (closestPoint2 == ZERO_VECTOR))
    {
        return 0.0;
    }

    // Minimum and maximum edge length ratios for detecting and
    // blending high aspect ratio.
    // TODO: Optimize values and test further?
    const double minRatio = 1.5;
    const double maxRatio = 3.0;

    // Blending fraction is applied if two closest points have similar
    // length, and if the third closest point is clearly farther away.
    const double lengthRatio1 = mag(closestPoint2) / mag(closestPoint1);
    const double lengthRatio2 = mag(closestPoint3) / mag(closestPoint2);

    if ((lengthRatio1 < minRatio) and (lengthRatio2 > minRatio))
    {
        const double frac = (lengthRatio2 - minRatio) / (maxRatio - minRatio);
        const double blendFrac = min(1.0, max(0.0, frac));
        return blendFrac;
    }

    return 0.0;
}

// Special smoothing algorithm for high aspect ratio cells. Blends
// resulting coordinates with centroidal point coordinates.

Foam::tmp<Foam::pointField> aspectRatioSmoothing
(
    const fvMesh& mesh,
    const boolList& isInternalPoint,
    const pointField& centroidalPoints,
    const labelListList& pointNeighPoints
)
{
    // Storage for surrounding point locations (relative to current
    // point) of closest edge points
    const label nPoints = mesh.nPoints();
    vectorList closestPoints1(nPoints, Zero);
    vectorList closestPoints2(nPoints, Zero);
    vectorList closestPoints3(nPoints, Zero);

    // Boolean to mark that the two shortest local edge points share a cell
    List<bool> hasCommonCell(nPoints, false);

    // New point locations storage
    tmp<pointField> tNewPoints(new pointField(nPoints, Zero));
    auto& newPoints = tNewPoints.ref();

    // Initialize with centroidal coordinates
    forAll(centroidalPoints, pointI)
    {
        newPoints[pointI] = centroidalPoints[pointI];
    }

    // Find closest points among processors
    findClosestPoints(mesh, isInternalPoint, closestPoints1, closestPoints2, closestPoints3, hasCommonCell, pointNeighPoints);

    // Calculate closest middle point coordinates and blend with centroidal coordinates
    forAll(mesh.points(), pointI)
    {
        // Process only internal points
        if (! isInternalPoint[pointI])
            continue;

        const double blendFrac = calcARSmoothingRatio(closestPoints1[pointI], closestPoints2[pointI], closestPoints3[pointI], hasCommonCell[pointI]);

        if (blendFrac > 0.0)
        {
            const vector aCoords = mesh.points()[pointI] + (closestPoints1[pointI] + closestPoints2[pointI]) / 2.0;
            const vector newCoords = (1.0 - blendFrac) * centroidalPoints[pointI] + blendFrac * aCoords;
            newPoints[pointI] = newCoords;
        }
    }

    return tNewPoints;
}


// Quality control function which prohibits the decrease of edge
// length by freezing points to current locations, if edge length is
// below minEdgeLength and length is decreasing. This feature is used
// to limit the squishing or compression of cells near concave
// features.

int restrictEdgeShortening
(
    const fvMesh& mesh,
    pointField& origPoints,
    const double minEdgeLength,
    const bool totalMinFreeze,
    boolList& isFrozenPoint
)
{
    forAll(origPoints, pointI)
    {
        if (isFrozenPoint[pointI])
            continue;

        const vector cCoords = mesh.points()[pointI];
        const vector nCoords = origPoints[pointI];

        // Calculate shortest edge length from current mesh point and
        // new mesh point
        double shortestCurrentEdgeLength = GREAT;
        double shortestNewEdgeLength = GREAT;
        forAll(mesh.pointPoints(pointI), pointPpI)
        {
            const label neighI = mesh.pointPoints(pointI)[pointPpI];
            const double testCurrentLength = getPointDistance(mesh.points()[neighI], cCoords);
            if (testCurrentLength < shortestCurrentEdgeLength)
                shortestCurrentEdgeLength = testCurrentLength;
            const double testNewLength = getPointDistance(mesh.points()[neighI], nCoords);
            if (testNewLength < shortestNewEdgeLength)
                shortestNewEdgeLength = testNewLength;
        }

        // Unconditional freeze below minEdgeLength (totalMinFreeze option)
        const double shortestLength = min(shortestNewEdgeLength, shortestCurrentEdgeLength);
        if (totalMinFreeze and (shortestLength < minEdgeLength))
        {
            isFrozenPoint[pointI] = true;
        }

        // Otherwise freeze only if edge length decreases and length
        // is below threshold value
        else if ((shortestNewEdgeLength < minEdgeLength) and
                 (shortestNewEdgeLength < shortestCurrentEdgeLength))
        {
            isFrozenPoint[pointI] = true;
        }

    }

    return 0;
}


// Calculate and return the maximum step length of the proposed new points

double getProposedMaxStepLength
(
    const fvMesh& mesh,
    const pointField& newPoints
)
{
    double maxLength = 0.0;

    forAll(mesh.points(), pointI)
    {
        const double length = mag(vector(newPoints[pointI] - mesh.points()[pointI]));

        if (length > maxLength)
            maxLength = length;
    }

    const double proposedMaxLength = returnReduce(maxLength, maxOp<double>());

    return proposedMaxLength;
}


// Quality control function to constrain the length of a step jump to
// new coordinates by an absolute length value. This increases the
// stability of the smoothing process, in case target coordinates are
// far off.

int constrainMaxStepLength
(
    const fvMesh& mesh,
    pointField& origPoints,
    const double maxStepLength,
    const double relStepFrac,
    const bool doGlobalScaling
)
{
    // Copy original points for temporary working point field
    tmp<pointField> tNewPoints(new pointField(mesh.nPoints(), Zero));
    auto& newPoints = tNewPoints.ref();

    forAll(newPoints, pointI)
        newPoints[pointI] = origPoints[pointI];

    // Global scaling factor is a means to scale down all point
    // movements with a common factor. It is meant to allow centroidal
    // smoothing to accelerate movement of points which need larger
    // step sizes than surrounding points for stable smoothing. If
    // global scaling factor is 1.0, then movement is scaled down by
    // the relative step fraction only.

    // Calculate global scaling factor from maximum proposed step size
    // and maximum allowed step size.
    const double proposedMaxLength = getProposedMaxStepLength(mesh, origPoints);

    double globalScale = 1.0;
    if (doGlobalScaling)
    {
        globalScale = min(1.0, maxStepLength / (proposedMaxLength * relStepFrac));
    }

    // Info << "Proposed maximum step length = " << proposedMaxLength << endl
    //      << "Maximum allowed step length = " << maxStepLength << endl
    //      << "Global scaling factor = " << globalScale << endl
    //      << endl << endl;

    forAll(origPoints, pointI)
    {
        // Scale down the length of the jump from current coordinates
        // towards new coordinates if jump would be too long
        const vector cCoords = mesh.points()[pointI];
        const vector stepDir = origPoints[pointI] - cCoords;

        // Max step length constraining if global scaling is disabled
        if (! doGlobalScaling)
        {
            if (mag(stepDir) > maxStepLength)
            {
                globalScale = maxStepLength / (mag(stepDir) * relStepFrac);
            }
            else
            {
                globalScale = 1.0;
            }
        }

        // Save the constrained point
        const vector nCoords = cCoords + relStepFrac * globalScale * stepDir;
        newPoints[pointI] = nCoords;
    }

    // Save new point coordinates to original point field
    forAll(origPoints, pointI)
    {
        origPoints[pointI] = newPoints[pointI];
    }

    return 0;
}


/////////////////////////////////////////////////////
// Help functions for restrictMinEdgeAngleDecrease //
/////////////////////////////////////////////////////

// Calculate and return the edge-edge angle (in radians,
// 0 < angle < pi) of two edges which share a common point
// at coordinate cCoords. The two edge end point coordinates are
// p1Coords and p2Coords.

double edgeEdgeAngle
(
    const vector cCoords,
    const vector p1Coords,
    const vector p2Coords
)
{
    vector vec1 = (p1Coords - cCoords);
    vector vec2 = (p2Coords - cCoords);
    vec1 /= mag(vec1);
    vec2 /= mag(vec2);

    const double cosA = vec1 & vec2;

    // Ensure cos angle is in sane range before calling arc cos
    const double MAX = 0.99999;
    const double cosAlpha = std::max(-MAX, std::min(MAX, cosA));
    const double angle = std::acos(cosAlpha);

    return angle;
}


// Find the neighbour mesh point indices neighPI1 and neighPI2 for the
// point with index pointI, where all points are part of face index
// faceI.

int getNeighbourPoints
(
    const fvMesh& mesh,
    const label centerPointI,
    const label faceI,
    label *neighPI1,
    label *neighPI2
)
{
    int i = 0;
    const face facePoints = mesh.faces()[faceI];
    const label nPoints = facePoints.size();

    forAll(facePoints, pointI)
    {
        if (facePoints[pointI] == centerPointI)
        {
            // Handle index wrapping at beginning and at end
            label prevI = pointI - 1;
            if (pointI == 0)
                prevI = nPoints - 1;

            label nextI = pointI + 1;
            if (pointI == nPoints - 1)
                nextI = 0;

            *neighPI1 = facePoints[prevI];
            *neighPI2 = facePoints[nextI];
            return 0;
        }
        i++;
    }

    FatalError << "Sanity broken, didn't find neighbour points for center point "
               << centerPointI << " and face " << faceI << endl
               << abort(FatalError);

    return 0;
}

// Calculate and store the minimum edge angles of point with index
// pointI for current mesh (minCAngleStorage) and a minimum for new
// mesh points (minNAngleStorage).

int calc_min_edge_angles
(
    const fvMesh& mesh,
    const pointField& newPoints,
    const label pointI,
    double *minCAngleStorage,
    double *minNAngleStorage
)
{
    // Current and new minimum angles
    double minCAngle = DBL_MAX;
    double minNAngle = DBL_MAX;

    forAll(mesh.pointFaces()[pointI], pointFI)
    {
        // For both edges connected at this point, find the other point index
        label neighPI1 = 0;
        label neighPI2 = 0;
        const label faceI = mesh.pointFaces()[pointI][pointFI];
        getNeighbourPoints(mesh, pointI, faceI, &neighPI1, &neighPI2);

        // Current angle from current mesh
        const vector cp0 = mesh.points()[pointI];
        const vector cp1 = mesh.points()[neighPI1];
        const vector cp2 = mesh.points()[neighPI2];
        const double cAngle = edgeEdgeAngle(cp0, cp1, cp2);

        // New angle from current mesh
        const vector np0 = newPoints[pointI];
        const double nAngle0 = edgeEdgeAngle(np0, cp1, cp2);

        // New angle from new mesh points
        const vector np1 = newPoints[neighPI1];
        const vector np2 = newPoints[neighPI2];
        const double nAngle1 = edgeEdgeAngle(np0, np1, np2);

        // New angle from combination of current and new mesh points
        // Note: These might not be really necessary, but do it for now
        // for completeness sake.
        const double nAngle2 = edgeEdgeAngle(np0, cp1, np2);
        const double nAngle3 = edgeEdgeAngle(np0, np1, cp2);

        // Use smallest of new angles.
        const double nAngle = min(min(min(nAngle0, nAngle1), nAngle2), nAngle3);

        // Update minimum values
        if (cAngle < minCAngle)
            minCAngle = cAngle;
        if (nAngle < minNAngle)
            minNAngle = nAngle;
    }

    // Save minimum angles to storage
    *minCAngleStorage = minCAngle;
    *minNAngleStorage = minNAngle;

    return 0;
}

// Restrict decrease of smallest edge-edge angles when angle is below
// minAngle (in degrees). This is meant to avoid creation of
// self-intersections for concave features.

int restrictMinEdgeAngleDecrease
(
    const fvMesh& mesh,
    pointField& origPoints,
    const double minEdgeAngleInDegrees,
    boolList& isFrozenPoint
)
{
    forAll(origPoints, pointI)
    {
        if (isFrozenPoint[pointI])
            continue;

        // Calculate current minimum angle and new minimum angle
        double minCAngle;
        double minNAngle;
        calc_min_edge_angles(mesh, origPoints, pointI, &minCAngle, &minNAngle);

        // If minimum angle is below threshold and would decrease,
        // then freeze the point. Note: This allows angle to
        // increase, so points are not frozen permanently.
        const double smallAngle = M_PI * minEdgeAngleInDegrees / 180.0;

        if ((minNAngle < smallAngle) and (minNAngle < minCAngle))
        {
            isFrozenPoint[pointI] = true;
        }
    }

    return 0;
}

///////////////////////////////////////////////////////
// Help functions for restrictFaceAngleDeterioration //
///////////////////////////////////////////////////////

// Map face angle information from edges to points

int mapCurrentMinMaxFaceAnglesToPoints
(
    const fvMesh& mesh,
    const List<double>& minFaceAnglesForEdges,
    const List<double>& maxFaceAnglesForEdges,
    List<double>& minFaceAnglesForPoints,
    List<double>& maxFaceAnglesForPoints
)
{
    forAll(mesh.points(), pointI)
    {
        minFaceAnglesForPoints[pointI] = 2.0 * M_PI;
        maxFaceAnglesForPoints[pointI] = 0.0;
    }

    forAll(mesh.edges(), edgeI)
    {
        const edge e = mesh.edges()[edgeI];

        {
            const label pointI = e.start();
            if (minFaceAnglesForPoints[pointI] > minFaceAnglesForEdges[edgeI])
                minFaceAnglesForPoints[pointI] = minFaceAnglesForEdges[edgeI];
            if (maxFaceAnglesForPoints[pointI] < maxFaceAnglesForEdges[edgeI])
                maxFaceAnglesForPoints[pointI] = maxFaceAnglesForEdges[edgeI];
        }

        {
            const label pointI = e.end();
            if (minFaceAnglesForPoints[pointI] > minFaceAnglesForEdges[edgeI])
                minFaceAnglesForPoints[pointI] = minFaceAnglesForEdges[edgeI];
            if (maxFaceAnglesForPoints[pointI] < maxFaceAnglesForEdges[edgeI])
                maxFaceAnglesForPoints[pointI] = maxFaceAnglesForEdges[edgeI];
        }
    }

    return 0;
}

// Calculate sum of angles between two edge unit vectors and a center
// unit vector. Assumes unit length for all vectors.

double calcEdgeCenterEdgeAngle
(
    const vector p0,
    const vector cC,
    const vector p1
)
{
    const double cosA0 = p0 & cC;
    const double cosA1 = cC & p1;

    // Ensure cos angle is in sane range before calling arc cos
    static const double MAX = 0.99999;
    const double cosAlpha0 = std::max(-MAX, std::min(MAX, cosA0));
    const double angle0 = std::acos(cosAlpha0);
    const double cosAlpha1 = std::max(-MAX, std::min(MAX, cosA1));
    const double angle1 = std::acos(cosAlpha1);

    return angle0 + angle1;
}

// Calculate minimum and maximum face angles from the collection of
// precalculated point data projected on the edge normal plane

int calcMinMaxFinalProjectedAngle
(
    const label nCells,
    const vectorList& pVecs,
    const vectorList& cVecs,
    const labelList& faceIs,
    const labelList& f0Is,
    const labelList& f1Is,
    double& minFaceAngle,
    double& maxFaceAngle
)
{
    double minAngle = 2.0 * M_PI;
    double maxAngle = 0.0;

    for (label i = 0; i < nCells; ++i)
    {
        const label f0I = f0Is[i];
        const label f1I = f1Is[i];
        const vector p0 = pVecs[f0I];
        const vector p1 = pVecs[f1I];
        const vector cC = cVecs[i];
        const double angle = calcEdgeCenterEdgeAngle(p0, cC, p1);

        if (angle < minAngle)
            minAngle = angle;
        if (angle > maxAngle)
            maxAngle = angle;
    }

    minFaceAngle = minAngle;
    maxFaceAngle = maxAngle;

    return 0;
}

// Find the pair of faces which are part of the cell cellI, and store
// the face label pair.

int findCellFacePair
(
    const fvMesh& mesh,
    const label cellI,
    const labelList& faceIs,
    const label nFaces,
    label& f0I,
    label& f1I
)
{
    label face0I = -1;
    label face1I = -1;
    static const int nInternalFaces = mesh.owner().size();

    for (label i = 0; i < nFaces; ++i)
    {
        const label faceI = faceIs[i];
        const label ownerI = mesh.owner()[faceI];

        // Guard that neighbour index search to below nInternalFaces
        // to avoid "Invalid read of size 4" error from valgrind

        label neighI = -1;
        if (faceI < nInternalFaces)
            neighI = mesh.neighbour()[faceI];

        if (ownerI == cellI)
        {
            if (face0I == -1)
                face0I = i;
            else
                face1I = i;
        }

        // neighI is not defined for boundary faces, therefore it
        // needs an extra check

        if ((faceI < nInternalFaces) and (neighI == cellI))
        {
            if (face0I == -1)
                face0I = i;
            else
                face1I = i;
        }
    }

    // Sanity check
    if ((face0I == -1) or (face1I == -1) or (face0I == face1I))
    {
        Info << "faceIs";
        for (int i = 0; i < nFaces; ++i)
            Info << " " << faceIs[i];
        Info << " nInternalFaces " << nInternalFaces << endl;
        FatalError << "Sanity broken, didn't find face pairs for cell "
                   << cellI << ". Face indices: " << face0I << " " << face1I << endl
                   << abort(FatalError);
    }

    // Save the face pair
    f0I = face0I;
    f1I = face1I;

    return 0;
}

// Calculate face center for face faceI as a weighted average of point coordinates.
// If point label pointI1 > -1, then that point is assumed to be at coords1.
// If point label pointI2 > -1, then that point is assumed to be at coords2.

vector calcFaceCenter
(
    const fvMesh& mesh,
    const label faceI,
    const label pointI1,
    const vector coords1,
    const label pointI2,
    const vector coords2
)
{
    vector center = vector(0, 0, 0);

    forAll(mesh.faces()[faceI], facePointI)
    {
        const label pointI = mesh.faces()[faceI][facePointI];

        if ((pointI1 >= 0) and (pointI == pointI1))
            center += coords1;
        else if ((pointI2 >= 0) and (pointI == pointI2))
            center += coords2;
        else
            center += mesh.points()[pointI];
    }

    center /= double(mesh.faces()[faceI].size());

    return center;
}

// Calculate minimum and maximum face angles for a single edge, with
// optional move of pointI1 to coords1 and pointI2 to coords2.

int calcMinMaxFaceAngleForEdge
(
    const fvMesh& mesh,
    const label edgeI,
    double& minFaceAngle,
    double& maxFaceAngle,
    const int pointI1,
    const vector coords1,
    const int pointI2,
    const vector coords2
)
{
    // Find all faces of this edge
    const labelList edgeFaces = mesh.edgeFaces(edgeI);
    const label nFaces = edgeFaces.size();

    const edge e = mesh.edges()[edgeI];

    // Edge start point coordinates
    const label e0I = e.start();
    vector e0 = mesh.points()[e0I];
    if ((pointI1 >= 0) and (e0I == pointI1))
        e0 = coords1;
    else if ((pointI2 >= 0) and (e0I == pointI2))
        e0 = coords2;

    // Edge end point coordinates
    const label e1I = e.end();
    vector e1 = mesh.points()[e1I];
    if ((pointI1 >= 0) and (e1I == pointI1))
        e1 = coords1;
    else if ((pointI2 >= 0) and (e1I == pointI2))
        e1 = coords2;

    // Edge center coordinates
    const vector cCoords = 0.5 * (e0 + e1);

    // Edge normal vector
    const vector eVec = (e1 - e0) / mag(e1 - e0);

    // Note: Edge center and edge normal vector defines the plane
    // where points are projected to, prior to angle calculation.

    // Calculate the projected face center vectors
    vectorList pVecs(nFaces, UNDEF_VECTOR);  // projected coordinates
    labelList faceIs(nFaces, UNDEF_LABEL);  // face indices of edge faces

    forAll(edgeFaces, edgeFaceI)
    {
        // Get face center coordinates
        const label faceI = edgeFaces[edgeFaceI];
        const vector fCoords = calcFaceCenter(mesh, faceI, pointI1, coords1, pointI2, coords2);

        // Project face center to the plane
        const vector cf = cCoords - fCoords;
        const double dotProd = cf & eVec;
        const vector pCoords = fCoords + dotProd * eVec;

        // Save projected face center coordinate vector
        const vector cp = (pCoords - cCoords) / mag(pCoords - cCoords);
        pVecs[edgeFaceI] = cp;

        // Save face index for pairing below
        faceIs[edgeFaceI] = faceI;
    }

    // Face-face angles are to be calculated for face pairs
    // belonging to all cells of the edge. There are always
    // exactly two faces for each cell, connected at the
    // edge. Search face indices for the face pairs.

    const labelList edgeCells = mesh.edgeCells(edgeI);
    const int nCells = edgeCells.size();
    labelList f0Is(nCells, UNDEF_LABEL);  // index of first face of the cell
    labelList f1Is(nCells, UNDEF_LABEL);  // index of second face of the cell
    vectorList cVecs(nCells, UNDEF_VECTOR);  // projected cell center coordinates

    // Calculate, project and save projected cell center coordinates
    forAll(edgeCells, i)
    {
        findCellFacePair(mesh, edgeCells[i], faceIs, nFaces, f0Is[i], f1Is[i]);
        const label cellI = edgeCells[i];
        const vector cellCenter = mesh.C()[cellI];
        const vector cf = cCoords - cellCenter;
        const double dotProd = cf & eVec;
        const vector pCoords = cellCenter + dotProd * eVec;
        const vector cp = (pCoords - cCoords) / mag(pCoords - cCoords);
        cVecs[i] = cp;
    }

    // Finally calculate the minimum and maximum angles using the
    // projected point data and save to storage
    calcMinMaxFinalProjectedAngle(nCells, pVecs, cVecs, faceIs, f0Is, f1Is, minFaceAngle, maxFaceAngle);

    return 0;
}

// Wrapper for calcMinMaxEdgeAngle for current mesh points

int calcMinMaxFaceAngleForCurrentMeshEdge
(
    const fvMesh& mesh,
    const label edgeI,
    double& minFaceAngle,
    double& maxFaceAngle
)
{
    static const vector dummy = vector(0, 0, 0);
    calcMinMaxFaceAngleForEdge(mesh, edgeI, minFaceAngle, maxFaceAngle, -1, dummy, -1, dummy);
    return 0;
}

// Calculate and store the minimum face angles for all mesh
// edges

int calcCurrentMinMaxFaceAnglesForEdges
(
    const fvMesh& mesh,
    List<double>& minFaceAnglesForEdges,
    List<double>& maxFaceAnglesForEdges
)
{
    forAll(mesh.edges(), edgeI)
    {
        double minAngle;
        double maxAngle;
        calcMinMaxFaceAngleForCurrentMeshEdge(mesh, edgeI, minAngle, maxAngle);
        minFaceAnglesForEdges[edgeI] = minAngle;
        maxFaceAnglesForEdges[edgeI] = maxAngle;
    }

    return 0;
}

// Calculate and store minimum and maximum face angle for argument
// first point pointI1, with optional movement of pointI2 to given
// coordinates.

int calcMinMaxFaceAngleForPoint
(
    const fvMesh& mesh,
    const int pointI1,
    const vector coords1,
    const int pointI2,
    const vector coords2,
    double& minFaceAngle,
    double& maxFaceAngle
)
{
    // Initialize
    minFaceAngle = 2.0 * M_PI;
    maxFaceAngle = 0.0;

    // Loop over all point edges to find min and max angle

    forAll(mesh.pointEdges()[pointI1], pointEdgeI)
    {
        double minAngle;
        double maxAngle;
        const label edgeI = mesh.pointEdges()[pointI1][pointEdgeI];
        calcMinMaxFaceAngleForEdge(mesh, edgeI, minAngle, maxAngle, pointI1, coords1, pointI2, coords2);

        if (minFaceAngle > minAngle)
            minFaceAngle = minAngle;
        if (maxFaceAngle < maxAngle)
            maxFaceAngle = maxAngle;
    }

    return 0;
}


// Quality control function to restrict the decrease of smallest
// face-face angles when angle is below minAngle (in degrees). This is
// meant to avoid creation of self-intersections. The face-face angle
// can be thought of a measure which quantifies squishing or folding
// of cell faces. In addition to prohibiting the deterioration of
// face-face angles of each point, it is also necessary to freeze
// movement of neighbouring points if their movement would deteriorate
// the angles at the current point.

int restrictFaceAngleDeterioration
(
    const fvMesh& mesh,
    pointField& origPoints,
    const double minFaceAngleInDegrees,
    const double maxFaceAngleInDegrees,
    boolList& isFrozenPoint
)
{
    // Calculate minimum and maximum face angles for all edges from
    // current mesh
    static const size_t nEdges = size_t(mesh.nEdges());
    List<double> currentMinAnglesForEdges(nEdges, GREAT);
    List<double> currentMaxAnglesForEdges(nEdges, GREAT);
    calcCurrentMinMaxFaceAnglesForEdges(mesh, currentMinAnglesForEdges, currentMaxAnglesForEdges);

    // Map face angle information from edges to points
    static const size_t nPoints = size_t(mesh.nPoints());
    List<double> currentMinAnglesForPoints(nPoints, GREAT);
    List<double> currentMaxAnglesForPoints(nPoints, GREAT);
    mapCurrentMinMaxFaceAnglesToPoints(mesh, currentMinAnglesForEdges, currentMaxAnglesForEdges, currentMinAnglesForPoints, currentMaxAnglesForPoints);

    // Use a stack to walk through all points. If a point is frozen by
    // it's neighbour then the point must be processed again to allow
    // recursive neighbour freezing. Use of stack is needed, since
    // list of items to be processed is modified on the run.
    std::stack<label> pointStack;

    // Initialize stack with all points. Angle calculation must be
    // made for all points, also boundary points, in order to be able
    // to stop angle deterioration at boundary by neighbour point
    // movement
    forAll(origPoints, pointI)
        pointStack.push(pointI);

    while (! pointStack.empty())
    {
        // Get and remove a point label from stack.
        const label pointI = pointStack.top();
        pointStack.pop();

        // 1. Check nothing for points whose face angles are in good range

        const double smallAngle = M_PI * minFaceAngleInDegrees / 180.0;
        const double largeAngle = M_PI * maxFaceAngleInDegrees / 180.0;

        if ((currentMinAnglesForPoints[pointI] > smallAngle) and
            (currentMaxAnglesForPoints[pointI] < largeAngle))
            continue;

        const vector cCoords = mesh.points()[pointI];
        vector nCoords = origPoints[pointI];

        // If this point is already frozen, set coordinates to current
        // mesh coordinates
        if (isFrozenPoint[pointI])
            nCoords = cCoords;

        // 2. Calculate new face angles for this point, if it's not
        // frozen. This point is hypothetically moved to new
        // coordinates, while surrounding points are kept at current
        // locations. Freeze this point if angle change with the move
        // of this point deteriorates the angles.

        if (nCoords != cCoords)
        {
            double newMinFaceAngle;
            double newMaxFaceAngle;
            calcMinMaxFaceAngleForPoint(mesh, pointI, nCoords, -1, nCoords, newMinFaceAngle, newMaxFaceAngle);

            if (((newMinFaceAngle < smallAngle) and
                (newMinFaceAngle < currentMinAnglesForPoints[pointI])) or
                ((newMaxFaceAngle > largeAngle) and
                (newMaxFaceAngle > currentMaxAnglesForPoints[pointI])))
            {
                // Freeze this point (self freeze)
                nCoords = cCoords;
                isFrozenPoint[pointI] = true;
            }
        }

        // 3. Calculate the effect from all neighbouring point movements
        // to face angles at this point. Freeze the neighbour point if
        // angle becomes worse than it currently is.

        forAll(mesh.pointPoints()[pointI], pointNI)
        {
            const label neighPointI = mesh.pointPoints()[pointI][pointNI];
            const vector neighCoords = origPoints[neighPointI];

            // Skip this neighbour point if it's not moving
            if (isFrozenPoint[neighPointI])
                continue;
            if (neighCoords == mesh.points()[neighPointI])
                continue;

            double newMinFaceAngle;
            double newMaxFaceAngle;
            calcMinMaxFaceAngleForPoint(mesh, pointI, nCoords, neighPointI, neighCoords, newMinFaceAngle, newMaxFaceAngle);

            if (((newMinFaceAngle < smallAngle) and
                (newMinFaceAngle < currentMinAnglesForPoints[pointI])) or
                ((newMaxFaceAngle > largeAngle) and
                (newMaxFaceAngle > currentMaxAnglesForPoints[pointI])))
            {
                // Freeze the neighbour point (neighbour freeze)
                isFrozenPoint[neighPointI] = true;

                // Add neighbour point index to stack list, as it
                // needs to be (re)checked after neighbour freezing
                pointStack.push(neighPointI);
            }
        }
    }

    return 0;
}


// Help function to get patch numbers from a given command line option

labelList getPatchIdsForOption
(
    const fvMesh& mesh,
    const Foam::argList& args,
    const word optionName
)
{
    labelList patchIds;

    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    labelHashSet patchesHashSet(patches.size());

    if (args.optionFound(optionName))
    {
        patchesHashSet = patches.patchSet
        (
            wordReList(args.optionLookup(optionName)())
        );

        forAllConstIter(labelHashSet, patchesHashSet, iter)
        {
            const label id = iter.key();
            if (findIndex(patchIds, id) == -1)
            {
                patchIds.append(id);
            }
        }
    }
    return patchIds;
}


// Calculate and return the minimum edge length of the current mesh

int getMeshMinMaxEdgeLength
(
    const fvMesh& mesh,
    double& minEdgeLength,
    double& maxEdgeLength
)
{
    double minLength = VGREAT;
    double maxLength = 0.0;

    forAll(mesh.edges(), edgeI)
    {
        const edge e = mesh.edges()[edgeI];

        const label startI = e.start();
        const vector startCoords = mesh.points()[startI];
        const label endI = e.end();
        const vector endCoords = mesh.points()[endI];

        const double length = mag(vector(endCoords - startCoords));

        if (length < minLength)
            minLength = length;
        if (length > maxLength)
            maxLength = length;
    }

    const double meshMinLength = returnReduce(minLength, minOp<double>());
    const double meshMaxLength = returnReduce(maxLength, maxOp<double>());

    minEdgeLength = meshMinLength;
    maxEdgeLength = meshMaxLength;
    return 0;
}

// Calculate the maximum relative change in the point positions, to be
// used as a measure of residual of smoothing

double calculateResidual
(
    const fvMesh& mesh,
    const pointField& newPoints,
    const boolList& isInternalPoint,
    const double maxStepLength
)
{
    scalar maxStep = ZERO;

    forAll(isInternalPoint, pointI)
    {
        const scalar distance =
            mag(newPoints[pointI] - mesh.points()[pointI]) / maxStepLength;

        if (distance > maxStep)
        {
            maxStep = distance;
        }
    }

    const scalar maxMaxStep = returnReduce(maxStep, maxOp<scalar>());

    return maxMaxStep;
}


// Help function to check does a file exists on the system

bool fileExists
(
    const std::string& fileName
)
{
    std::ifstream file(fileName);
    return file.good();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Move internal mesh points to increase mesh quality"
    );

    #include "addRegionOption.H"
    #include "addOverwriteOption.H"

    argList::addOption
    (
        "time",
        "time",
        "Specify the time (default is latest)"
    );


    argList::addOption
    (
        "centroidalIters",
        "label",
        "Maximum number of centroidal smoothing iterations (default 1000)"
    );

    argList::addOption
    (
        "maxStepLength",
        "double",
        "Maximum absolute step length applied in smoothing (default 0.01)"
    );

    argList::addOption
    (
        "edgeAngleConstraint",
        "bool",
        "Option to apply the minimum edge angle control constraint (default: true)"
    );

    argList::addOption
    (
        "faceAngleConstraint",
        "bool",
        "Option to apply the minimum and maximum face angle constraint (default: true)"
    );

    argList::addOption
    (
        "minEdgeLength",
        "double",
        "Edge length below which edge points are frozen if edge length would decrease in smoothing (default 0.05)"
    );

    argList::addOption
    (
        "totalMinFreeze",
        "bool",
        "Makes minEdgeLength an absolute requirement, freezing short edges even if edge length would increase in smoothing (default false)"
    );

    argList::addOption
    (
        "minAngle",
        "double",
        "Face-face angle below which points are frozen (in degrees, default: 35)"
    );

    argList::addOption
    (
        "maxAngle",
        "double",
        "Face-face angle above which points are frozen (in degrees, default: 160)"
    );

    argList::addOption
    (
        "layerMaxBlendingFraction",
        "double",
        "Maximum blending fraction to force prismatic boundary layer treatment on edges (default: 0)"
    );

    argList::addOption
    (
        "layerEdgeLength",
        "double",
        "Target thickness for first boundary layer (default: 0.05)"
    );

    argList::addOption
    (
        "layerExpansionRatio",
        "double",
        "The expansion ratio for the increase in boundary layer thickness (default: 1.3)"
    );

    argList::addOption
    (
        "minLayers",
        "label",
        "Number of outermost boundary layers that receive maximum blending (default: 1)"
    );

    argList::addOption
    (
        "maxLayers",
        "label",
        "Number of boundary layers affected by the boundary layer treatment (default: 4)"
    );

    argList::addOption
    (
        "layerPatches",
        "wordRe",
        "Specify single patch or multiple patches for the boundary layer treatment."
        " No patches are included by default."
        " For example 'walls' or '( stator \"rotor.*\" )'"
    );

    argList::addOption
    (
        "smoothingPatches",
        "wordRe",
        "Specify single patch or multiple patches for boundary point smoothing."
        " All patches are included by default."
        " For example 'walls' or '( stator \"rotor.*\" )'"
    );

    argList::addOption
    (
        "relTol",
        "double",
        "Relative tolerance for stopping the smoothing iterations (default: 0.02)"
    );

    argList::addOption
    (
        "writeInterval",
        "label",
        "Interval to write mesh during iterations (default value 0)"
    );

    #include "setRootCase.H"
    #include "createTime.H"

    const bool overwrite = args.optionFound("overwrite");

    // Handle time
    if (args.optionFound("time"))
    {
        if (args["time"] == "constant")
        {
            runTime.setTime(instant(0, "constant"), 0);
        }
        else
        {
            const scalar timeValue = args.optionRead<scalar>("time");
            runTime.setTime(instant(timeValue), 0);
        }
    }

    #ifdef OPENFOAM_ORG
        #include "createMesh.H"
    #else
        #include "createMeshNoClear.H"
    #endif

    const word oldInstance = mesh.pointsInstance();

    // Get patch ids for boundary layer treatment
    labelList layerPatchIds = getPatchIdsForOption(mesh, args, "layerPatches");
    if (layerPatchIds.size() > 0)
    {
        Info<< "Patches for boundary layer treatment: "
            << args["layerPatches"] << endl;
    }
    else
    {
        Info<< "Patches for boundary layer treatment: none"
            << endl;
    }

    // Add all patches as default for smoothing patches if nothing is
    // specified
    if (! args.optionFound("smoothingPatches"))
    {
        args.setOption("smoothingPatches", "(\".*\")");
    }

    labelList smoothingPatchIds = getPatchIdsForOption(mesh, args, "smoothingPatches");
    if (smoothingPatchIds.size() > 0)
    {
        Info<< "Patches for boundary point smoothing: "
            << args["smoothingPatches"] << endl;
    }
    else
    {
        Info<< "Patches for boundary point smoothing: none"
            << endl;
    }

    // Default minimum edge length is smaller than the minimum edge
    // length from initial mesh to allow good boundary smoothing
    double meshMinEdgeLength;
    double meshMaxEdgeLength;
    getMeshMinMaxEdgeLength(mesh, meshMinEdgeLength, meshMaxEdgeLength);

    double minEdgeLength =
        args.optionLookupOrDefault("minEdgeLength", 0.5 * meshMinEdgeLength);

    double maxStepLength =
        args.optionLookupOrDefault("maxStepLength", 0.3 * minEdgeLength);

    if (maxStepLength > 0.5 * minEdgeLength)
    {
        Info << "WARNING: The maximum allowed step length is more "
             << "than half of the minimum edge length! This may "
             << "cause unstability in smoothing." << endl << endl;
    }

    double relStepFrac =
        args.optionLookupOrDefault("relStepFrac", 0.5);

    bool totalMinFreeze =
        args.optionLookupOrDefault("totalMinFreeze", false);

    double minAngle =
        args.optionLookupOrDefault("minAngle", 35.0);

    double maxAngle =
        args.optionLookupOrDefault("maxAngle", 160.0);

    bool edgeAngleConstraint =
        args.optionLookupOrDefault("edgeAngleConstraint", true);

    bool faceAngleConstraint =
        args.optionLookupOrDefault("faceAngleConstraint", true);

    double layerMaxBlendingFraction =
        args.optionLookupOrDefault("layerMaxBlendingFraction", 0.5);

    double layerEdgeLength =
        args.optionLookupOrDefault("layerEdgeLength", minEdgeLength);

    double layerExpansionRatio =
        args.optionLookupOrDefault("layerExpansionRatio", 1.3);

    label minLayers =
        args.optionLookupOrDefault("minLayers", 1);

    label maxLayers =
        args.optionLookupOrDefault("maxLayers", 4);

    double relTol =
        args.optionLookupOrDefault("relTol", 0.02);

    label centroidalIters =
        args.optionLookupOrDefault("centroidalIters", 1000);

    label writeInterval =
        args.optionLookupOrDefault("writeInterval", 0);

    // Boundary point smoothing edge and surface meshes
    const string initEdgesFileString("constant/geometry/initEdges.obj");
    const string targetEdgesFileString("constant/geometry/targetEdges.obj");
    const string targetSurfacesFileString("constant/geometry/targetSurfaces.obj");

    const fileName initEdgesFileName(initEdgesFileString);
    const fileName targetEdgesFileName(targetEdgesFileString);
    const fileName targetSurfacesFileName(targetSurfacesFileString);

    // Print out applied parameter values
    Info << "Applying following parameter values in smoothing:" << endl;
    Info << "    centroidalIters        " << centroidalIters << endl;
    Info << "    relTol                 " << relTol << endl;
    Info << "    minEdgeLength          " << minEdgeLength << endl;
    Info << "    maxStepLength          " << maxStepLength << endl;
    Info << "    relStepFrac            " << relStepFrac << endl;
    Info << "    totalMinFreeze         " << totalMinFreeze << endl;

    if (edgeAngleConstraint)
    {
        Info << "    edgeAngleConstraint    true" << endl;
        Info << "    minAngle               " << minAngle << endl;
    }
    else
    {
        Info << "    edgeAngleConstraint    false (edge min angle quality constraint is NOT applied)" << endl;
    }

    if (faceAngleConstraint)
    {
        Info << "    faceAngleConstraint    true" << endl;
        Info << "    minAngle               " << minAngle << endl;
        Info << "    maxAngle               " << maxAngle << endl;
    }
    else
    {
        Info << "    faceAngleConstraint    false (face angle quality constraints are NOT applied)" << endl;
    }

    if (layerMaxBlendingFraction > SMALL)
    {
        Info << "    layerMaxBlendingFraction " << layerMaxBlendingFraction << endl;
        Info << "    layerEdgeLength          " << layerEdgeLength << endl;
        Info << "    layerExpansionRatio      " << layerExpansionRatio << endl;
        Info << "    minLayers                " << minLayers << endl;
        Info << "    maxLayers                " << maxLayers << endl;
    }
    else
    {
        Info << "    layerMaxBlendingFraction 0 (boundary layer treatment is NOT applied)" << endl;
    }

    Info << endl;


    // Storage for markers for internal points
    boolList isInternalPoint(mesh.nPoints(), false);
    const label nInternalPoints = findInternalMeshPoints(mesh, isInternalPoint);

    // Storage for number of edge hops to reach layer and free
    // boundaries for mesh points (for boundary layer treatment)
    labelList pointHopsToLayerBoundary(mesh.nPoints(), UNDEF_LABEL);
    labelList pointHopsToSmoothingBoundary(mesh.nPoints(), UNDEF_LABEL);

    // Storage for point normals (for boundary layer treatment)
    tmp<pointField> tPointNormals(new pointField(mesh.nPoints(), Zero));
    pointField& pointNormals = tPointNormals.ref();

    // Storage for neighbour point locations (for boundary layer treatment)
    tmp<pointField> tInnerNeighCoords(new pointField(mesh.nPoints(), UNDEF_VECTOR));
    pointField& innerNeighCoords = tInnerNeighCoords.ref();
    tmp<pointField> tOuterNeighCoords(new pointField(mesh.nPoints(), UNDEF_VECTOR));
    pointField& outerNeighCoords = tOuterNeighCoords.ref();

    // Storage for marking neighbour point being inside same processor
    // domain (for boundary layer treatment)
    boolList isInnerNeighInProc(mesh.nPoints(), false);
    boolList isOuterNeighInProc(mesh.nPoints(), false);

    // Storage for index map from point to neighbour point inside same
    // processor domain. One map points towards inner mesh, and
    // another map towards outer boundary. (for boundary layer treatment)
    labelList pointToInnerPointMap(mesh.nPoints(), UNDEF_LABEL);
    labelList pointToOuterPointMap(mesh.nPoints(), UNDEF_LABEL);

    // Boolean list for marking frozen points (points not allowed to
    // move during smoothing). This list is synced among processors.
    boolList isFrozenPoint(mesh.nPoints(), false);

    // A list of point label lists to indicate cell sharing. Used in
    // aspectRatioSmoothing.
    Info << "Starting to build pointNeighPoints (this may take some time)" << endl;
    labelListList pointNeighPoints(mesh.nPoints());
    generatePointNeighPoints(mesh, pointNeighPoints);
    Info << "Done building pointNeighPoints" << endl << endl;

    // Check prerequisites for carrying out boundary layer treatment
    bool doLayerTreatment = false;
    if ((layerPatchIds.size() > 0) and (layerMaxBlendingFraction > SMALL))
    {
        doLayerTreatment = true;
        Info << "Enabled boundary layer treatment" << endl << endl;
    }
    else
    {
        Info << "Boundary layer treatment is disabled. Either no layerPatches were specified or boundaryMaxBlendingFraction is zero" << endl << endl;
    }

    // Check prerequisites for carrying out boundary point smoothing
    bool doBoundarySmoothing = false;
    if ((fileExists(targetSurfacesFileString)) and
        (fileExists(initEdgesFileString)) and
        (smoothingPatchIds.size() > 0))
    {
        doBoundarySmoothing = true;
        Info << "Enabled boundary point smoothing" << endl << endl;
    }
    else
    {
        Info << "Boundary point smoothing is disabled. Missing smoothingPatches, or one or both of files:" << endl
             << targetSurfacesFileString << endl
             << initEdgesFileString << endl << endl;
    }

    if ((doLayerTreatment) and (! doBoundarySmoothing))
    {
        Info << "WARNING: Boundary layer treatment will be done without boundary point smoothing. This can result in distorted boundary cells." << endl << endl;
    }

    // Objects for boundary point snapping to surfaces
    autoPtr<triSurface> surf(nullptr);
    autoPtr<triSurfaceSearch> searchSurfaces(nullptr);
    autoPtr<indexedOctree<treeDataTriSurface>> tree(nullptr);

    // Objects for boundary point snapping to feature edges
    autoPtr<edgeMesh> initEdges(nullptr);
    autoPtr<edgeMesh> targetEdges(nullptr);

    // Point classification lists and corner target locations
    boolList isConnectedToInternalPoint(mesh.nPoints(), false);
    boolList isFeatureEdgePoint(mesh.nPoints(), false);
    boolList isLayerSurfacePoint(mesh.nPoints(), false);
    boolList isSmoothingSurfacePoint(mesh.nPoints(), false);
    boolList isFrozenSurfacePoint(mesh.nPoints(), false);
    boolList isProcessorPoint(mesh.nPoints(), false);
    boolList isCornerPoint(mesh.nPoints(), false);
    vectorList cornerPoints(mesh.nPoints(), UNDEF_VECTOR);

    // Closest edge mesh point indices (for feature edge points)
    labelList closestEdgePointIs(mesh.nPoints(), UNDEF_LABEL);

    // Preparations for boundary point smoothing
    if (doBoundarySmoothing)
    {
        // Target surface mesh, build search tree
        surf.reset(new triSurface(targetSurfacesFileName));
        Info << "Target surfaces file " << targetSurfacesFileName << " stats:" << endl;
        surf().writeStats(Info);
        Info << endl;

        searchSurfaces.reset(new triSurfaceSearch(surf));
        tree.reset(new indexedOctree<treeDataTriSurface>(searchSurfaces().tree()));

        // Initial and target feature edge meshes
        initEdges.reset(new edgeMesh(initEdgesFileName));
        Info << "Initial feature edges file " << initEdgesFileName << " stats:" << endl;
        initEdges().writeStats(Info);
        Info << endl;

        if (fileExists(targetEdgesFileString))
        {
            targetEdges.reset(new edgeMesh(targetEdgesFileName));
            Info << "Target feature edges file " << targetEdgesFileName << " stats:" << endl;
            targetEdges().writeStats(Info);
        }
        else
        {
            targetEdges.reset(new edgeMesh(initEdgesFileName));
            Info << "Warning: Initial feature edges will be used also as target edges, because" << endl

            << "did not find file " << targetEdgesFileString << "." << endl;
        }
        Info << endl;
    }
    else
    {
        // Initialize empty edge mesh, to avoid issues in
        // classifyBoundaryPoints()
        initEdges.reset(new edgeMesh());
        targetEdges.reset(new edgeMesh());
    }

    const label nPoints = returnReduce(mesh.nPoints(), sumOp<label>());
    Info << "Mesh includes a total of " << nPoints << " points:" << endl
         << "  - " << nInternalPoints << " internal (non-boundary) points" << endl
         << "  - " << nPoints - nInternalPoints << " boundary points" << endl
         << "Mesh minimum edge length = " << meshMinEdgeLength << endl
         << "Mesh maximum edge length = " << meshMaxEdgeLength << endl << endl;

    // Classify boundary points and find target corner points for feature edge snapping
    classifyBoundaryPoints
    (
        mesh,
        initEdges,
        targetEdges,
        layerPatchIds,
        smoothingPatchIds,
        isInternalPoint,
        isProcessorPoint,
        isConnectedToInternalPoint,
        isFeatureEdgePoint,
        isCornerPoint,
        cornerPoints,
        isLayerSurfacePoint,
        isSmoothingSurfacePoint,
        isFrozenSurfacePoint
    );

    // Preparations for optional smoothing and treatment
    if ((doBoundarySmoothing) or (doLayerTreatment))
    {
        calculatePointHopsToBoundary(mesh, layerPatchIds, isInternalPoint, isConnectedToInternalPoint, pointHopsToLayerBoundary, maxLayers + 1);
        calculatePointHopsToBoundary(mesh, smoothingPatchIds, isInternalPoint, isConnectedToInternalPoint, pointHopsToSmoothingBoundary, 2);
        calculateBoundaryPointNormals(mesh, pointNormals);
        propagateOuterNeighInfo(mesh, isInternalPoint, isLayerSurfacePoint, isOuterNeighInProc, pointToOuterPointMap, pointNormals, pointHopsToLayerBoundary, maxLayers + 1);
        propagateInnerNeighInfo(mesh, isSmoothingSurfacePoint, isConnectedToInternalPoint, isInnerNeighInProc, pointToInnerPointMap, pointHopsToSmoothingBoundary);
    }

    // Find initial closest edge points (for feature edge snapping)
    if (doBoundarySmoothing)
    {
        forAll(mesh.points(), pointI)
        {
            if (isFeatureEdgePoint[pointI])
            {
                const label newClosestEdgePointI = findClosestEdgeMeshPointIndex(mesh.points()[pointI], targetEdges, false, true);
                closestEdgePointIs[pointI] = newClosestEdgePointI;
            }
        }
    }

    // Carry out smoothing iterations
    // ------------------------------

    for (label i = 0; i < centroidalIters; ++i)
    {
        Info << "Starting iteration " << i << endl;

        // Reset frozen points
        forAll(isFrozenPoint, pointI)
            isFrozenPoint[pointI] = false;

        // Recalculate point normals
        calculateBoundaryPointNormals(mesh, pointNormals);

        // Calculate new point locations using centroidal smoothing
        tmp<pointField> tCentroidalPoints = centroidalSmoothing(mesh, isInternalPoint, doBoundarySmoothing);
        pointField& centroidalPoints = tCentroidalPoints.ref();

        // Constrain absolute length of jump to new coordinates, to stabilize smoothing
        constrainMaxStepLength(mesh, centroidalPoints, maxStepLength, relStepFrac, false);

        // Blend centroidal points with points from aspect ratio smoothing
        tmp<pointField> tNewPoints = aspectRatioSmoothing(mesh, isInternalPoint, centroidalPoints, pointNeighPoints);
        pointField& newPoints = tNewPoints.ref();

        // Optional boundary layer treatment
        if (doLayerTreatment)
        {
            // Update neighbour coordinates and synchronize among processors
            updateNeighCoords(mesh, isOuterNeighInProc, pointToOuterPointMap, outerNeighCoords);
            // Blend orthogonal and centroidal coordinates to newPoints
            blendWithOrthogonalPoints
            (
                 mesh,
                 newPoints,
                 isInternalPoint,
                 pointHopsToLayerBoundary,
                 pointNormals,
                 outerNeighCoords,
                 layerMaxBlendingFraction,
                 layerEdgeLength,
                 layerExpansionRatio,
                 minLayers,
                 maxLayers + 1  // +1 for correct number of layers
            );

            // Constrain absolute length of jump to new coordinates, to stabilize smoothing
            constrainMaxStepLength(mesh, newPoints, maxStepLength, relStepFrac, false);
        }

        if (doBoundarySmoothing)
        {
            // Update neighbour coordinates and synchronize among processors
            updateNeighCoords(mesh, isInnerNeighInProc, pointToInnerPointMap, innerNeighCoords);
            // Project boundary points
            projectBoundaryPointsToEdgesAndSurfaces
            (
                mesh,
                newPoints,
                pointNormals,
                isInternalPoint,
                isSmoothingSurfacePoint,
                isFeatureEdgePoint,
                isCornerPoint,
                cornerPoints,
                targetEdges,
                closestEdgePointIs,
                surf,
                tree,
                meshMaxEdgeLength
            );

            // Constrain absolute length of jump to new coordinates, to stabilize smoothing
            // constrainMaxStepLength(mesh, newPoints, maxStepLength, relStepFrac, false);

            // Use the locations of first cell layer points for
            // projecting points to boundary surfaces
            projectFreeBoundaryPointsToSurfaces
            (
                mesh,
                newPoints,
                pointHopsToSmoothingBoundary,
                pointNormals,
                isSmoothingSurfacePoint,
                isConnectedToInternalPoint,
                isFeatureEdgePoint,
                isCornerPoint,
                innerNeighCoords
            );

            // Constrain absolute length of jump to new coordinates, to stabilize smoothing
            constrainMaxStepLength(mesh, newPoints, maxStepLength, relStepFrac, false);
        }

        // Avoid shortening of short edge length
        restrictEdgeShortening(mesh, newPoints, minEdgeLength, totalMinFreeze, isFrozenPoint);

        if (edgeAngleConstraint)
        {
            // Restrict decrease of smallest edge-edge angle
            restrictMinEdgeAngleDecrease(mesh, newPoints, minAngle, isFrozenPoint);
        }

        if (faceAngleConstraint)
        {
            // Restrict deterioration of face-face angles
            restrictFaceAngleDeterioration(mesh, newPoints, minAngle, maxAngle, isFrozenPoint);
        }

        // Synchronize and combine the list of frozen points
        syncTools::syncPointList
        (
            mesh,
            isFrozenPoint,
            orEqOp<bool>(),
            false               // null value
        );

        // Restore original coordinates for frozen points
        label nFrozenPoints = 0;
        forAll(newPoints, pointI)
        {
            if (isFrozenPoint[pointI])
            {
                newPoints[pointI] = mesh.points()[pointI];
                ++nFrozenPoints;
            }
        }

        // Calculate and print residual
        const double res = calculateResidual(mesh, newPoints, isInternalPoint, maxStepLength);
        Info << "Smoothing iteration=" << (i + 1) << " nFrozenPoints=" << returnReduce(nFrozenPoints, sumOp<label>()) << " residual=" << res << endl;

        // Push coordinates to mesh
        mesh.movePoints(tNewPoints);

        if (res < relTol)
        {
            Info << "Residual reached relTol, stopping." << endl;
            break;
        }

        if (i == (centroidalIters - 1))
        {
            Info << "Maximum centroidalIters reached, stopping." << endl;
        }

        // Increase time
        runTime++;

        if ((writeInterval > 0) and ((i % writeInterval) == 0))
        {

            // Save mesh
            if (overwrite)
            {
                mesh.setInstance(oldInstance);
            }

            // Set the precision of the points data to 10
            IOstream::defaultPrecision(max(10u, IOstream::defaultPrecision()));

            Info << "Writing new mesh to time " << runTime.name()
                 << endl << endl;

            mesh.write();
        }
    }

    // Save mesh
    {
        if (overwrite)
        {
            mesh.setInstance(oldInstance);
        }

        // Set the precision of the points data to 10
        IOstream::defaultPrecision(max(10u, IOstream::defaultPrecision()));

        Info << "Writing new mesh to time " << runTime.name()
             << endl << endl;

        mesh.write();
    }

    Info<< "ClockTime = "
        << runTime.elapsedClockTime() << " s." << endl;

    Info << endl << "End" << endl;

    return 0;
}

// ************************************************************************* //
