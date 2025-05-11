/*---------------------------------------------------------------------------*\
Application
    smoothMesh

Description
    Smooth internal mesh points to improve mesh quality
\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "ReadFields.H"
#include "regionProperties.H"
#include "syncTools.H"
#include "weightedPosition.H"
#include "meshTools.H"
#include <bits/stdc++.h> // For std::stack
#include "vectorList.H"
#include "stringListOps.H"  // For stringListOps::findMatching()
#include "SortList.H"

// Macros for value definitions
#define UNDEF_LABEL -1
#define UNDEF_VECTOR vector(GREAT, GREAT, GREAT)
#define ZERO_VECTOR vector(0, 0, 0)

// Boolean for developer mode. Set to false to use bleeding edge
// work in progress features.
#define USE_STABLE_FEATURES_ONLY true

// #include <typeinfo>
// Typeinfo is needed only for getting types while debugging, for example:
// Info << "Type is " << typeid(x).name() << endl;
// Use terminal command like this to demangle the mangled name:
// c++filt -t N4Foam4faceE

#include "orthogonalBoundaryBlending.C"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Help function to find internal points of the argument mesh.
// Updates argument bitSet accordingly: true for internal points
// (including processor points) and false for boundary points.

int findInternalMeshPoints
(
    const fvMesh& mesh,
    bitSet& isInternalPoint
)
{
    // Start from all points in
    forAll(isInternalPoint, pointI)
        isInternalPoint.set(pointI);

    // Remove points on boundary patches from bit set, except not
    // processor patches nor empty patches
    const faceList& faces = mesh.faces();
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
            forAll (f, pointI)
            {
                const label i = faces[faceI][pointI];
                isInternalPoint.unset(i);
            }
        }
    }

    label nPoints = 0;
    forAll(isInternalPoint, pointI)
    {
        if (isInternalPoint.test(pointI))
            ++nPoints;
    }

    const label sumNPoints = returnReduce(nPoints, sumOp<label>());
    return sumNPoints;
}

// Function for centroidal smoothing of internal mesh points.
// Adapted from Foam::snappySnapDriver::smoothInternalDisplacement in
// https://develop.openfoam.com/Development/openfoam/-/blob/OpenFOAM-v2312/src/mesh/snappyHexMesh/snappyHexMeshDriver/snappySnapDriver.C

Foam::tmp<Foam::pointField> centroidalSmoothing
(
    const fvMesh& mesh,
    const bitSet isInternalPoint
)
{
    // Centroidal smoothing algorithm

    // Calculate number and sum of surrounding cell center
    // coordinates using weightedPosition class for data storage

    Field<weightedPosition> wps
    (
        mesh.nPoints(),
        pTraits<weightedPosition>::zero
    );

    forAll(isInternalPoint, pointI)
    {
        if (isInternalPoint.test(pointI))
        {
            const labelList& pCells = mesh.pointCells(pointI);

            // First element of Tuple2 stores number of entries
            wps[pointI].first() = pCells.size();

            // Second element of Tuple2 stores the sum of coordinates
            for (const label celli : pCells)
            {
                wps[pointI].second() += mesh.cellCentres()[celli];
            }
        }
    }

    // Synchronize among processors
    weightedPosition::syncPoints(mesh, wps);

    // Calculate new point locations
    tmp<pointField> tnewPoints(new pointField(mesh.nPoints(), Zero));
    pointField& newPoints = tnewPoints.ref();

    label nPoints = 0;
    forAll(newPoints, pointI)
    {
        const weightedPosition& wp = wps[pointI];

        // internal point
        if (mag(wp.first()) > VSMALL)
        {
            newPoints[pointI] =
                wp.second()/wp.first();
            nPoints++;
        }

        // boundary point
        else
        {
            newPoints[pointI] = mesh.points()[pointI];
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

// Generate a set of neighbour point pairs from pointPoints where the
// point pair share a cell. Used in aspectRatioSmoothing.

int generatePointNeighSet
(
    const fvMesh& mesh,
    std::set<std::pair<label, label>>& pointNeighSet
)
{
    // TODO: Building the set is super slow, fix it

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
                const std::pair<label, label> pointPair = {pointI, pointPointI};
                if (pointNeighSet.count(pointPair) == 1)
                    continue;
                pointNeighSet.insert(pointPair);
            }
        }
    }

    return 0;
}

// Help function for aspectRatioSmoothing to find out closest
// three points (connected by edges) to each point

int findClosestProcPoints
(
    const fvMesh& mesh,
    const bitSet isInternalPoint,
    vectorList& closestProcPoints1,
    vectorList& closestProcPoints2,
    vectorList& closestProcPoints3,
    vectorList& closestSyncPoints1,
    vectorList& closestSyncPoints2,
    boolList& hasCommonCell,
    const std::set<std::pair<label, label>>& pointNeighSet
)
{
    forAll (mesh.points(), pointI)
    {
        if (! isInternalPoint.test(pointI))
        {
            closestProcPoints1[pointI] = ZERO_VECTOR;
            closestProcPoints2[pointI] = ZERO_VECTOR;
            closestProcPoints3[pointI] = ZERO_VECTOR;
            closestSyncPoints1[pointI] = ZERO_VECTOR;
            closestSyncPoints2[pointI] = ZERO_VECTOR;
            continue;
        }

        const vector cCoords = mesh.points()[pointI];
        const labelList pointPoints = mesh.pointPoints()[pointI];
        const label nPoints = pointPoints.size();
        List<scalar> edgeLengths(nPoints, 0.0);

        // Calculate edge lengths and put into list
        forAll(pointPoints, neighPointI)
        {
            const label pointI = pointPoints[neighPointI];
            const double edgeLength = getPointDistance(mesh.points()[pointI], cCoords);
            edgeLengths[neighPointI] = edgeLength;
        }

        // Get labels of sorted length list
        SortList<scalar> sortedEdgeLengths(edgeLengths);
        sortedEdgeLengths.sort();
        const labelList sLabels = sortedEdgeLengths.indices();

        // Save values
        closestProcPoints1[pointI] = mesh.points()[pointPoints[sLabels[0]]] - cCoords;
        closestProcPoints2[pointI] = mesh.points()[pointPoints[sLabels[1]]] - cCoords;
        closestProcPoints3[pointI] = mesh.points()[pointPoints[sLabels[2]]] - cCoords;

        closestSyncPoints1[pointI] = mesh.points()[pointPoints[sLabels[0]]] - cCoords;
        closestSyncPoints2[pointI] = mesh.points()[pointPoints[sLabels[1]]] - cCoords;

        // Check if closest two points share a cell
        const std::pair<label, label> pointPair = {pointPoints[sLabels[0]], pointPoints[sLabels[1]]};
        hasCommonCell[pointI] = (pointNeighSet.count(pointPair) == 1);

        // Info << pointI << " at " << mesh.points()[pointI] << " closest1 " << closestProcPoints1[pointI] << " closest2 " << closestProcPoints2[pointI] << " closest3 " << closestProcPoints3[pointI] << " sLabel0 " << sLabels[0] << " sLabel1 " << sLabels[1] << " pointPair " << pointPair << " pointNeighSet " << pointNeighSet.count(pointPair) << " hasCommonCell " << hasCommonCell[pointI] << endl;
    }

    return 0;
}

// Help function for aspectRatioSmoothing to identify point locations
// of two short edge points which don't share cells.  Calculates also
// a blending fraction used for weighting the coordinates with
// centroidal smoothing coordinates. Blending fraction is set to zero
// if points share a cell, to prevent smoothing in that case.

int deduceClosestPoints
(
    const fvMesh& mesh,
    const label pointI,
    const vector& closestProcPoint1,
    const vector& closestProcPoint2,
    const vector& closestProcPoint3,
    const vector& closestSyncPoint1,
    const vector& closestSyncPoint2,
    vector& closestPoint1,
    vector& closestPoint2,
    double& blendFrac,
    bool hasCommonCell
)
{
    // Minimum and maximum edge length ratios for detecting and
    // blending high aspect ratio.
    // TODO: Optimize values and test further?
    const double minRatio = 1.5;
    const double maxRatio = 3.0;

    // Initialize
    closestPoint1 = ZERO_VECTOR;
    closestPoint2 = ZERO_VECTOR;
    blendFrac = 0.0;

    const vectorList points = {closestProcPoint1, closestProcPoint2, closestProcPoint3, closestSyncPoint1, closestSyncPoint2};

    // Calculate edge lengths and put into list
    List<scalar> edgeLengths(points.size(), 0.0);
    forAll(points, neighPointI)
    {
        const double edgeLength = mag(points[neighPointI]);
        edgeLengths[neighPointI] = edgeLength;
    }

    // Get labels of uniquely sorted length list
    SortList<scalar> sortedEdgeLengths(edgeLengths);
    sortedEdgeLengths.uniqueSort();

    // Exit if there is no data for points
    if (sortedEdgeLengths.size() == 1)
    {
        return 0;
    }

    const labelList sLabels = sortedEdgeLengths.indices();

    // Info << "pointI " << pointI << " points " << points << endl;
    // Info << "  sortedEdgeLengths " << sortedEdgeLengths << endl;
    closestPoint1 = points[sLabels[0]];
    closestPoint2 = points[sLabels[1]];
    const vector closestPoint3 = points[sLabels[2]];
    // Info << "  closestPoints " << closestPoint1 << " " << closestPoint2 << " " << closestPoint3 << endl;

    // Make note if synchronization changed the closest point
    bool syncAffectedClosest = false;
    if ((closestProcPoint1 != closestPoint1) or
        (closestProcPoint2 != closestPoint2))
        {
            syncAffectedClosest = true;
        }

    // Info << "  hasCommonCell " << hasCommonCell << " syncAffectedClosest " << syncAffectedClosest << endl;

    // Blending fraction is kept at zero if this processor contains
    // both closest points and if the points are part of a common cell
    if ((! syncAffectedClosest) and (hasCommonCell))
    {
        return 0;
    }

    // Blending fraction is applied if two closest points have similar
    // length, and if the third closest point is clearly farther away.
    const double lengthRatio1 = mag(closestPoint2) / mag(closestPoint1);
    const double lengthRatio2 = mag(closestPoint3) / mag(closestPoint2);

    if ((lengthRatio1 < minRatio) and (lengthRatio2 > minRatio))
    {
        const double frac = (lengthRatio2 - minRatio) / (maxRatio - minRatio);
        blendFrac = min(1.0, max(0.0, frac));
    }

    // Info << "  lengthRatio1 " << lengthRatio1 << " lengthRatio2 " << lengthRatio2 << " blendFrac " << blendFrac << endl;

    return 0;
}

// Special smoothing algorithm for high aspect ratio cells. Blends
// resulting coordinates with centroidal point coordinates.

Foam::tmp<Foam::pointField> aspectRatioSmoothing
(
    const fvMesh& mesh,
    const bitSet isInternalPoint,
    const pointField& centroidalPoints,
    const std::set<std::pair<label, label>>& pointNeighSet
)
{
    // Storage for surrounding point locations (relative to current
    // point) of closest edge points
    const label nPoints = mesh.nPoints();
    vectorList closestProcPoints1(nPoints, ZERO_VECTOR);
    vectorList closestProcPoints2(nPoints, ZERO_VECTOR);
    vectorList closestProcPoints3(nPoints, ZERO_VECTOR);

    // Storage for sharing closest edge relative point locations
    // among processors
    tmp<vectorField> tClosestSyncPoints1(new vectorField(nPoints, ZERO_VECTOR));
    vectorField& closestSyncPoints1 = tClosestSyncPoints1.ref();
    tmp<vectorField> tClosestSyncPoints2(new vectorField(nPoints, ZERO_VECTOR));
    vectorField& closestSyncPoints2 = tClosestSyncPoints2.ref();

    // Boolean to mark that the two shortest local edge points share a cell
    boolList hasCommonCell(nPoints, false);

    // New point locations storage
    tmp<pointField> tNewPoints(new pointField(nPoints, Zero));
    pointField& newPoints = tNewPoints.ref();

    // Initialize with centroidal coordinates
    forAll(isInternalPoint, pointI)
    {
        newPoints[pointI] = centroidalPoints[pointI];
    }

    if (USE_STABLE_FEATURES_ONLY)
        return tNewPoints;

    // Find closest local edge points and initialize closest
    // point info for synchronization among processors
    findClosestProcPoints(mesh, isInternalPoint, closestProcPoints1, closestProcPoints2, closestProcPoints3, closestSyncPoints1, closestSyncPoints2, hasCommonCell, pointNeighSet);

    // Synchronize closest points
    syncTools::syncPointList
    (
         mesh,
         closestSyncPoints1,
         minMagSqrEqOp<vector>(),
         UNDEF_VECTOR               // null value
    );

    syncTools::syncPointList
    (
         mesh,
         closestSyncPoints2,
         minMagSqrEqOp<vector>(),
         UNDEF_VECTOR               // null value
    );

    // Deduce the closest point coordinates and blending fraction,
    // then calculate new point coordinates
    forAll(isInternalPoint, pointI)
    {
        double blendFrac = 0.0;
        vector closestPoint1(ZERO_VECTOR);
        vector closestPoint2(ZERO_VECTOR);

        deduceClosestPoints(mesh, pointI, closestProcPoints1[pointI], closestProcPoints2[pointI], closestProcPoints3[pointI], closestSyncPoints1[pointI], closestSyncPoints2[pointI], closestPoint1, closestPoint2, blendFrac, hasCommonCell[pointI]);

        if (blendFrac > 0.0)
        {
            const vector aCoords = mesh.points()[pointI] + (closestPoint1 + closestPoint2) / 2.0;
            const vector newCoords = (1.0 - blendFrac) * newPoints[pointI] + blendFrac * aCoords;
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

// Quality control function to constrain the length of a step jump to
// new coordinates by an absolute length value. This increases the
// stability of the smoothing process, in case target coordinates are
// far off.

int constrainMaxStepLength
(
    const fvMesh& mesh,
    pointField& origPoints,
    const double maxStepLength
)
{
    // Copy original points for temporary working point field
    tmp<pointField> tNewPoints(new pointField(mesh.nPoints(), Zero));
    pointField& newPoints = tNewPoints.ref();
    forAll(newPoints, pointI)
        newPoints[pointI] = origPoints[pointI];

    forAll(origPoints, pointI)
    {
        const vector cCoords = mesh.points()[pointI];
        vector nCoords = origPoints[pointI];

        // Scale down the length of the jump from current coordinates
        // towards new coordinates if jump would be too long
        const vector stepDir = nCoords - cCoords;
        const double stepLength = mag(stepDir);
        if (stepLength > maxStepLength)
        {
            const double scale = maxStepLength / stepLength;
            nCoords = cCoords + scale * stepDir;
        }

        // Save the constrained point
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
    vec1.normalise();
    vec2.normalise();

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
    const vector eVec = (e1 - e0).normalise();

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
        const vector cp = (pCoords - cCoords).normalise();
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
        const vector cp = (pCoords - cCoords).normalise();
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


// Help function to get patch numbers from regex. Copied from
// https://develop.openfoam.com/Development/openfoam/-/blob/OpenFOAM-v2412/applications/utilities/surface/surfaceMeshExtract/surfaceMeshExtract.C

labelList getSelectedPatches
(
    const polyBoundaryMesh& patches,
    const wordRes& allow,
    const wordRes& deny
)
{
    // Name-based selection
    labelList indices
    (
        stringListOps::findMatching
        (
            patches,
            allow,
            deny,
            nameOp<polyPatch>()
        )
    );


    // Remove undesirable patches

    label count = 0;
    for (const label patchi : indices)
    {
        const polyPatch& pp = patches[patchi];

        if (isType<emptyPolyPatch>(pp))
        {
            continue;
        }
        else if (Pstream::parRun() && bool(isA<processorPolyPatch>(pp)))
        {
            break; // No processor patches for parallel output
        }

        indices[count] = patchi;
        ++count;
    }

    indices.resize(count);

    return indices;
}

// Calculate and return the minimum edge length of the current mesh

double getMeshMinEdgeLength
(
    const fvMesh& mesh
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

    Info << "Mesh minimum edge length = " << meshMinLength << endl;
    Info << "Mesh maximum edge length = " << meshMaxLength << endl << endl;

    return meshMinLength;
}

// Calculate the average relative change in the internal point
// positions, to be used as a measure of residual of the effect of
// smoothing

double calculateResidual
(
    const fvMesh& mesh,
    const pointField& newPoints,
    const bitSet isInternalPoint,
    const double maxStepLength
)
{
    label nPoints = 0;
    double res = 0.0;

    forAll(isInternalPoint, pointI)
    {
        if (! isInternalPoint.test(pointI))
            continue;

        ++nPoints;

        // Normalise the step length by the maximum step length
        res += mag(newPoints[pointI] - mesh.points()[pointI]) / maxStepLength;
    }

    const label sumNPoints = returnReduce(nPoints, sumOp<label>());
    const double sumRes = returnReduce(res, sumOp<double>());

    return (sumRes / sumNPoints);
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
        "boundaryMaxBlendingFraction",
        "double",
        "Maximum blending fraction to force prismatic boundary layer treatment on edges (default: 0)"
    );

    argList::addOption
    (
        "boundaryMaxPointBlendingFraction",
        "double",
        "Maximum blending fraction to move boundary points towards orthogonal position (default: 0)"
    );

    argList::addOption
    (
        "boundaryEdgeLength",
        "double",
        "Target thickness for first boundary layer (default: 0.05)"
    );

    argList::addOption
    (
        "boundaryExpansionRatio",
        "double",
        "The expansion ratio for the increase in boundary layer thickness (default: 1.3)"
    );

    argList::addOption
    (
        "boundaryMinLayers",
        "label",
        "Number of outermost boundary layers that receive maximum blending (default: 1)"
    );

    argList::addOption
    (
        "boundaryMaxLayers",
        "label",
        "Number of boundary layers affected by the boundary layer treatment (default: 4)"
    );

    argList::addOption
    (
        "patches",
        "wordRes",
        "Specify single patch or multiple patches for the boundary layer treatment.\n"
        "All patches are included by default.\n"
        "For example 'walls' or '( stator \"rotor.*\" )'"
    );

    argList::addOption
    (
        "boundaryPointSmoothingPatches",
        "wordRes",
        "Specify single patch or multiple patches for boundary point smoothing.\n"
        "All patches are included by default.\n"
        "For example 'walls' or '( stator \"rotor.*\" )'"
    );

    argList::addOption
    (
        "relTol",
        "double",
        "Relative tolerance for stopping the smoothing iterations (default: 0.02)"
    );

    argList::addOption
    (
        "boundaryPointSmoothing",
        "bool",
        "Boolean option to allow smoothing of boundary points (default false)"
    );


    #include "addOverwriteOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedMesh.H"

    const bool overwrite = args.found("overwrite");
    const word oldInstance = mesh.pointsInstance();

    // Handle time
    if (args.found("time"))
    {
        if (args["time"] == "constant")
        {
            runTime.setTime(instant(0, "constant"), 0);
        }
        else
        {
            const scalar timeValue = args.get<scalar>("time");
            runTime.setTime(instant(timeValue), 0);
        }
    }

    // Read in patches for boundary layer treatment
    wordRes includePatches, excludePatches;
    if (args.readListIfPresent<wordRe>("patches", includePatches))
    {
        Info<< "Patches for boundary layer treatment limited to: "
            << flatOutput(includePatches) << endl;
    }
    else
    {
        Info<< "Patches for boundary layer treatment limited to: all patches"
            << endl;
    }

    // Read in patches for boundary point smoothing
    wordRes includeBoundaryPointSmoothingPatches, excludeBoundaryPointSmoothingPatches;
    if (args.readListIfPresent<wordRe>("boundaryPointSmoothingPatches", includeBoundaryPointSmoothingPatches))
    {
        Info<< "Patches for boundary point smoothing limited to: "
            << flatOutput(includeBoundaryPointSmoothingPatches) << nl << endl;
    }
    else
    {
        Info<< "Patches for boundary point smoothing limited to: all patches"
            << nl << endl;
    }

    // Generate a list of patches for boundary layer treatment
    const polyBoundaryMesh& bMesh = mesh.boundaryMesh();
    const labelList patchIds = getSelectedPatches(bMesh, includePatches, excludePatches);
    const labelList boundaryPointSmoothingPatchIds = getSelectedPatches(bMesh, includeBoundaryPointSmoothingPatches, excludeBoundaryPointSmoothingPatches);

    // Default minimum edge length is smaller than the minimum edge
    // length from initial mesh to allow good boundary smoothing
    double minEdgeLength(0.5 * getMeshMinEdgeLength(mesh));
    args.readIfPresent("minEdgeLength", minEdgeLength);

    double maxStepLength(0.3 * minEdgeLength);
    args.readIfPresent("maxStepLength", maxStepLength);

    if (maxStepLength > 0.5 * minEdgeLength)
    {
        Info << "WARNING: The maximum allowed step length is more "
             << "than half of the minimum edge length! This may "
             << "cause unstability in smoothing." << endl << endl;
    }

    bool totalMinFreeze(false);
    args.readIfPresent("totalMinFreeze", totalMinFreeze);

    double minAngle(35);
    args.readIfPresent("minAngle", minAngle);

    double maxAngle(160);
    args.readIfPresent("maxAngle", maxAngle);

    bool edgeAngleConstraint(true);
    args.readIfPresent("edgeAngleConstraint", edgeAngleConstraint);

    bool faceAngleConstraint(true);
    args.readIfPresent("faceAngleConstraint", faceAngleConstraint);

    double boundaryMaxBlendingFraction(0.0);
    args.readIfPresent("boundaryMaxBlendingFraction", boundaryMaxBlendingFraction);

    double boundaryMaxPointBlendingFraction(0.0);
    args.readIfPresent("boundaryMaxPointBlendingFraction", boundaryMaxPointBlendingFraction);

    double boundaryEdgeLength(minEdgeLength);
    args.readIfPresent("boundaryEdgeLength", boundaryEdgeLength);

    double boundaryExpansionRatio(1.3);
    args.readIfPresent("boundaryExpansionRatio", boundaryExpansionRatio);

    label boundaryMinLayers(1);
    args.readIfPresent("boundaryMinLayers", boundaryMinLayers);

    label boundaryMaxLayers(4);
    args.readIfPresent("boundaryMaxLayers", boundaryMaxLayers);

    double relTol(0.02);
    args.readIfPresent("relTol", relTol);

    label centroidalIters(1000);
    args.readIfPresent("centroidalIters", centroidalIters);

    // Print out applied parameter values
    Info << "Applying following parameter values in smoothing:" << endl;
    Info << "    centroidalIters        " << centroidalIters << endl;
    Info << "    relTol                 " << relTol << endl;
    Info << "    minEdgeLength          " << minEdgeLength << endl;
    Info << "    maxStepLength          " << maxStepLength << endl;
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

    if ((boundaryMaxBlendingFraction > SMALL) or
        (boundaryMaxPointBlendingFraction > SMALL))
    {
        Info << "    boundaryMaxBlendingFraction " << boundaryMaxBlendingFraction << endl;
        Info << "    boundaryMaxPointBlendingFraction " << boundaryMaxPointBlendingFraction << endl;
        Info << "    boundaryEdgeLength     " << boundaryEdgeLength << endl;
        Info << "    boundaryExpansionRatio " << boundaryExpansionRatio << endl;
        Info << "    boundaryMinLayers      " << boundaryMinLayers << endl;
        Info << "    boundaryMaxLayers      " << boundaryMaxLayers << endl;
    }
    else
    {
        Info << "    boundaryMaxBlendingFraction 0 (boundary layer treatment is NOT applied)" << endl;
        Info << "    boundaryMaxPointBlendingFraction 0 (boundary point smoothing is NOT applied)" << endl;
    }

    Info << endl;


    // Storage for markers for internal points
    bitSet isInternalPoint(mesh.nPoints(), false);
    const label nPoints = findInternalMeshPoints(mesh, isInternalPoint);
    Info << "Starting smoothing of " << nPoints << " internal mesh points" << endl << endl;

    // Storage for number of edge hops to reach boundary for all mesh
    // points (for boundary layer treatment)
    labelList pointHopsToBoundary(mesh.nPoints(), -1);

    // Storage for point normals (for boundary layer treatment)
    tmp<pointField> tPointNormals(new pointField(mesh.nPoints(), ZERO_VECTOR));
    pointField& pointNormals = tPointNormals.ref();

    // Storage for neighbour point locations (for boundary layer treatment)
    tmp<pointField> tInnerNeighCoords(new pointField(mesh.nPoints(), UNDEF_VECTOR));
    pointField& innerNeighCoords = tInnerNeighCoords.ref();
    tmp<pointField> tOuterNeighCoords(new pointField(mesh.nPoints(), UNDEF_VECTOR));
    pointField& outerNeighCoords = tOuterNeighCoords.ref();

    // Storage for marking neighbour point being inside same processor
    // domain (for boundary layer treatment)
    bitSet isInnerNeighInProc(mesh.nPoints(), false);
    bitSet isOuterNeighInProc(mesh.nPoints(), false);

    // Storage for index map from point to neighbour point inside same
    // processor domain. One map points towards inner mesh, and
    // another map towards boundary. (for boundary layer treatment)
    labelList pointToInnerPointMap(mesh.nPoints(), UNDEF_LABEL);
    labelList pointToOuterPointMap(mesh.nPoints(), UNDEF_LABEL);

    // Storage for marking points on a flat (not curving) boundary patch
    boolList isFlatPatchPoint(mesh.nPoints(), false);

    // Boolean list for marking frozen points. This list is synced among processors.
    boolList isFrozenPoint(mesh.nPoints(), false);

    // A set for marking pointPoints which share a cell. Used in
    // aspectRatioSmoothing. Uses two point labels as a key.
    Info << "Starting to build pointNeighSet" << endl;
    std::set<std::pair<label, label>> pointNeighSet;
    generatePointNeighSet(mesh, pointNeighSet);
    Info << "Done building pointNeighSet" << endl;

    // Preparations for optional orthogonal boundary layer treatment
    if ((boundaryMaxBlendingFraction > SMALL) or
        (boundaryMaxPointBlendingFraction > SMALL))
    {
        calculatePointHopsToBoundary(mesh, pointHopsToBoundary, boundaryMaxLayers + 1);
        calculateBoundaryPointNormals(mesh, pointNormals, isFlatPatchPoint);
        propagateOuterNeighInfo(mesh, patchIds, isOuterNeighInProc, pointToOuterPointMap, pointNormals, pointHopsToBoundary, boundaryMaxLayers + 1);
        propagateInnerNeighInfo(mesh, boundaryPointSmoothingPatchIds, isInnerNeighInProc, pointToInnerPointMap, pointHopsToBoundary);
    }

    // Carry out smoothing iterations
    for (label i = 0; i < centroidalIters; ++i)
    {
        // Reset frozen points
        forAll(isFrozenPoint, pointI)
            isFrozenPoint[pointI] = false;

        // Calculate new point locations using centroidal smoothing
        tmp<pointField> tCentroidalPoints = centroidalSmoothing(mesh, isInternalPoint);
        pointField& centroidalPoints = tCentroidalPoints.ref();

        // Blend centroidal points with points from aspect ratio smoothing
        tmp<pointField> tNewPoints = aspectRatioSmoothing(mesh, isInternalPoint, centroidalPoints, pointNeighSet);
        pointField& newPoints = tNewPoints.ref();


        // Optional orthogonal boundary layer treatment
        if ((boundaryMaxBlendingFraction > SMALL) or
            (boundaryMaxPointBlendingFraction > SMALL))
        {
            // Update neighbour coordinates and synchronize among processors
            updateNeighCoords(mesh, isInnerNeighInProc, pointToInnerPointMap, innerNeighCoords);
            updateNeighCoords(mesh, isOuterNeighInProc, pointToOuterPointMap, outerNeighCoords);

            // Blend orthogonal and centroidal coordinates to newPoints
            blendWithOrthogonalPoints
            (
                 mesh,
                 newPoints,
                 isInternalPoint,
                 pointHopsToBoundary,
                 pointNormals,
                 outerNeighCoords,
                 boundaryMaxBlendingFraction,
                 boundaryEdgeLength,
                 boundaryExpansionRatio,
                 boundaryMinLayers,
                 boundaryMaxLayers + 1  // +1 for correct number of layers
            );

            projectBoundaryPoints
            (
                mesh,
                newPoints,
                isFlatPatchPoint,
                pointHopsToBoundary,
                pointNormals,
                innerNeighCoords,
                boundaryPointSmoothingPatchIds,
                boundaryMaxPointBlendingFraction
            );
        }

        // Constrain absolute length of jump to new coordinates, to stabilize smoothing
        constrainMaxStepLength(mesh, newPoints, maxStepLength);

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
        forAll(newPoints, pointI)
        {
            if (isFrozenPoint[pointI])
                newPoints[pointI] = mesh.points()[pointI];
        }

        // Calculate and print residual
        const double res = calculateResidual(mesh, newPoints, isInternalPoint, maxStepLength);
        Info << "Smoothing iteration=" << (i + 1) << " residual=" << res << endl;

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
    }

    // Save mesh
    {
        if (!overwrite)
        {
            ++runTime;
            mesh.setInstance(runTime.timeName());
        }
        else
        {
            mesh.setInstance(oldInstance);
        }

        // Set the precision of the points data to 10
        IOstream::defaultPrecision(max(10u, IOstream::defaultPrecision()));

        Info << "Writing new mesh to time " << runTime.timeName()
             << nl << endl;

        mesh.write();
        runTime.printExecutionTime(Info);
    }

    Info << nl << "End" << nl << endl;

    return 0;
}

// ************************************************************************* //
