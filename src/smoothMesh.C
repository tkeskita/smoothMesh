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
#include <float.h>
#include <bits/stdc++.h> // For std::stack

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
    // processor patches
    const faceList& faces = mesh.faces();
    forAll(mesh.boundary(), patchI)
    {
        const polyPatch& pp = mesh.boundaryMesh()[patchI];
        if ((! isA<processorPolyPatch>(pp)) and (! isA<emptyPolyPatch>(pp)))
        {
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
    }

    return 0;
}

// Function for centroidal smoothing of internal mesh points.
// Adapted from Foam::snappySnapDriver::smoothInternalDisplacement in
// https://develop.openfoam.com/Development/openfoam/-/blob/OpenFOAM-v2312/src/mesh/snappyHexMesh/snappyHexMeshDriver/snappySnapDriver.C

Foam::tmp<Foam::pointField> centroidalSmoothing
(
    const fvMesh& mesh,
    const label nIter,
    const bitSet isMovingPoint
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

    forAll(isMovingPoint, pointI)
    {
        if (isMovingPoint.test(pointI))
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

    Info << "Iteration " << nIter << ": centroidal smoothing of "
         << returnReduce(nPoints, sumOp<label>())
         << " points" << endl;

    return tnewPoints;
}

// Help function to return distance to current mesh point with index
// pointI from argument coords.

double getPointDistance
(
    const fvMesh& mesh,
    const label pointI,
    const vector coords
)
{
    const vector v = mesh.points()[pointI] - coords;
    return mag(v);
}

// Prohibits decrease of edge length by freezing points to current
// locations, if edge length is below minEdgeLength and length is
// decreasing. This feature is used to limit the squishing or
// compression of cells near concave features.

int restrictEdgeShortening
(
    const fvMesh& mesh,
    pointField& origPoints,
    const bitSet isMovingPoint,
    const double minEdgeLength,
    const bool totalMinFreeze,
    boolList& isFrozenPoint
)
{
    forAll(origPoints, pointI)
    {
        if (! isMovingPoint.test(pointI))
            continue;
        if (isFrozenPoint[pointI])
            continue;

        const vector cCoords = mesh.points()[pointI];
        const vector nCoords = origPoints[pointI];

        // Calculate shortest edge length from current mesh point and
        // new mesh point
        double shortestCurrentEdgeLength = DBL_MAX;
        double shortestNewEdgeLength = DBL_MAX;
        forAll(mesh.pointPoints(pointI), pointPpI)
        {
            const label neighI = mesh.pointPoints(pointI)[pointPpI];
            const double testCurrentLength = getPointDistance(mesh, neighI, cCoords);
            if (testCurrentLength < shortestCurrentEdgeLength)
                shortestCurrentEdgeLength = testCurrentLength;
            const double testNewLength = getPointDistance(mesh, neighI, nCoords);
            if (testNewLength < shortestNewEdgeLength)
                shortestNewEdgeLength = testNewLength;
        }

        // Unconditional freeze below minEdgeLength (totalMinFreeze option)
        const double shortestLength = min(shortestNewEdgeLength, shortestCurrentEdgeLength);
        if (totalMinFreeze and (shortestLength < minEdgeLength))
        {
            isFrozenPoint[pointI] = true;
        }

        // Freezeif edge length decreases and length is below threshold value
        else if ((shortestNewEdgeLength < minEdgeLength) and
                 (shortestNewEdgeLength < shortestCurrentEdgeLength))
        {
            isFrozenPoint[pointI] = true;
        }

    }

    return 0;
}

// Constrain the length of a step jump to new coordinates by an
// absolute length value. This increases the stability of the
// smoothing process, in case target coordinates are far off.

int constrainMaxStepLength
(
    const fvMesh& mesh,
    pointField& origPoints,
    const bitSet isMovingPoint,
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
        if (! isMovingPoint.test(pointI))
            continue;

        const vector cCoords = mesh.points()[pointI];
        vector nCoords = origPoints[pointI];

        // Scale down the length of the jump from current coordinates
        // towards new coordinates if jump would be otherwise too long
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


///////////////////////////////////////////////////////
// Help functions for restrictFaceAngleDeterioration //
///////////////////////////////////////////////////////

// Map face angle information from edges to points

int mapCurrentMinMaxFaceAnglesToPoints
(
    const fvMesh& mesh,
    const List<double> &minFaceAnglesForEdges,
    const List<double> &maxFaceAnglesForEdges,
    List<double> &minFaceAnglesForPoints,
    List<double> &maxFaceAnglesForPoints
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

// Calculate angle between two edge points using a midpoint.
// Assumes unit length for all point vectors.

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
    const vector *pVecs,
    const vector *cVecs,
    const label *faceIs,
    const label *f0Is,
    const label *f1Is,
    double *minFaceAngle,
    double *maxFaceAngle
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
        // Info << "i " << i << " f0I=" << f0I << " f1I=" << f1I << " -- " << p0 << p1 << cC << " angle " << angle << endl;

        if (angle < minAngle)
            minAngle = angle;
        if (angle > maxAngle)
            maxAngle = angle;
    }

    *minFaceAngle = minAngle;
    *maxFaceAngle = maxAngle;
    return 0;
}

// Finds the pair of faces which are part of the cellI, and stores the
// face labels.

int findCellFacePair
(
    const fvMesh& mesh,
    const label cellI,
    const label *faceIs,
    const label nFaces,
    label *f0I,
    label *f1I
)
{
    label face0I = -1;
    label face1I = -1;
    static const int nInternalFaces = mesh.owner().size();

    for (int i = 0; i < nFaces; ++i)
    {
        const label faceI = faceIs[i];
        const label ownerI = mesh.owner()[faceI];
        const label neighI = mesh.neighbour()[faceI];
        // Info << "faceI " << faceI << " owner " << ownerI << " neigh " << neighI << endl;

        if (ownerI == cellI)
        {
            if (face0I == -1)
                face0I = i;
            else
                face1I = i;
        }

        // neighI is zero for boundary faces, needs extra check
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
    *f0I = face0I;
    *f1I = face1I;
    // Info << "cellI " << cellI << " face pair: " << face0I << " " << face1I << endl;

    return 0;
}

// Calculate face center for face face as a weighted average of point coordinates.
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
// optional move of pointI1 and pointI2 to given coordinates.

int calcMinMaxFaceAngleForEdge
(
    const fvMesh &mesh,
    const label edgeI,
    double *minFaceAngle,
    double *maxFaceAngle,
    const int pointI1,
    const vector coords1,
    const int pointI2,
    const vector coords2
)
{
    // Find all faces of this edge
    const List<label> edgeFaces = mesh.edgeFaces(edgeI);
    const label nFaces = edgeFaces.size();
    // Info << "edge " << edgeI << " faces " << nFaces << endl;

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

    // Edge center and edge normal vector defines the plane where
    // points are projected to, prior to angle calculation.

    // Calculate the projected face center vectors
    vector pVecs[nFaces];  // projected coordinates
    label faceIs[nFaces];  // face indices of edge faces

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
    const List<label> edgeCells = mesh.edgeCells(edgeI);
    const int nCells = edgeCells.size();
    label f0Is[nCells];  // index of first face of the cell
    label f1Is[nCells];  // index of second face of the cell
    vector cVecs[nCells];  // projected cell center

    forAll(edgeCells, i)
    {
        findCellFacePair(mesh, edgeCells[i], faceIs, nFaces, &f0Is[i], &f1Is[i]);

        // Calculate, project and save projected cell center coordinates
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

    // Info << "     edgeI " << edgeI << " minFaceAngle " << *minFaceAngle << " maxFaceAngle " << *maxFaceAngle << endl;
    return 0;
}

// Wrapper for calcMinMaxEdgeAngle for current mesh points

int calcMinMaxFaceAngleForCurrentMeshEdge
(
    const fvMesh &mesh,
    const label edgeI,
    double *minFaceAngle,
    double *maxFaceAngle
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
    List<double> &minFaceAnglesForEdges,
    List<double> &maxFaceAnglesForEdges
)
{
    forAll(mesh.edges(), edgeI)
    {
        double minAngle;
        double maxAngle;
        calcMinMaxFaceAngleForCurrentMeshEdge(mesh, edgeI, &minAngle, &maxAngle);
        minFaceAnglesForEdges[edgeI] = minAngle;
        maxFaceAnglesForEdges[edgeI] = maxAngle;
    }

    return 0;
}

// Calculate and store minimum and maximum face angle for argument
// first point, with optional movement of points to given coordinates.

int calcMinMaxFaceAngleForPoint
(
    const fvMesh &mesh,
    const int pointI1,
    const vector coords1,
    const int pointI2,
    const vector coords2,
    double *minFaceAngle,
    double *maxFaceAngle
)
{
    // Initialize
    *minFaceAngle = 2.0 * M_PI;
    *maxFaceAngle = 0.0;

    // Loop over all point edges to find min and max angle
    forAll(mesh.pointEdges()[pointI1], pointEdgeI)
    {
        double minAngle;
        double maxAngle;
        const label edgeI = mesh.pointEdges()[pointI1][pointEdgeI];
        calcMinMaxFaceAngleForEdge(mesh, edgeI, &minAngle, &maxAngle, pointI1, coords1, pointI2, coords2);

        // Info << "    -- edgeI " << edgeI << " minAngle " << minAngle << " maxAngle " << maxAngle << endl;

        if (*minFaceAngle > minAngle)
            *minFaceAngle = minAngle;
        if (*maxFaceAngle < maxAngle)
            *maxFaceAngle = maxAngle;
    }

    return 0;
}


// Restrict decrease of smallest face-face angles when angle is below
// minAngle (in degrees). This is another angle heuristic (in addition
// to edge-edge angle heuristic), meant to avoid creation of
// self-intersections for concave features. The edge-edge angle
// version is not enough for cases where cell is flat and warped, and
// edge directions are highly non-orthogonal, e.g. due to edge length
// decrease. Face-face angle can be thought of a measure which is
// correlated with squished or folded cell shapes. In addition to
// prohibit the decrease of minimum face-face angle of a point, it is
// also necessary to freeze movement of neighbouring points if their
// movement would decrease the minimum angle at the current point.

int restrictFaceAngleDeterioration
(
    const fvMesh& mesh,
    pointField& origPoints,
    const bitSet isMovingPoint,
    const double minFaceAngleInDegrees,
    const double maxFaceAngleInDegrees,
    boolList& isFrozenPoint
)
{
    // Calculate minimum and maximum face angles for all edges from
    // current mesh
    static const size_t nEdges = size_t(mesh.nEdges());
    List<double> currentMinAnglesForEdges(nEdges);
    List<double> currentMaxAnglesForEdges(nEdges);
    calcCurrentMinMaxFaceAnglesForEdges(mesh, currentMinAnglesForEdges, currentMaxAnglesForEdges);

    // Map face angle information from edges to points
    static const size_t nPoints = size_t(mesh.nPoints());
    List<double> currentMinAnglesForPoints(nPoints);
    List<double> currentMaxAnglesForPoints(nPoints);
    mapCurrentMinMaxFaceAnglesToPoints(mesh, currentMinAnglesForEdges, currentMaxAnglesForEdges, currentMinAnglesForPoints, currentMaxAnglesForPoints);

    // Debug
    // labelList test = mesh.cellPoints()[7702];
    // Info << "Cell 7702 cellCenter " << mesh.C()[7702] << endl;
    // forAll (test, i)
    // {
    //     const label pointI = test[i];
    //     Info << "pointI " << pointI
    //          << " min angle " << currentMinAnglesForPoints[pointI]
    //          << " max angle " << currentMaxAnglesForPoints[pointI]
    //          << endl;
    // }

    // Use a stack to walk through all points. If a point is frozen by
    // it's neighbour then the point must be processed again to allow
    // recursive neighbour freezing. Use of stack is needed, since
    // list of items to be processed is modified on the run.
    std::stack<label> pointStack;

    // Initialize stack with all points. Angle calculation must be
    // made for all points, also boundary points, in order to be able
    // to stop angle deterioration at boundary by neighbour point
    // movement
    forAll(origPoints, pointi)
        pointStack.push(pointi);

    while (! pointStack.empty())
    {
        // Get and remove a point label from stack.
        const label pointI = pointStack.top();
        pointStack.pop();

        // Info << "===== Processing point " << pointI << " with currentMinAngle " << currentMinAnglesForPoints[pointI] << " and currentMaxAngle " << currentMaxAnglesForPoints[pointI] <<  endl;

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
            calcMinMaxFaceAngleForPoint(mesh, pointI, nCoords, -1, nCoords, &newMinFaceAngle, &newMaxFaceAngle);
            // Info << "   - own newMinFaceAngle " << newMinFaceAngle << " newMaxFaceAngle " << newMaxFaceAngle << endl;

            // Debug
            // if (pointI == 6252)
            // {
            //     Info << " MinFaceAngle c&n:" << currentMinAnglesForPoints[pointI] << " " << newMinFaceAngle << " MaxFaceAngle c&n:" << currentMaxAnglesForPoints[pointI] << " " << newMaxFaceAngle << endl;
            // }

            if (((newMinFaceAngle < smallAngle) and
                (newMinFaceAngle < currentMinAnglesForPoints[pointI])) or
                ((newMaxFaceAngle > largeAngle) and
                (newMaxFaceAngle > currentMaxAnglesForPoints[pointI])))
            {
                // Freeze this point (self freeze)
                nCoords = cCoords;
                isFrozenPoint[pointI] = true;
                // Info << "-- Self-froze point " << pointI << " MinFaceAngle c&n:" << currentMinAnglesForPoints[pointI] << " " << newMinFaceAngle << " MaxFaceAngle c&n:" << currentMaxAnglesForPoints[pointI] << " " << newMaxFaceAngle << endl;
            }
        }

        // Debug
        // if (pointI == 6252)
        // {
        //     Info << "pointI 6252 cCoords " << cCoords
        //          << " nCoords " << nCoords << endl;
        // }

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
            // Info << "   - checking moving neighbour point " << neighPointI << endl;

            double newMinFaceAngle;
            double newMaxFaceAngle;
            calcMinMaxFaceAngleForPoint(mesh, pointI, nCoords, neighPointI, neighCoords, &newMinFaceAngle, &newMaxFaceAngle);
            // Info << "   - neighbour " << neighPointI << " newMinFaceAngle " << newMinFaceAngle << " newMaxFaceAngle " << newMaxFaceAngle << endl;
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
                // Info << "-- point " << pointI << " froze neighbour " << neighPointI << " MinFaceAngle c&n:" << currentMinAnglesForPoints[pointI] << " " << newMinFaceAngle << " MaxFaceAngle c&n:" << currentMaxAnglesForPoints[pointI] << " " << newMaxFaceAngle << endl;
            }

            // Debug
            // if (pointI == 6252)
            // {
            //     Info << "pointI 6252 neighPointI " << neighPointI << " oldMaxFaceAngle " << currentMaxAnglesForPoints[pointI];
            //     Info << " newMaxFaceAngle " << newMaxFaceAngle << endl;
            // }
        }
    }

    return 0;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Move mesh points to increase mesh quality"
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
        "Number of centroidal smoothing iterations (default 20)"
    );

    argList::addOption
    (
        "maxStepLength",
        "double",
        "Maximum absolute step length applied in smoothing (default 0.01)"
    );

    argList::addOption
    (
        "orthogonalBlendingFraction",
        "double",
        "Fraction to force orthogonal side edges on the boundary (default 0 (=no orthogonality is forced))"
    );

    argList::addOption
    (
        "faceAngleConstraint",
        "bool",
        "Option to apply the minimum and maximum face  angle control constraint (default: true)"
    );

    argList::addOption
    (
        "minEdgeLength",
        "double",
        "A quality control constraint: Edge length below which edge vertices are fully frozen, but only if edge length would decrease in smoothing (default 0.05)"
    );

    argList::addOption
    (
        "totalMinFreeze",
        "bool",
        "A quality control constraint: Make minEdgeLength an absolute requirement, freezing short edges even if edge length would increase in smoothing (default false)"
    );

    argList::addOption
    (
        "minAngle",
        "double",
        "A quality control constraint: Face-face angle below which vertices are fully frozen (in degrees, default: 35)"
    );

    argList::addOption
    (
        "maxAngle",
        "double",
        "A quality control constraint: Face-face angle above which vertices are fully frozen (in degrees, default: 170)"
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

    double orthogonalBlendingFraction(0.0);
    args.readIfPresent("orthogonalBlendingFraction", orthogonalBlendingFraction);

    double maxStepLength(0.01);
    args.readIfPresent("maxStepLength", maxStepLength);

    double minEdgeLength(0.05);
    args.readIfPresent("minEdgeLength", minEdgeLength);

    if (maxStepLength > 0.5 * minEdgeLength)
    {
        Pout << "WARNING: The maximum allowed step length is more "
             << "than half of the minimum edge length! This may "
             << "cause unstability in smoothing." << endl << endl;
    }

    bool totalMinFreeze(false);
    args.readIfPresent("totalMinFreeze", totalMinFreeze);

    double minAngle(35);
    args.readIfPresent("minAngle", minAngle);

    double maxAngle(170);
    args.readIfPresent("maxAngle", maxAngle);

    bool faceAngleConstraint(true);
    args.readIfPresent("faceAngleConstraint", faceAngleConstraint);

    // Storage for markers for internal points
    bitSet isInternalPoint(mesh.nPoints());

    // Storage for point normals (for boundary smoothing)
    tmp<pointField> tPointNormals(new pointField(mesh.nPoints(), Zero));
    pointField& pointNormals = tPointNormals.ref();

    // Storage for markers for existence of  point normals (for boundary smoothing)
    bitSet hasPointNormals(mesh.nPoints(), false);

    // Storage for point-to-boundary-point map (for boundary smoothing)
    labelList uniValenceBoundaryMap(mesh.nPoints());

    findInternalMeshPoints(mesh, isInternalPoint);
    calculatePointNormals(mesh, pointNormals, hasPointNormals);
    calculateUniValenceBoundaryMap(mesh, uniValenceBoundaryMap, isInternalPoint);

    // Boolean list for marking frozen points. This list is synced among processors.
    boolList isFrozenPoint(mesh.nPoints(), false);

    label centroidalIters(20);
    args.readIfPresent("centroidalIters", centroidalIters);

    // Carry out centroidal smoothing iterations
    for (label i = 0; i < centroidalIters; ++i)
    {
        // Reset frozen points
        forAll(isFrozenPoint, pointI)
            isFrozenPoint[pointI] = false;

        tmp<pointField> tNewPoints = centroidalSmoothing(mesh, i, isInternalPoint);
        pointField& newPoints = tNewPoints.ref();

        // Orthogonal point blending
        if (orthogonalBlendingFraction > SMALL)
        {
            blendWithOrthogonalPoints(mesh, newPoints, uniValenceBoundaryMap, hasPointNormals, pointNormals, orthogonalBlendingFraction);
        }

        // Constrain absolute length of jump to new coordinates, to stabilize smoothing
        constrainMaxStepLength(mesh, newPoints, isInternalPoint, maxStepLength);

        // Avoid shortening of short edge length
        restrictEdgeShortening(mesh, newPoints, isInternalPoint, minEdgeLength, totalMinFreeze, isFrozenPoint);

        if (faceAngleConstraint)
        {
            // Restrict deterioration of face-face angles (WIP)
            restrictFaceAngleDeterioration(mesh, newPoints, isInternalPoint, minAngle, maxAngle, isFrozenPoint);
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

        // Push coordinates to mesh
        mesh.movePoints(tNewPoints);

        Info << endl;
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
