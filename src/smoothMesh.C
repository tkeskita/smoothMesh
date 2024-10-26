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

// Help function to return distance to mesh point from argument
// coordinates

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
// locations. Edge points are fully frozen below minEdgeLength, and
// fully free to move above maxEdgeLength, with linear interpolation
// in between. This feature can be used to limit the squishing or
// compression of cells near concave features.

int restrictEdgeShortening
(
    const fvMesh& mesh,
    pointField& origPoints,
    const bitSet isMovingPoint,
    const double minEdgeLength,
    const double maxEdgeLength
)
{
    // Copy original points for temporary working point field
    tmp<pointField> tNewPoints(new pointField(mesh.nPoints(), Zero));
    // tmp<pointField> tNewPoints(new pointField(origPoints)); // TODO: test, not good?
    pointField& newPoints = tNewPoints.ref();

    forAll(newPoints, pointI)
        newPoints[pointI] = origPoints[pointI];

    forAll(origPoints, pointI)
    {
        if (! isMovingPoint.test(pointI))
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

        // Blending fraction of current (0.0) and new (1.0) point coordinates
        double frac = 1.0;

        // Consider freezing only if edge length decreases
        if (shortestNewEdgeLength < shortestCurrentEdgeLength)
        {
            // Full freeze below minEdgeLength
            if (shortestNewEdgeLength < minEdgeLength)
            {
                frac = 0.0;
            }

            // Partial freeze below maxEdgeLength
            else if (shortestNewEdgeLength < maxEdgeLength)
            {
                frac = (shortestNewEdgeLength - minEdgeLength) / (maxEdgeLength - minEdgeLength);
            }

            // No freeze above maxEdgeLength
            else
            {
                frac = 1.0;
            }
        }

        // Blend current and new coordinates
        const vector newCoords = ((1.0 - frac) * cCoords) + (frac * nCoords);

        // Save the constrained point
        newPoints[pointI] = newCoords;
    }

    // Save new point coordinates to original point field
    forAll(origPoints, pointI)
    {
        origPoints[pointI] = newPoints[pointI];
    }

    return 0;
}

// Constrain the length of a step jump to new coordinates by fraction
// of minimum edge length. This increases the stability of the
// relaxation process, in case target coordinates are far off.

int constrainLocalStepLength
(
    const fvMesh& mesh,
    pointField& origPoints,
    const bitSet isMovingPoint,
    const double relaxationFactor
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

        // Calculate shortest edge length
        double shortestEdgeLength = DBL_MAX;
        forAll(mesh.pointPoints(pointI), pointPpI)
        {
            const label neighI = mesh.pointPoints(pointI)[pointPpI];
            const double testLength = getPointDistance(mesh, neighI, cCoords);
            if (testLength < shortestEdgeLength)
            {
                shortestEdgeLength = testLength;
            }
        }

        // Scale down the length of the jump from current coordinates
        // towards new coordinates if jump would be otherwise too long
        const vector stepDir = nCoords - cCoords;
        const double stepLength = mag(stepDir);
        const double maxLength = relaxationFactor * shortestEdgeLength;
        if (stepLength > maxLength)
        {
            const double scale = maxLength / stepLength;
            nCoords = cCoords + scale * stepDir;
            // Info << "pointI " << pointI << " maxLength " << maxLength << " stepLength "
            //      << stepLength << " scale " << scale << endl;
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

// Constrain the length of a step jump to new coordinates by an
// absolute length value. This increases the stability of the
// relaxation process, in case target coordinates are far off.

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
            // Info << "pointI " << pointI << " maxLength " << maxLength << " stepLength "
            //      << stepLength << " scale " << scale << endl;
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

// Calculate and return the edge-edge angle (in radians) of two edges
// which share a common point cCoords. The two edge end point
// coordinates are p1Coords and p2Coords. fCoords is the face center
// coordinates, which is needed for deducing if the edge-edge angle is
// convex (<180 deg) or concave (>180 deg).
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

    // Info << "   Edge angle " << angle << " coords: " << cCoords << " - "
    //     << p1Coords << " - " << p2Coords << " - " << fCoords << endl;

    return angle;
}

// Find the neighbour point indices neighPI1 and neighPI2 for the
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
    // Info << "facePoints are " << facePoints << endl;
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

// Compute and return variance for argument vals (array of doubles).
// Array length is n.
double calcVariance
(
    double* vals,
    label n
)
{
    if (n == 0)
        FatalError << "n is zero" << abort(FatalError);

    double mean_val = 0.0;
    for (label i = 0; i < n; i++)
    {
        mean_val += vals[i];
    }
    mean_val /= static_cast<double>(n);

    double variance = 0.0;
    for (label i = 0; i < n; i++)
    {
        const double d = vals[i] - mean_val;
        variance += d * d;
    }
    variance /= static_cast<double>(n);

    return variance;
}


// Calculate the edge angles of point with index pointI using the
// current mesh point locations
bool calc_edge_angles0
(
    const fvMesh& mesh,
    const label pointI,
    bitSet &isFrozenPoint,
    double *angles0
)
{
    bool isConcavePoint = false; // marker for concave points
    const vector cCoords = mesh.points()[pointI];

    forAll(mesh.pointFaces()[pointI], pointFI)
    {
        // For each edge connected at this point, find the other point
        label neighPI1 = 0;
        label neighPI2 = 0;
        const label faceI = mesh.pointFaces()[pointI][pointFI];
        getNeighbourPoints(mesh, pointI, faceI, &neighPI1, &neighPI2);

        // Use current mesh point locations to calculate current angle variance
        const double edgeAngle = edgeEdgeAngle(cCoords, mesh.points()[neighPI1], mesh.points()[neighPI2]);

        angles0[pointFI] = edgeAngle;

        // If any edge-edge angle > pi radians, then the edges
        // form a concave angle, and the point can't be
        // frozen, otherwise concavity will just increase
        // during smoothing.
        if (edgeAngle > 3.14)
            isConcavePoint = true;

        // Info << "Indices " << pointI << ":" << neighPI1 << "-"
        //    << neighPI2 << " - " << edgeAngle << endl;
    }

    return isConcavePoint;
}

// Calculate the edge angles of point with index pointI using the
// current mesh point locations
int calc_edge_angles1
(
    const fvMesh& mesh,
    const pointField& newPoints,
    const label pointI,
    bitSet &isFrozenPoint,
    double *angles1
)
{
    forAll(mesh.pointFaces()[pointI], pointFI)
    {
        // For both edges connected at this point, find the other point index
        label neighPI1 = 0;
        label neighPI2 = 0;
        const label faceI = mesh.pointFaces()[pointI][pointFI];
        getNeighbourPoints(mesh, pointI, faceI, &neighPI1, &neighPI2);

        // Use neighbour point's new point coordinates, unless
        // the neighbour point is frozen, in which case use it's current
        // mesh point coordinates.
        vector point1 = newPoints[neighPI1];
        if (isFrozenPoint.test(neighPI1))
            point1 = mesh.points()[neighPI1];

        vector point2 = newPoints[neighPI2];
        if (isFrozenPoint.test(neighPI2))
            point2 = mesh.points()[neighPI2];

        const vector nCoords = newPoints[pointI];
        const double edgeAngle = edgeEdgeAngle(nCoords, point1, point2);

        angles1[pointFI] = edgeAngle;
    }

    return 0;
}


// Restrict the increase of variance of edge-edge angles of internal
// faces. This is a mesh quality constraint aimed at avoiding
// decrease of point quality during the smoothing process.
// Returns normalized sum of final angle variances as a return value.
double restrictAngleVarianceIncrease
(
    const fvMesh& mesh,
    pointField& origPoints,
    const bitSet isMovingPoint
)
{
    double totVariance = 0.0;
    int nMovingPoints = 0;
    int nFrozenPoints = 0;
    int nFrozenPointsOld = -1;
    int nFrozenPointIters = 0;

    // Storage of frozen points
    bitSet isFrozenPoint(mesh.nPoints(), false);

    // Copy original points for temporary working point field
    tmp<pointField> tNewPoints(new pointField(mesh.nPoints(), Zero));
    pointField& newPoints = tNewPoints.ref();
    forAll(newPoints, pointI)
        newPoints[pointI] = origPoints[pointI];

    // Loop the quality checking and point freezing procedure until no
    // new points are frozen
    while (nFrozenPoints > nFrozenPointsOld)
    {
        nMovingPoints = 0;
        nFrozenPointIters++;
        nFrozenPointsOld = nFrozenPoints;

        forAll(origPoints, pointI)
        {
            if (! isMovingPoint.test(pointI))
                continue;
            if (isFrozenPoint.test(pointI))
                continue;

            nMovingPoints++;
            const vector cCoords = mesh.points()[pointI];
            vector nCoords = newPoints[pointI];

            const label nPointFaces = mesh.pointFaces()[pointI].size();
            bool isConcavePoint = false; // marker for concave points

            // Calculate edge-edge angles and their variance for the
            // current point coordinates
            double angles0[nPointFaces];
            isConcavePoint = calc_edge_angles0(mesh, pointI, isFrozenPoint, angles0);
            const double variance0 = calcVariance(angles0, nPointFaces);

            // Calculate edge-edge cosine angles and variance for the
            // proposed new point coordinates
            double angles1[nPointFaces];
            calc_edge_angles1(mesh, newPoints, pointI, isFrozenPoint, angles1);
            const double variance1 = calcVariance(angles1, nPointFaces);

            // Concave points are never frozen. For other points:
            // If variance increases when moving point to new coordinates,
            // then freeze the point to current coordinates, otherwise leave
            // the new point coordinates unchanged.
            if ((! isConcavePoint) and (variance1 > 0.1) and (variance1 > variance0))
            {
                nFrozenPoints++;
                isFrozenPoint.set(pointI);
                nCoords = cCoords;
                totVariance += variance0;
            }
            else
            {
                totVariance += variance1;
            }

            // Save the constrained point
            newPoints[pointI] = nCoords;
        }
    }

    Info << "Froze " << nFrozenPoints << "/" << mesh.nPoints() << " points in "
         << nFrozenPointIters << " internal loops" << endl;

    // Save new point coordinates to original point field
    forAll(origPoints, pointI)
    {
        if (origPoints[pointI] != newPoints[pointI])
            Info << "origPointCoords!=newPointCoords at index " << pointI << endl;
        origPoints[pointI] = newPoints[pointI];
    }

    if (nMovingPoints > 0)
        return (totVariance / static_cast<double>(nMovingPoints));
    else
        return 0.0;
}

// Calculate and store the minimum edge angles of point with index
// pointI for current mesh (minCAngleStorage) and a minimum for new
// mesh points (minNAngleStorage)
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


// Restrict decrease of smallest face angles when angle is below
// minAngle (in degrees). This is meant to avoid creation of
// self-intersections for concave features.
int restrictMinAngleDecrease
(
    const fvMesh& mesh,
    pointField& origPoints,
    const bitSet isMovingPoint,
    const double minAngleInDegrees
)
{
    // Copy original points to temporary point field
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

        // Calculate current minimum angle and new minimum angle
        double minCAngle;
        double minNAngle;
        calc_min_edge_angles(mesh, newPoints, pointI, &minCAngle, &minNAngle);

        // If minimum angle is below threshold and would decrease,
        // then don't move the point. Note: This allows angle to
        // increase, so points are not frozen permanently.
        const double smallAngle = M_PI * minAngleInDegrees / 180.0;

        if ((minNAngle < smallAngle) and (minNAngle < minCAngle))
        {
            nCoords = cCoords;
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
        "qualityControl",
        "bool",
        "Enable or disable all quality control features (default: enabled)"
    );

    argList::addOption
    (
        "minEdgeLength",
        "double",
        "A quality control feature: Edge length below which edge vertices are fully frozen (default 0.002)"
    );

    argList::addOption
    (
        "maxEdgeLength",
        "double",
        "A quality control feature: Edge length above which edge vertices are fully free to move (default: 1.001 * minEdegeLength)"
    );

    argList::addOption
    (
        "minAngle",
        "double",
        "A quality control feature: Edge-edge angle for internal faces below which vertices are fully frozen (in degrees, default: 45)"
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

    double minEdgeLength(0.002);
    args.readIfPresent("minEdgeLength", minEdgeLength);

    double maxEdgeLength(1.001 * minEdgeLength);
    args.readIfPresent("maxEdgeLength", maxEdgeLength);

    double minAngle(45);
    args.readIfPresent("minAngle", minAngle);

    bool qualityControl(true);
    args.readIfPresent("qualityControl", qualityControl);

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

    label centroidalIters(20);
    args.readIfPresent("centroidalIters", centroidalIters);

    if (centroidalIters == 0)
    {
        Info << "Use -centroidalIters option to specify the number "
             << "of iteration rounds. Doing nothing."
             << nl << endl
             << "End" << nl << endl;
        return 0;
    }

    // Carry out centroidal smoothing iterations
    for (label i = 0; i < centroidalIters; ++i)
    {
        tmp<pointField> tNewPoints = centroidalSmoothing(mesh, i, isInternalPoint);
        pointField& newPoints = tNewPoints.ref();

        // Orthogonal point blending
        if (orthogonalBlendingFraction > SMALL)
        {
            blendWithOrthogonalPoints(mesh, newPoints, uniValenceBoundaryMap, hasPointNormals, pointNormals, orthogonalBlendingFraction);
        }

        // Constrain absolute length of jump to new coordinates, to stabilize smoothing
        constrainMaxStepLength(mesh, newPoints, isInternalPoint, maxStepLength);

        // Additional constraints aiming to avoid quality issues near concave features
        if (qualityControl)
        {
            // const double totVariance = restrictAngleVarianceIncrease(mesh, newPoints, isInternalPoint);
            // Info << "Total normalized angle variance = " << totVariance << endl;

            // Old step length control disabled: Constrain the local step
            // length by fraction of shortest edge length
            // constrainLocalStepLength(mesh, newPoints, isInternalPoint, 0.05);

            // Avoid shortening of short edge length
            restrictEdgeShortening(mesh, newPoints, isInternalPoint, minEdgeLength, maxEdgeLength);

            // Restrict decrease of smallest face angle
            restrictMinAngleDecrease(mesh, newPoints, isInternalPoint, minAngle);
        }


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
