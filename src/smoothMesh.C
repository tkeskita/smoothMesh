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

// Prohibits decrease of shortest edge length. This is necessary to
// limit squishing or compression of cells near concave features.

int constrainEdgeLength
(
    const fvMesh& mesh,
    pointField& origPoints,
    const bitSet isMovingPoint,
    const double minimumEdgeLength
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
        vector nCoords = origPoints[pointI];

        // Calculate shortest edge length
        double shortestEdgeLength = DBL_MAX;
        forAll(mesh.pointPoints(pointI), pointPpI)
        {
            const label neighI = mesh.pointPoints(pointI)[pointPpI];
            const double testLength = getPointDistance(mesh, neighI, cCoords);
            if (testLength < shortestEdgeLength)
                shortestEdgeLength = testLength;
        }

        // Project the new coordinate position away from neighbour
        // points until the length to all neighbour points is bigger
        // than current shortest edge length.
        label i = 0;
        while (true)
        {
            // Info << "Constraining point " << pointI << " iter " << i << endl;
            bool innerConvergence = true;
            i++;
            forAll(mesh.pointPoints(pointI), pointPpI)
            {
                const label neighI = mesh.pointPoints(pointI)[pointPpI];
                const double testLength = getPointDistance(mesh, neighI, nCoords);
                if ((testLength < shortestEdgeLength) and
                    (testLength < minimumEdgeLength))
                {
                    // Info << "Constraining point " << pointI << " length "
                    //      << testLength << "<" << shortestEdgeLength << endl;
                    innerConvergence = false;
                    const vector pushBackDir = nCoords - mesh.points()[neighI];
                    const double scale = 1.01 * shortestEdgeLength / testLength;
                    nCoords = mesh.points()[neighI] + scale * pushBackDir;
                }
            }
            // Stop after a few iterations or if converged
            if (i > 2)
                break;
            if (innerConvergence == true)
                break;
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
        "Number of centroidal smoothing iterations"
    );

    argList::addOption
    (
        "maxStepLength",
        "double",
        "Maximum absolute step length applied in smoothing (default 0.001)"
    );

    argList::addOption
    (
        "orthogonalBlendingFraction",
        "double",
        "Fraction to force orthogonal side edges on the boundary (default 0.3)"
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

    double orthogonalBlendingFraction(0.3);
    args.readIfPresent("orthogonalBlendingFraction", orthogonalBlendingFraction);

    double maxStepLength(0.001);
    args.readIfPresent("maxStepLength", maxStepLength);

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

    // return 0;

    // Calculate new point locations and apply to mesh
    label centroidalIters(0);
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

        // Constrain absolute length of jump to new coordinates
        constrainMaxStepLength(mesh, newPoints, isInternalPoint, maxStepLength);

        // WIP disabled: constrain local step length by fraction of shortest edge length
        // constrainLocalStepLength(mesh, newPoints, isInternalPoint, 0.05);

        // WIP disabled: Constrain smoothing by shortest edge length
        // constrainEdgeLength(mesh, newPoints, isInternalPoint, minimumEdgeLength);

        mesh.movePoints(tNewPoints);
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
