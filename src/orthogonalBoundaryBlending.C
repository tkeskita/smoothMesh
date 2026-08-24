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
#include "smoothMeshCommon.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Help function which returns index of the slot for the only
// labelList with a defined value in the given point index

label getUniqueSlotIndex
(
    const label pointI,
    const labelList& pointLabels1,
    const labelList& pointLabels2,
    const labelList& pointLabels3
)
{
    label n = 0;
    label slotI = 0;

    if (pointLabels1[pointI] != UNDEF_LABEL)
    {
        ++n;
        slotI = 1;
    }
    if (pointLabels2[pointI] != UNDEF_LABEL)
    {
        ++n;
        slotI = 2;
    }
    if (pointLabels3[pointI] != UNDEF_LABEL)
    {
        ++n;
        slotI = 3;
    }

    if (n == 1)
    {
        return slotI;
    }

    return UNDEF_LABEL;
}


// Help function to identify prisms slots for orthogonal treatment.
// Only unique prismatic edges are considered, so there is maximum one
// orthogonal point label per point.

label identifyOrthogonalPrismSlots
(
    const fvMesh& mesh,
    const label i,
    labelList& orthogonalPrismSlots,
    const labelList& innerPrismPointLabels1,
    const labelList& innerPrismPointLabels2,
    const labelList& innerPrismPointLabels3,
    const labelList& outerPrismPointLabels1,
    const labelList& outerPrismPointLabels2,
    const labelList& outerPrismPointLabels3
)
{
    label n = 0;

    // Debug option to set true for printing edges as a STL file
    // Best visualized as wireframe in Paraview
    const bool exportEdgesAsStl = true;
    std::ofstream myfile;
    if (exportEdgesAsStl)
    {
        myfile.open("debugOrthogonalPrismEdgesAsStl_proc" + std::to_string(Pstream::myProcNo()) + "_layer_" + std::to_string(i) + ".stl");
        myfile << "solid orthoEdgesAsStl\n";
    }

    forAll (mesh.points(), pointI)
    {
        // Skip point if it is processed already, or if it does not
        // have exactly one internal point
        const label slotI = getUniqueSlotIndex(pointI, innerPrismPointLabels1, innerPrismPointLabels2, innerPrismPointLabels3);

        if ((orthogonalPrismSlots[pointI] != UNDEF_LABEL) or
            (slotI == UNDEF_LABEL))
        {
            continue;
        }

        orthogonalPrismSlots[pointI] = slotI;
        ++n;

        if (exportEdgesAsStl)
        {
            const label i1 = pointI;
            label i2;
            if (slotI == 1)
                i2 = innerPrismPointLabels1[pointI];
            else if (slotI == 2)
                i2 = innerPrismPointLabels2[pointI];
            else if (slotI == 3)
                i2 = innerPrismPointLabels3[pointI];

            if (i1 == i2)
                FatalError << "inner point " << i1 << " is same as outer point" << endl << abort(FatalError);
            myfile << edgeToStl(mesh.points()[i1], mesh.points()[i2]);
        }
    }

    if (exportEdgesAsStl)
    {
        myfile << "endsolid\n";
        myfile.close();
    }

    return n;
}


// Help function to project point orthogonally to a plane

point projectPointToPlane
(
    const point p,
    const point planePoint,
    const vector normalVec
)
{
    if (mag(p - planePoint) < ABS_TOL)
    {
        FatalError << "Point " << p << " and plane point " << planePoint << " are too close to each other" << endl << abort(FatalError);
    }

    const double c1 = (p - planePoint) & normalVec;
    if (fabs(c1) < ABS_TOL)
    {
        FatalError << "Line vector for point " << p << " and plane point " << planePoint << ", and normal vector " << normalVec << " are almost the same" << endl << abort(FatalError);
    }

    const vector target = p - c1 * normalVec;
    return target;
}


// Help function for testing does a point lie within a 3D triangle of
// points. Reference: W. Heidrich, Journal of Graphics, GPU, and Game
// Tools, Volume 10, Issue 3, 2005.

bool isPointWithinTriangle
(
    const point& p,
    const point& p1,
    const point& p2,
    const point& p3
)
{
    const vector u = p2 - p1;
    const vector v = p3 - p1;
    const vector w = p - p1;

    // In OpenFOAM, cross product is ^ and inner product (dot product) is &

    const vector n = u ^ v;
    if (magSqr(n) < VSMALL)
    {
        FatalError << "magSqr(n) < VSMALL" << endl << abort(FatalError);
    }
    const vector npn = n / magSqr(n);

    // Barycentric coordinates
    const double gamma = (u ^ w) & npn;
    const double beta = (w ^ v) & npn;
    const double alpha = 1.0 - beta - gamma;

    // The test
    return
    (
        (alpha >= 0) and (alpha <= 1.0) and
        (beta >= 0) and (beta <= 1.0) and
        (gamma >= 0) and (gamma <= 1.0)
    );
}

// Help function to calculate normalized normal vector from triangle
// points

vector calcTriangleNormal
(
    const point& p1,
    const point& p2,
    const point& p3
)
{
    const vector i = p2 - p1;
    const vector j = p3 - p1;
    const vector n = i ^ j;
    if (mag(n) < VSMALL)
    {
        FatalError << "mag(n) < VSMALL" << endl << abort(FatalError);
    }

    return n / mag(n);
}

// Help function to project a point to a triangle of points

point projectPointToTriangle
(
    const point& p,
    const point& p1,
    const point& p2,
    const point& p3
)
{
    // Check for point equality
    if (p == p1) return p;
    if (p == p2) return p;
    if (p == p3) return p;

    // Project point to triangle plane
    const vector n = calcTriangleNormal(p1, p2, p3);
    const point pproj = projectPointToPlane(p, p1, n);

    // If point is on triangle, return it
    if (isPointWithinTriangle(pproj, p1, p2, p3))
    {
        return pproj;
    }

    // Otherwise project to triangle edges and choose closest
    const point pproj1 = projectPointToLine(p, p1, p2);
    const point pproj2 = projectPointToLine(p, p2, p3);
    const point pproj3 = projectPointToLine(p, p3, p1);

    const double len1 = mag(pproj1 - p);
    const double len2 = mag(pproj2 - p);
    const double len3 = mag(pproj3 - p);

    if ((len1 <= len2) and (len1 <= len3)) return pproj1;
    if ((len2 <= len1) and (len2 <= len3)) return pproj2;
    if ((len3 <= len2) and (len3 <= len1)) return pproj3;

    FatalError << "Sanity failed in projectPointToTriangle" << endl << abort(FatalError);

    return p;
}

// Function to improve orthogonality of prismatic boundary edges

int blendWithOrthoPoints
(
    const fvMesh& mesh,
    pointField& newPoints,
    const labelList& orthogonalPrismSlots,
    const vectorList& innerPrismPoints1,
    const vectorList& innerPrismPoints2,
    const vectorList& innerPrismPoints3,
    const vectorList& pointNormals1,
    const vectorList& pointNormals2,
    const vectorList& pointNormals3,
    const double orthoBlendingFraction,
    const double orthoMinLayers,
    const double orthoMaxLayers
)
{
    forAll(mesh.points(), pointI)
    {
        if (orthogonalPrismSlots[pointI] == UNDEF_LABEL)
            continue;

        vector orthoPoint(UNDEF_VECTOR);

        // Projection of inner point is done only for the unique
        // non-overlapping prismatic edges

        if (orthogonalPrismSlots[pointI] == 1)
        {
            orthoPoint = projectPointToPlane(innerPrismPoints1[pointI], newPoints[pointI], pointNormals1[pointI]);
        }

        else if (orthogonalPrismSlots[pointI] == 2)
        {
            orthoPoint = projectPointToPlane(innerPrismPoints2[pointI], newPoints[pointI], pointNormals2[pointI]);
        }

        else if (orthogonalPrismSlots[pointI] == 3)
        {
            orthoPoint = projectPointToPlane(innerPrismPoints3[pointI], newPoints[pointI], pointNormals3[pointI]);
        }

        // if (pointI == 109)
        //     Info<< "pointI109"
        //         << " orthogonalPrismSlots " << orthogonalPrismSlots[pointI]
        //         << " innerPrismPoints1 " << innerPrismPoints1[pointI]
        //         << " innerPrismPoints2 " << innerPrismPoints2[pointI]
        //         << " innerPrismPoints3 " << innerPrismPoints3[pointI]
        //         << " pointNormals1 " << pointNormals1[pointI]
        //         << " pointNormals2 " << pointNormals2[pointI]
        //         << " pointNormals3 " << pointNormals3[pointI]
        //         << " orthoPoint " << orthoPoint
        //         << endl;
                
        // Update point coordinates
        if (orthoPoint != UNDEF_VECTOR)
        {
            const point newPoint = newPoints[pointI];
            const vector blendedPoint = orthoBlendingFraction * orthoPoint +
                (1.0 - orthoBlendingFraction) * newPoint;
            newPoints[pointI] = blendedPoint;
        }
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
    vectorList& pointNormals,
    boolList& isSharpEdgePoint
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

        const fvPatch& patch = mesh.boundary()[patchI];
        const label startI = patch.start();
        const label endI = patch.size();

        // Add inversed unit normal vectors of patch faces to all
        // pointNormals of the face points
        for (label faceI = 0; faceI < endI; faceI++)
        {
            // Sf is unit normal vector multiplied by surface area, so
            // need to normalise it before use
            const vector cSf = patch.Sf()[faceI];
            vector Sf = cSf / patch.magSf()[faceI];

            const face& f = mesh.faces()[startI + faceI];
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


    // Classify and mark very sharp edge points
    forAll(pointNormals, pointI)
    {
        if (nFaces[pointI] < 1)
            continue;

        const double magNorm = mag(pointNormals[pointI]);

        // Zero the point normal for very sharp edge points
        if (magNorm < 0.1)
        {
            pointNormals[pointI] = ZERO_VECTOR;
            isSharpEdgePoint[pointI] = true;
        }
        else
        {
            isSharpEdgePoint[pointI] = false;
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

// ************************************************************************* //
