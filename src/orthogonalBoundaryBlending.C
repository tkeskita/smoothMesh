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

// Help function to find an store the orthogonal prism slot number
// which points from the outer point towards the inner point

void setOuterPointOrthogonalSlot
(
    const fvMesh& mesh,
    const label innerPointI,
    const label outerPointI,
    labelList& orthogonalPrismSlots,
    label& n,
    const labelList& innerPrismPointLabels1,
    const labelList& innerPrismPointLabels2,
    const labelList& innerPrismPointLabels3,
    std::ofstream& myfile,
    const bool exportEdgesAsStl
)
{
    if ((outerPointI != UNDEF_LABEL) and
        (orthogonalPrismSlots[outerPointI] == UNDEF_LABEL))
    {
        bool writeDebugEdge = false;

        if (innerPrismPointLabels1[outerPointI] == innerPointI)
        {
            orthogonalPrismSlots[outerPointI] = 1;
            ++n;
            writeDebugEdge = true;
        }
        if (innerPrismPointLabels2[outerPointI] == innerPointI)
        {
            orthogonalPrismSlots[outerPointI] = 2;
            ++n;
            writeDebugEdge = true;
        }
        if (innerPrismPointLabels3[outerPointI] == innerPointI)
        {
            orthogonalPrismSlots[outerPointI] = 3;
            ++n;
            writeDebugEdge = true;
        }

        if ((writeDebugEdge) and (exportEdgesAsStl))
        {
            const label i1 = innerPointI;
            const label i2 = outerPointI;
            if (i1 == i2)
                FatalError << "inner point " << i1 << " is same as outer point" << endl << abort(FatalError);
            myfile << edgeToStl(mesh.points()[i1], mesh.points()[i2]);
        }
    }
}


// Help function to identify prisms slots for orthogonal treatment

label identifyOrthogonalPrismSlots
(
    const fvMesh& mesh,
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
        myfile.open("debugOrthogonalPrismEdgesAsStl_proc" + std::to_string(Pstream::myProcNo()) + ".stl");
        myfile << "solid orthoEdgesAsStl\n";
    }

    forAll (mesh.points(), pointI)
    {
        // Skip non-front points
        if ((innerPrismPointLabels1[pointI] != UNDEF_LABEL) or
            (innerPrismPointLabels2[pointI] != UNDEF_LABEL) or
            (innerPrismPointLabels3[pointI] != UNDEF_LABEL))
        {
            continue;
        }

        // Go through outer points. If outer point is pointing to this
        // point, add the correct slot number for orthogonal
        // treatment.
        const label outerI1 = outerPrismPointLabels1[pointI];
        const label outerI2 = outerPrismPointLabels2[pointI];
        const label outerI3 = outerPrismPointLabels3[pointI];

        setOuterPointOrthogonalSlot(mesh, pointI, outerI1, orthogonalPrismSlots, n, innerPrismPointLabels1, innerPrismPointLabels2, innerPrismPointLabels3, myfile, exportEdgesAsStl);
        setOuterPointOrthogonalSlot(mesh, pointI, outerI2, orthogonalPrismSlots, n, innerPrismPointLabels1, innerPrismPointLabels2, innerPrismPointLabels3, myfile, exportEdgesAsStl);
        setOuterPointOrthogonalSlot(mesh, pointI, outerI3, orthogonalPrismSlots, n, innerPrismPointLabels1, innerPrismPointLabels2, innerPrismPointLabels3, myfile, exportEdgesAsStl);
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
