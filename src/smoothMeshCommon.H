/*---------------------------------------------------------------------------*\
Library
    Common definitions

Description
    Common definitions for source files

\*---------------------------------------------------------------------------*/

#ifndef smoothMeshCommon_H
#define smoothMeshCommon_H

// Macros for value definitions
#define UNDEF_LABEL -1
#define UNDEF_VECTOR vector(GREAT, GREAT, GREAT)
#define ZERO_VECTOR vector(0, 0, 0)
#define ZERO 0.0

// General safety factors or tolerances
#define REL_TOL 1e-4
#define ABS_TOL 1e-6

// Option to write extra debugging files
#define WRITE_DEBUG_FILES true

// Boolean for developer mode. Set to false to use bleeding edge
// work in progress features.
// #define USE_STABLE_FEATURES_ONLY true

// #include <typeinfo>
// Typeinfo is needed only for getting types while debugging, for example:
// Info << "Type is " << typeid(x).name() << endl;
// Use terminal command like this to demangle the mangled name:
// c++filt -t N4Foam4faceE

// #include "smoothMeshCommon.C"

using namespace Foam;

// Help function to project point pt to a line defined by startPoint
// and endPoint

vector projectPointToLine
(
    const point pt,
    const point startPoint,
    const point endPoint
)
{
    const double edgeLength = mag(endPoint - startPoint);
    if (edgeLength < VSMALL)
    {
        FatalError << "startPoint " << startPoint << " and endPoint " << endPoint << " overlap" << endl << abort(FatalError);
    }

    const point c2pt = pt - startPoint;
    const point edgeVec = endPoint - startPoint;
    const double normalizedDotProd = (c2pt & edgeVec) / sqr(edgeLength);
    const point projPoint = startPoint + normalizedDotProd * edgeVec;

    return projPoint;
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


// Help function to get boundary face normal vector. This is needed
// because boundary face info must be accessed through patches.

vector getBoundaryNf
(
    const fvMesh& mesh,
    const label faceI
)
{
    forAll(mesh.boundaryMesh(), patchI)
    {
        const fvPatch& patch = mesh.boundary()[patchI];
        const label startI = patch.start();
        const label endI = startI + patch.size() - 1;

        if ((faceI >= startI) and (faceI <= endI))
        {
            const label patchFaceI = faceI - startI;
            const vector Sf = patch.Sf()[patchFaceI];
            const vector nf = -Sf / patch.magSf()[patchFaceI];
            return nf;
        }
    }

    FatalError << "Did not find faceI " << faceI << " in patches" << endl << abort(FatalError);

    return vector(Zero);
}


// Help function to get boundary face center coordinates. This is
// needed because boundary face info must be accessed through patches.

vector getBoundaryCf
(
    const fvMesh& mesh,
    const label faceI
)

{
    forAll(mesh.boundaryMesh(), patchI)
    {
        const fvPatch& patch = mesh.boundary()[patchI];
        const label startI = patch.start();
        const label endI = startI + patch.size() - 1;

        if ((faceI >= startI) and (faceI <= endI))
        {
            const label patchFaceI = faceI - startI;
            return patch.Cf()[patchFaceI];
        }
    }

    FatalError << "getBoundaryCf did not find faceI " << faceI << " in patches" << endl << abort(FatalError);

    return UNDEF_VECTOR;
}


// Help function to export an edge, defined as the two argument
// points, as a string defining a STL facet. Since STL only handles
// triangles the edge is represented as a sliver triangle.

string edgeToStl
(
    const point p1,
    const point p2
)
{
    std::stringstream facet;
    facet
        << "facet normal 0 0 0\n"
        << " outer loop\n"
        << "  vertex "
        << p1[0] << " "
        << p1[1] << " "
        << p1[2] << "\n"
        << "  vertex "
        << p2[0] << " "
        << p2[1] << " "
        << p2[2] << "\n"
        << "  vertex "
        << p1[0] * (1.0 + ABS_TOL) + ABS_TOL << " "
        << p1[1] * (1.0 + ABS_TOL) + ABS_TOL << " "
        << p1[2] * (1.0 + ABS_TOL) + ABS_TOL << "\n"
        << " endloop\n"
        << "endfacet\n";

    return facet.str();
}

// Help function to export edges using given vector values (point
// locations) and edge end point labels to an STL file

void writeDebugStlEdges
(
    const string& filename,
    const vectorList& vectorValues,
    const labelList& pointLabels1,
    const labelList& pointLabels2,
    const labelList& pointLabels3
)
{
    if (! WRITE_DEBUG_FILES)
        return;

    std::ofstream myfile;
    myfile.open(filename + "_proc" + std::to_string(Pstream::myProcNo()) + ".stl");
    myfile << "solid " << filename << "\n";

    forAll(vectorValues, pointI)
    {
        const point p = vectorValues[pointI];

        if (pointLabels1[pointI] != UNDEF_LABEL)
        {
            myfile << edgeToStl(p, vectorValues[pointLabels1[pointI]]);
        }

        if (pointLabels2[pointI] != UNDEF_LABEL)
        {
            myfile << edgeToStl(p, vectorValues[pointLabels2[pointI]]);
        }

        if (pointLabels3[pointI] != UNDEF_LABEL)
        {
            myfile << edgeToStl(p, vectorValues[pointLabels3[pointI]]);
        }
    }

    myfile << "endsolid\n";
    myfile.close();
}


// Help function to count boolean true values for point lists

scalar countBoolListValues
(
    const boolList bList,
    const labelList nProcessorsOnPoint
)
{
    scalar nVals = 0.0;

    forAll (bList, i)
    {
        if (bList[i])
        {
            nVals += 1.0 / nProcessorsOnPoint[i];
        }
    }

    return returnReduce(nVals, sumOp<scalar>());
}

#endif
