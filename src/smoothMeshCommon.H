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

// Help function to return a list of vector elements from a vector list

vectorList getVecListElements
(
    const vectorList& vecs,
    const labelList& labelIs
)
{
    vectorList elems;
    forAll (labelIs, i)
    {
        elems.append(vecs[labelIs[i]]);
    }
    return elems;
}


// Help function to calculate average vector of a list of vectors

vector calcMeanVec
(
    const vectorList& vecs
)
{
    vector r(Zero);
    if (vecs.size() == 0)
    {
        return r;
    }

    forAll (vecs, i)
    {
        r += vecs[i];
    }
    r /= vecs.size();

    return r;
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



// Help function for syncVecList to find and synchronize next vector
// value not existing in global vector list for all mesh points

int updateNextSyncVec
(
    const fvMesh& mesh,
    const List<vectorList>& localVecs,
    vectorList& syncVecs,
    const List<vectorList>& globalVecs
)
{
    forAll (mesh.points(), pointI)
    {
        // Default initialization: nothing to synchronize
        syncVecs[pointI] = UNDEF_VECTOR;

        const vectorList& globals = globalVecs[pointI];
        const vectorList& locals = localVecs[pointI];

        forAll(locals, localI)
        {
            const vector vec = locals[localI];

            // Write vector to syncVecs if it's not already in globals
            if (findIndex(globals, vec) == -1)
            {
                syncVecs[pointI] = vec;
                break;
            }
        }
    }

    return 0;
}


// Help function for syncVecList to count number of vector entries not
// yet in globals list for all mesh points.

int countRemainingNs
(
    const fvMesh& mesh,
    const List<vectorList>& localVecs,
    const List<vectorList>& globalVecs,
    labelList& ns
)
{
    forAll (mesh.points(), pointI)
    {
        // Default initialization: nothing to synchronize
        ns[pointI] = 0;

        const vectorList& globals = globalVecs[pointI];
        const vectorList& locals = localVecs[pointI];

        forAll(locals, localI)
        {
            const vector vec = locals[localI];

            if (findIndex(globals, vec) == -1)
            {
                ns[pointI] = locals.size() - localI;
                break;
            }
        }
    }

    return 0;
}

// Help function for syncVecList for updating the global vector list

int updateGlobalVecs
(
    const fvMesh& mesh,
    const vectorList& priorSyncVecs,
    const vectorList& syncVecs,
    List<vectorList>& globalVecs
)
{
    forAll (mesh.points(), pointI)
    {
        const vector& priorSyncVec = priorSyncVecs[pointI];
        const vector& syncVec = syncVecs[pointI];

        if (isSmallerByVectorElements(priorSyncVec, syncVec))
        {
            globalVecs[pointI].append(priorSyncVec);
        }
        else if (isSmallerByVectorElements(syncVec, priorSyncVec))
        {
            globalVecs[pointI].append(syncVec);
        }
    }

    return 0;
}


// Function to synchronize a local list of vectors at mesh points
// among parallel processes, to produce a global lists of unique
// vectors at mesh points. List elements are assumed to be ordered in
// smallest-by-component order.

int syncVecList
(
    const fvMesh& mesh,
    const List<vectorList>& localVecs,
    List<vectorList>& globalVecs
)
{
    // Synchronization work vectors
    vectorList syncVecs(mesh.nPoints(), UNDEF_VECTOR);

    // Number of local elements remaining to be synchronized for points
    labelList ns(mesh.nPoints(), UNDEF_LABEL);

    // Maximum number of local elements remaining to be synchronized
    label maxNs = 1;

    while (maxNs > 0)
    {
        Pout << "  syncVecList maxNs " << maxNs << endl;
        // Find and synchronize next unsynced vector for mesh points
        updateNextSyncVec(mesh, localVecs, syncVecs, globalVecs);
        const vectorList priorSyncVecs(syncVecs);

        syncTools::syncPointList
        (
            mesh,
            syncVecs,
            eqOp<vector>(),
            UNDEF_VECTOR               // null value
        );

        updateGlobalVecs(mesh, priorSyncVecs, syncVecs, globalVecs);
        countRemainingNs(mesh, localVecs, globalVecs, ns);
        maxNs = returnReduce(max(ns), maxOp<label>());
    }

    return 0;
}

#endif
