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

#endif
