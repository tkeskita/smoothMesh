# smoothMesh

OpenFOAM mesh smoothing tool to improve mesh quality.

Warning: Work in progress.

## Current restrictions

- Smoothing of internal mesh points only
- Applies centroidal smoothing algorithm
- Developed on OpenFOAM.com v2312

## Compilation instructions

```
. /usr/lib/openfoam/openfoam2312/etc/bashrc
cd smoothMesh/src
wclean
wmake
```

## Usage examples

- Parallel run example: `mpirun -np 3 smoothMesh -centroidalIters 50 -parallel`
- Serial run example: smoothMesh -centroidalIters 50`
