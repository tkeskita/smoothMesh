/*---------------------------------------------------------------------------*\
Library
    Orthogonal Boundary Blending

Description
    Special treatment of boundary cell layer, with aim to increase
    orthogonality of internal side faces on the cell layer.

    Idea is to find internal mesh points which have only one edge
    connection to the boundary points. The ideal direction for those
    edges is opposite of the boundary point vertex normal direction.
    The edge length is calculated by trigonometric projection from
    the current point location to the (inverted) vertex normal vector.
    The resulting point is called Orthogonal point coordinate.

    To avoid intersections caused by non-orthogonal meshes, the new
    point location is a weighted average of the point coordinate from
    main smoothing method and the Orthogonal point coordinate using
    a blending weight factor (orthogonalBlendingFraction).

IDEA: Extend orthogonality to also internal points with more than one
boundary edges. Orthogonality should be encouraged only to the
shortest edge. Possibly only if the edge length is 30%(?) or less of
the length of the second shortest edge.

IDEA: For multivalent edge points: Don't allow distance to boundary to
decrease. Distance calculation by assuming an edge for 2-valent and a
plane for 3-valent (or more) points.

\*---------------------------------------------------------------------------*/

#include "fvMesh.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Calculate point normals of boundary points starting from
// polyMesh. Stores point normals to pointNormals field, and marks the
// availability of point normal vector in hasPointNormal bit set.
//
// Uses primitivePatch.pointNormals() for calculating point normals,
// so points on the patch boundaries are excluded, as the normal
// direction can't be guaranteed to be correct when several patches
// share the point (e.g. corners can be non-trivial).

int calculatePointNormals
(
    const fvMesh& mesh,
    pointField& pointNormals,
    bitSet& hasPointNormal
)
{
    // Storage of already visited points
    bitSet visitedPoint(mesh.nPoints(), false);

    forAll(mesh.boundaryMesh(), patchI)
    {
        const polyPatch& pp = mesh.boundaryMesh()[patchI];
        forAll(pp.meshPoints(), pointPpI)
        {
            const label pointI = pp.meshPoints()[pointPpI];
            
            // If this point is shared by patches, remove it from
            // the set of points that have point normals.
            if (visitedPoint.test(pointI))
            {
                hasPointNormal.unset(pointI);
                pointNormals[pointI] = Zero;
            }
            else
            {
                visitedPoint.set(pointI);
                const vector pointNormal = -pp.pointNormals()[pointPpI];
                pointNormals[pointI] = pointNormal;
                hasPointNormal.set(pointI);
            }        
        }
    }

    // Make sure first point has no point normal, as it is used as
    // a no data marker.
    hasPointNormal.unset(0);
    pointNormals[0] = Zero;

    return 0;
}


// Calculate index label map from internal point to the boundary point
// for those internal mesh points which are directly connected to one
// and only one boundary point by an edge (points which have an
// univalent boundary edge).

int calculateUniValenceBoundaryMap
(
    const fvMesh& mesh,
    labelList& map,
    const bitSet& isInternalPoint
)
{
    forAll(map, pointI)
        // TODO: Can a label be -1? To distinguish from first point.
        map[pointI] = 0;

    forAll(mesh.points(), pointI)
    {
        if (isInternalPoint.test(pointI))
        {
            label valencePoints = 0;
            label boundaryPointLabel = 0;
            forAll(mesh.pointPoints(pointI), pointPpI)
            {
                const label neighI = mesh.pointPoints(pointI)[pointPpI];
                if (! isInternalPoint.test(neighI))
                {
                    valencePoints++;
                    boundaryPointLabel = neighI;
                }
            }
            
            if (valencePoints == 1)
            {
                map[pointI] = boundaryPointLabel;
            }
        }
    }

    return 0;
}

// Calculate the point coordinates for the orthogonally optimal point
// location and blend with the given new point coordinates

int blendWithOrthogonalPoints
(
    const polyMesh& mesh,
    pointField& newPoints,
    const labelList& uniValenceBoundaryMap,
    const bitSet& hasPointNormals,
    const pointField& pointNormals,
    const double orthogonalBlendingFraction
)
{
    for (label pointI = 0; pointI < mesh.nPoints(); pointI++)
    {
        const label boundaryPointI = uniValenceBoundaryMap[pointI];
        if (boundaryPointI == 0)
            continue;

        if (! hasPointNormals.test(boundaryPointI))
            FatalError << "Sanity broken, point " << boundaryPointI
                       << " should have a point normal" << endl
                       << abort(FatalError);

        const vector cCoords = mesh.points()[pointI];
        const vector bCoords = mesh.points()[boundaryPointI];
        const vector cb = cCoords - bCoords;
        const vector bNormal = pointNormals[boundaryPointI];
        const double dotProd = cb & bNormal;
        const vector orthoPoint = bCoords + dotProd * bNormal;
        const vector newPoint = newPoints[pointI];

        const vector blendedPoint = orthogonalBlendingFraction * orthoPoint +
            (1.0 - orthogonalBlendingFraction) * newPoint;

        // Update point coordinates
        newPoints[pointI] = blendedPoint;
    }

    return 0;
}

// ************************************************************************* //
