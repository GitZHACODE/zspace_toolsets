# zTs3DP Slicer Pseudo-Code

This is the first porting target for replacing the current placeholder slice generation with the mesh-walking method from the blend sketch and the print-path diagrams.

## Data

```cpp
struct z3DPSlicingFeatures
{
    zIntArray cornerVertexIds;       // boundary vertices with valence 3
    zIntArray boundaryLoopA;         // first boundary edge loop
    zIntArray boundaryLoopB;         // second boundary edge loop
    zIntArray shortestEdgeLoop;      // le in the paper diagram
    zIntArray topStripFaceIds;       // light-blue top quad strip
    zIntArray bottomStripFaceIds;    // light-blue bottom quad strip
    zIntArray slicingLoopHalfedges;  // halfedge ids used to produce slice meshes
};
```

## Feature Detection

```cpp
z3DPSlicingFeatures identifySlicingFeatures(zObjectMesh& mesh)
{
    zFnMesh fn(mesh);
    z3DPSlicingFeatures out;

    // 1. Find the four corner vertices.
    for (zItMeshVertex v(mesh); !v.end(); v++)
    {
        if (!v.onBoundary()) continue;
        if (v.getValence() == 3) out.cornerVertexIds.push_back(v.getId());
    }

    // 2. Walk boundary loops between valence-3 corners.
    //    Each candidate loop starts at a boundary halfedge whose start/end
    //    touches a corner vertex.
    for (zItMeshHalfEdge he(mesh); !he.end(); he++)
    {
        if (!he.onBoundary()) continue;
        if (!touchesCorner(he, out.cornerVertexIds)) continue;

        zIntArray loop = walkBoundaryCornerToCorner(he, out.cornerVertexIds);
        storeUniqueBoundaryLoop(loop, out.boundaryLoopA, out.boundaryLoopB);
    }

    // 3. Pick slicing direction from the shorter boundary loop.
    //    This is le in the diagram, where N = ceil(le / printLayerHeight).
    float lengthA = polylineLength(out.boundaryLoopA, mesh);
    float lengthB = polylineLength(out.boundaryLoopB, mesh);
    out.shortestEdgeLoop = (lengthA <= lengthB) ? out.boundaryLoopA : out.boundaryLoopB;

    // 4. Identify top and bottom quad strips.
    //    For each boundary halfedge on the two long side loops, collect the
    //    adjacent face. Classify the two strips by average normal/height or
    //    by which side of the shortest loop they lie on.
    for (int heId : boundaryHalfedges(out.boundaryLoopA, mesh))
    {
        zItMeshHalfEdge he(mesh, heId);
        if (!he.getSym().onBoundary()) addFaceIfQuad(he.getSym().getFace(), out.topStripFaceIds);
    }

    for (int heId : boundaryHalfedges(out.boundaryLoopB, mesh))
    {
        zItMeshHalfEdge he(mesh, heId);
        if (!he.getSym().onBoundary()) addFaceIfQuad(he.getSym().getFace(), out.bottomStripFaceIds);
    }

    orientTopBottomByAverageZOrNormal(out.topStripFaceIds, out.bottomStripFaceIds, mesh);
    return out;
}
```

## Boundary Walk

```cpp
zIntArray walkBoundaryCornerToCorner(zItMeshHalfEdge start, const zIntArray& cornerIds)
{
    zIntArray vertexLoop;
    zItMeshHalfEdge he = start;
    int guard = 0;

    while (guard++ < start.size())
    {
        int from = he.getStartVertex().getId();
        int to = he.getVertex().getId();

        if (vertexLoop.empty()) vertexLoop.push_back(from);
        vertexLoop.push_back(to);

        if (vertexLoop.size() > 1 && contains(cornerIds, to)) break;

        // Boundary traversal follows the next boundary halfedge around the rim.
        // Depending on orientation this is usually:
        // he = he.getNext()
        // or:
        // he = he.getPrev().getSym()
        he = nextBoundaryHalfedge(he);

        if (he.getId() == start.getId()) break;
    }

    return vertexLoop;
}
```

## Slice Mesh Generation

```cpp
void computeSliceMeshes(zObjectMesh& mesh, const z3DPSlicingFeatures& features, float printLayerHeight)
{
    float le = polylineLength(features.shortestEdgeLoop, mesh);
    int N = ceil(le / printLayerHeight);

    // 1. Split the shortest loop into N stations.
    zPointArray stations = dividePolylineByLength(features.shortestEdgeLoop, N, mesh);

    // 2. For every station, walk across the quad mesh in the transverse direction.
    //    This is the same halfedge family as zFnMesh::computeEdgeLoop:
    //    next = he.getNext().getSym().getNext()
    for (int i = 0; i <= N; i++)
    {
        zItMeshHalfEdge seed = halfedgeNearestStation(stations[i], features.shortestEdgeLoop, mesh);
        zItMeshHalfEdgeArray crossLoop;
        zFnMesh(mesh).computeEdgeLoop(seed, crossLoop);

        // 3. Intersect/interpolate that edge loop with the station parameter
        //    and collect a quad-dominant slice mesh.
        zObjectMesh sliceMesh;
        buildSliceMeshFromCrossLoop(crossLoop, stations[i], sliceMesh);
        sliceMeshes.push_back(sliceMesh);
    }
}
```

## Unroll Hook

```cpp
void computeUnrolledSlices()
{
    for each sliceMesh:
        triangulate if needed
        compute UV / planar parameterization
        store planar mesh in unrolledSliceMeshes
}
```

The next implementation step should expose `identifySlicingFeatures()` on `z3DPSlicer`, then draw the corner vertices plus top/bottom strips in the Alice sketch before generating real slice meshes.
