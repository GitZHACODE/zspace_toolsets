#ifndef EIGEN_DONT_PARALLELIZE
#define EIGEN_DONT_PARALLELIZE
#endif

#include <zspace/zToolsets/fabrication/slicing/utils/Parallel.h>
#include <zspace/zToolsets/fabrication/slicing/zSlicer/HarmonicSlicer.h>

#ifndef ZSPACE_INLINE
#define ZSPACE_INLINE ZSPACE_TOOLSETS_INLINE
#endif

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <utility>

namespace zSpace {

using namespace zTsUtils;

struct HarmonicSlicer::ContourData {
  zPointArray Positions, UV;
  zFloatArray Heights;
  zIntArray Counts, Connects;
  std::vector<std::pair<int, int>> TetEdges;
  std::unordered_map<unsigned long long, int> EdgeVertices;
};

ZSPACE_INLINE bool HarmonicSlicer::computeVertexLayerHeights(int numLayers, zFloatArray &heights) const {
  heights.clear();
  if (numLayers < 2 || !_Tet.FieldValid ||
      _Tet.VolumeFieldGradients.size() != _Tet.Volumes.size()) return false;
  // Recover a continuous nodal height from piecewise-constant tetrahedral
  // gradients. This is local normal spacing, not a global Z increment or a
  // bound imposed by LayerSpacing. Both caps use the same full iso interval.
  std::vector<double> weighted(_Tet.Vertices.size(), 0.0), weights(weighted.size(), 0.0);
  const double delta = 1.0 / (numLayers - 1);
  for (size_t t = 0; t < _Tet.Volumes.size(); ++t) {
    const auto &volume = _Tet.Volumes[t];
    if (!volume.Valid) continue;
    const double gradient = _Tet.VolumeFieldGradients[t].norm();
    if (!std::isfinite(gradient) || gradient <= 0.0) return false;
    const double height = delta / gradient;
    if (!std::isfinite(height) || height <= 0.0 ||
        !std::isfinite(volume.Volume) || volume.Volume <= 0.0) return false;
    for (int id : volume.VertexIds) {
      weighted[id] += volume.Volume * height;
      weights[id] += volume.Volume;
    }
  }
  heights.resize(weighted.size());
  for (size_t i = 0; i < heights.size(); ++i) {
    if (!std::isfinite(weights[i]) || weights[i] <= 0.0) {
      heights.clear();
      return false;
    }
    heights[i] = static_cast<float>(weighted[i] / weights[i]);
    if (!std::isfinite(heights[i]) || heights[i] <= 0.0f) {
      heights.clear();
      return false;
    }
  }
  return !heights.empty();
}

ZSPACE_INLINE int HarmonicSlicer::crossingVertex(int a, int b, float iso,
    const zFloatArray &vertexHeights, ContourData &slice) {
  // Classify and interpolate with the same precision as the nodal field.
  // Rounding phi/t to float can collapse two distinct edge intersections onto
  // one endpoint, discard a triangle and leave an untagged hairline boundary.
  const double aPhi = _Tet.UVW[a].z;
  const double bPhi = _Tet.UVW[b].z;
  const double denominator = bPhi - aPhi;
  const double t = denominator == 0 ? 0.5 : (iso - aPhi) / denominator;
  const double tc = std::clamp(t, 0.0, 1.0);
  const int sourceA = tc == 1.0f ? b : a;
  const int sourceB = tc == 0.0f ? a : b;
  const unsigned long long key = makeEdgeKey(sourceA, sourceB);
  const auto found = slice.EdgeVertices.find(key);
  if (found != slice.EdgeVertices.end()) return found->second;
  const int id = static_cast<int>(slice.Positions.size());
  slice.Positions.push_back(_Tet.Vertices[a].Position + (_Tet.Vertices[b].Position - _Tet.Vertices[a].Position) * tc);
  slice.UV.push_back(_Tet.UVW[a] + (_Tet.UVW[b] - _Tet.UVW[a]) * tc);
  slice.Heights.push_back((1.0f - tc) * vertexHeights[a] + tc * vertexHeights[b]);
  // An endpoint intersection belongs to the vertex, rather than the whole edge.
  slice.TetEdges.emplace_back(sourceA, sourceB);
  slice.EdgeVertices.emplace(key, id);
  return id;
}

ZSPACE_INLINE bool HarmonicSlicer::appendTriangle(zIntArray &connects, int a, int b, int c,
    const zPointArray &positions, const zVector &lowToHigh,
    std::unordered_set<std::array<int, 3>, TriangleKeyHash> &uniqueTriangles) {
  if (a == b || b == c || c == a) return false;

  // Keep small triangles: an absolute area cutoff can open holes in the slice.
  // Compute differences and the normal in double precision before testing degeneracy.
  const Eigen::Vector3d pa = toEigenVec(positions[a]);
  const Eigen::Vector3d pb = toEigenVec(positions[b]);
  const Eigen::Vector3d pc = toEigenVec(positions[c]);
  if (!pa.allFinite() || !pb.allFinite() || !pc.allFinite()) return false;
  const Eigen::Vector3d normal = (pb - pa).cross(pc - pa);
  if (!normal.allFinite() || normal.cwiseAbs().maxCoeff() == 0.0) return false;
  const Eigen::Vector3d direction(lowToHigh.x, lowToHigh.y, lowToHigh.z);
  if (normal.dot(direction) < 0.0) std::swap(b, c);

  std::array<int, 3> triangleKey = {a, b, c};
  std::sort(triangleKey.begin(), triangleKey.end());
  if (!uniqueTriangles.insert(triangleKey).second) {
    return false;
  }

  connects.push_back(a);
  connects.push_back(b);
  connects.push_back(c);
  return true;
}

ZSPACE_INLINE void HarmonicSlicer::extractIsoMesh(float iso, zObjMesh &mesh, zPointArray &uv,
    const zFloatArray &vertexHeights, zFloatArray &heights,
    const WallWidthLookup &lookup, EdgeWidthIds &edgeIds, std::vector<std::string> &debug) {
  ContourData slice;
  std::unordered_set<std::array<int, 3>, TriangleKeyHash> uniqueTriangles;
  for (const auto &volume : _Tet.Volumes) {
    const auto &tet = volume.VertexIds;
    std::array<int, 4> inside{}, outside{};
    int insideCount = 0, outsideCount = 0;
    for (int localId = 0; localId < 4; ++localId) {
      if (_Tet.UVW[tet[localId]].z < iso) inside[insideCount++] = tet[localId];
      else outside[outsideCount++] = tet[localId];
    }
    if (insideCount == 0 || outsideCount == 0) continue;

    zVector lowToHigh;
    for (int i = 0; i < insideCount; ++i) {
      for (int j = 0; j < outsideCount; ++j) lowToHigh += _Tet.Vertices[outside[j]].Position - _Tet.Vertices[inside[i]].Position;
    }

    if (insideCount == 1 || outsideCount == 1) {
      const bool oneInside = insideCount == 1;
      const int lone = oneInside ? inside[0] : outside[0];
      const auto &other = oneInside ? outside : inside;
      const int a = crossingVertex(lone, other[0], iso, vertexHeights, slice);
      const int b = crossingVertex(lone, other[1], iso, vertexHeights, slice);
      const int c = crossingVertex(lone, other[2], iso, vertexHeights, slice);
      if (appendTriangle(slice.Connects, a, oneInside ? b : c, oneInside ? c : b,
          slice.Positions, lowToHigh, uniqueTriangles)) slice.Counts.push_back(3);
    } else {
      const int a = crossingVertex(inside[0], outside[0], iso, vertexHeights, slice);
      const int b = crossingVertex(inside[0], outside[1], iso, vertexHeights, slice);
      const int c = crossingVertex(inside[1], outside[1], iso, vertexHeights, slice);
      const int d = crossingVertex(inside[1], outside[0], iso, vertexHeights, slice);
      if (appendTriangle(slice.Connects, a, b, c, slice.Positions, lowToHigh, uniqueTriangles)) slice.Counts.push_back(3);
      if (appendTriangle(slice.Connects, a, c, d, slice.Positions, lowToHigh, uniqueTriangles)) slice.Counts.push_back(3);
    }
  }

  zFnMesh fn(mesh);
  fn.clear();
  uv.clear();
  heights.clear();
  edgeIds.clear();
  if (!slice.Counts.empty()) {
    // A crossing may belong only to a discarded degenerate triangle. Keep
    // vertices supported by surviving faces, with all provenance aligned.
    std::vector<int> remap(slice.Positions.size(), -1);
    for (int id : slice.Connects) remap[id] = 0;
    int used = 0;
    for (size_t id = 0; id < remap.size(); ++id) if (remap[id] >= 0) {
      remap[id] = used;
      slice.Positions[used] = slice.Positions[id];
      slice.UV[used] = slice.UV[id];
      slice.Heights[used] = slice.Heights[id];
      slice.TetEdges[used] = slice.TetEdges[id];
      ++used;
    }
    for (int &id : slice.Connects) id = remap[id];
    slice.Positions.resize(used);
    slice.UV.resize(used);
    slice.Heights.resize(used);
    slice.TetEdges.resize(used);
    extractWidthId(slice.Connects, slice.Counts, slice.TetEdges, lookup, edgeIds, debug);
    fn.create(slice.Positions, slice.Counts, slice.Connects);
    uv.swap(slice.UV);
    heights.swap(slice.Heights);
  }
}

ZSPACE_INLINE void HarmonicSlicer::extractCapMesh(VertexTag tag, zObjMesh &mesh, zPointArray &uv,
                                                  const zFloatArray &vertexHeights, zFloatArray &heights,
                                                  const WallWidthLookup &lookup, EdgeWidthIds &edgeIds, std::vector<std::string> &debug) const {
  std::unordered_map<int, int> globalToLocal;
  zPointArray positions;
  zPointArray localUV;
  zFloatArray localHeights;
  std::vector<std::pair<int, int>> localTetEdges;
  zIntArray counts;
  zIntArray connects;
  const double w = (tag == VertexTag::Top) ? 1.0 : 0.0;
  const auto faceTag = tag == VertexTag::Top ? HarmonicField::TetFaceTag::Top : HarmonicField::TetFaceTag::Bottom;
  for (size_t faceId = 0; faceId < _Tet.Faces.size(); ++faceId) {
    const auto &face = _Tet.Faces[faceId];
    if (_Tet.FaceTags[faceId] != faceTag) continue;
    const auto &tri = face.VertexIds;
    const bool valid = std::all_of(tri.begin(), tri.end(), [&](int id) {
      return id >= 0 && static_cast<size_t>(id) < _Tet.UVW.size() &&
          static_cast<size_t>(id) < _Tet.Surface.Positions.size() &&
          toEigenVec(_Tet.Surface.Positions[id]).allFinite();
    });
    // UV trace validity is independent of cap geometry. Removing a valid
    // cap face here opens holes when refinement adds a failed UV trace;
    // unrollIsoMesh already falls back to planar projection for invalid UVs.
    if (!valid) continue;

    std::array<int, 3> localTri{};
    for (int local = 0; local < 3; ++local) {
      const int global = tri[local];
      const auto found = globalToLocal.find(global);
      if (found != globalToLocal.end()) {
        localTri[local] = found->second;
        continue;
      }
      const int localId = static_cast<int>(positions.size());
      globalToLocal[global] = localId;
      positions.push_back(_Tet.Surface.Positions[global]);
      localUV.push_back(zPoint(_Tet.UVW[global].x, _Tet.UVW[global].y, w));
      localHeights.push_back(vertexHeights[global]);
      localTetEdges.emplace_back(global, global);
      localTri[local] = localId;
    }

    counts.push_back(3);
    connects.push_back(localTri[0]);
    connects.push_back(localTri[1]);
    connects.push_back(localTri[2]);
  }

  zFnMesh fn(mesh);
  fn.clear();
  uv.clear();
  heights.clear();
  edgeIds.clear();
  if (positions.empty()) return;
  extractWidthId(connects, counts, localTetEdges, lookup, edgeIds, debug);
  fn.create(positions, counts, connects);
  uv.swap(localUV);
  heights.swap(localHeights);
}

ZSPACE_INLINE void HarmonicSlicer::extractWidthId(const zIntArray &connects, const zIntArray &counts,
    const std::vector<std::pair<int, int>> &vertexTetEdges,
    const WallWidthLookup &lookup, EdgeWidthIds &edgeIds, std::vector<std::string> &debug) const {
  const std::string tag = "[extractWidthId] ";
  edgeIds.clear();
  const auto edgeWidthId = [&](const std::pair<int, int> &a, const std::pair<int, int> &b) -> int {
    std::array<int, 4> ids = {a.first, a.second, b.first, b.second};
    std::sort(ids.begin(), ids.end());
    const int distinct = static_cast<int>(std::unique(ids.begin(), ids.end()) - ids.begin());
    if (distinct == 3) {
      const auto found = lookup.Triangles.find({ids[0], ids[1], ids[2]});
      if (found != lookup.Triangles.end()) return found->second;
    } else if (distinct == 2) {
      const auto found = lookup.Edges.find(makeEdgeKey(ids[0], ids[1]));
      if (found != lookup.Edges.end()) return found->second;
    }
    return -1;
  };

  if (connects.empty()) return;
  for (const auto &source : vertexTetEdges) {
    if (source.first < 0 || source.second < 0 ||
        source.first >= static_cast<int>(_Tet.Vertices.size()) ||
        source.second >= static_cast<int>(_Tet.Vertices.size())) {
      debug.push_back(tag + "InvalidSourceIndex: an intersection references an unavailable tetrahedral vertex; source provenance is invalid.");
      return;
    }
  }
  std::unordered_map<unsigned long long, int> uses;
  size_t offset = 0;
  for (int count : counts) {
    for (int i = 0; i < count; ++i) {
      const int a = connects[offset + i], b = connects[offset + (i + 1) % count];
      const unsigned long long key = makeEdgeKey(a, b);
      if (++uses[key] > 1) continue;
      const int id = edgeWidthId(vertexTetEdges[a], vertexTetEdges[b]);
      if (id >= 0) edgeIds[key] = id;
    }
    offset += count;
  }
  // Diagnose only unmatched singleton edges; interior edges need no side tag.
  bool hasRejected = false;
  for (const auto &entry : uses)
    if (entry.second == 1 && !edgeIds.count(entry.first)) { hasRejected = true; break; }
  if (!hasRejected) return;

  std::unordered_set<std::array<int, 3>, TriangleKeyHash> surfaceFaces, markedSideFaces;
  std::unordered_set<unsigned long long> surfaceEdges, markedSideEdges;
  for (size_t faceId = 0; faceId < _Tet.Faces.size(); ++faceId) {
    const auto &face = _Tet.Faces[faceId];
    if (face.VolumeIds[1] >= 0) continue;
    auto tri = face.VertexIds;
    std::sort(tri.begin(), tri.end());
    surfaceFaces.insert(tri);
    const bool marked = _Tet.FaceTags[faceId] == HarmonicField::TetFaceTag::Side && _Tet.FaceWidthIds[faceId] >= 0;
    if (marked) markedSideFaces.insert(tri);
    for (int e = 0; e < 3; ++e) {
      const auto key = makeEdgeKey(tri[e], tri[(e+1)%3]);
      surfaceEdges.insert(key);
      if (marked) markedSideEdges.insert(key);
    }
  }
  std::unordered_set<std::string> reported;
  const auto report = [&](const std::string &reason) {
    if (reported.insert(reason).second) debug.push_back(tag + reason);
  };
  offset = 0;
  for (int count : counts) {
    for (int i = 0; i < count; ++i) {
      const int a = connects[offset+i], b = connects[offset+(i+1)%count];
      const auto key = makeEdgeKey(a,b);
      if (uses[key] != 1 || edgeIds.count(key)) continue;
      const auto &sa = vertexTetEdges[a];
      const auto &sb = vertexTetEdges[b];
      std::array<int, 4> ids = {sa.first, sa.second, sb.first, sb.second};
      std::sort(ids.begin(), ids.end());
      const int distinct = static_cast<int>(std::unique(ids.begin(), ids.end()) - ids.begin());
      if (distinct == 3) {
        const std::array<int, 3> face = {ids[0], ids[1], ids[2]};
        if (markedSideFaces.count(face))
          report("SideLookupMissing: the source face has a side width marker but was not matched by the side lookup; inspect lookup construction.");
        else if (surfaceFaces.count(face))
          report("SurfaceFaceWithoutSideMarker: the source face is on the volume surface but lacks a side width marker; it may be a cap or a lost facet marker, not necessarily a misclassified side face.");
        else
          report("SourceFaceNotOnSurface: the boundary edge is supported by a non-surface tetrahedral face; inspect missing or rejected isosurface triangles and intersection-source correspondence.");
      } else if (distinct == 2) {
        const auto sourceKey = makeEdgeKey(ids[0], ids[1]);
        if (markedSideEdges.count(sourceKey))
          report("SideEdgeLookupMissing: the source edge belongs to a marked side face but the side lookup did not match it.");
        else if (surfaceEdges.count(sourceKey))
          report("SurfaceEdgeWithoutSideMarker: the source edge lies on the volume surface but no incident surface triangle has a side width marker; check cap classification and facet-marker transfer.");
        else
          report("SourceEdgeNotOnSurface: the rejected boundary follows an internal tetrahedral edge; check missing isosurface faces or incorrect vertex provenance.");
      } else if (distinct == 4) {
        report("SourceEdgesDoNotShareFace: the endpoint source edges cannot belong to one tetrahedral face. This can expose a triangulation diagonal after a neighbouring slice triangle was discarded, or indicate incorrect source correspondence.");
      } else {
        report("CollapsedSourceSupport: both endpoints reduce to the same source vertex; check duplicate intersections and degenerate slice edges.");
      }
    }
    offset += count;
  }

}

ZSPACE_INLINE bool HarmonicSlicer::unrollIsoMesh(zObjMesh &isoMesh,
    const zPointArray &uv, zObjMesh &mesh2D, zPointArray &positions3D,
    zPointArray &unrollPositions) const {
  zFnMesh(mesh2D).clear();
  const int iterations = (std::max)(0, i_Settings._HArapIterations);
  const bool useUVStart = (i_Settings._HUnrollStart == 1);
  zFnMesh fnIso(isoMesh);
  fnIso.getVertexPositions(positions3D);
  if (positions3D.empty() || fnIso.numPolygons() == 0 || uv.size() != positions3D.size()) {
    return false;
  }
  zIntArray connects, counts;
  fnIso.getPolygonData(connects, counts);

  zPointArray startPositions;
  bool validUVStart = useUVStart;
  if (useUVStart) {
    bool allValid = true;
    startPositions.reserve(uv.size());
    for (const zPoint &p : uv) {
      if (!std::isfinite(p.x) || !std::isfinite(p.y)) {
        allValid = false;
        break;
      }
      startPositions.emplace_back(p.x * _BottomUVScale.x + _BottomUVOrigin.x,
                                  p.y * _BottomUVScale.y + _BottomUVOrigin.y, 0.0);
    }
    validUVStart = allValid;
  }
  if (!validUVStart) zTsUnroll::planarProjection(isoMesh, startPositions);

  zTsUnroll::ArapStats arapStats;
  if (!zTsUnroll::arapUnroll(isoMesh, startPositions, iterations, unrollPositions, arapStats)) {
    unrollPositions = startPositions;
  }
  zTsUnroll::alignRigid2D(unrollPositions, startPositions);

  // Use the outer loop only to orient the chart, never to select output edges.
  zIntArray loop;
  outerBoundaryLoop(connects, counts, positions3D, loop);
  const int n = static_cast<int>(loop.size());
  if (n >= 3) {
    zPointArray loop2D(n), loop3D(n);
    for (int i = 0; i < n; ++i) {
      loop2D[i] = unrollPositions[loop[i]];
      loop3D[i] = positions3D[loop[i]];
    }
    const double area2D = signedArea2D(loop2D);
    const double area3D = signedArea2D(loop3D);
    if (std::abs(area3D) > 1e-3 * std::abs(area2D) && area2D * area3D < 0.0) {
      double cx = 0.0;
      for (const zPoint &p : unrollPositions) cx += p.x;
      cx /= static_cast<double>(unrollPositions.size());
      for (zPoint &p : unrollPositions) p.x = 2.0 * cx - p.x;
    }
  }
  zFnMesh(mesh2D).create(unrollPositions, counts, connects);
  return true;
}

ZSPACE_INLINE void HarmonicSlicer::createIsoBoundaries(const zPointArray &positions3D,
    const zPointArray &positions2D, const zIntArray &connects, const zIntArray &counts,
    const EdgeWidthIds &widthIds,
    zObjGraph &boundary3D, zObjGraph &boundary2D, zIntArray &boundaryWidthIds, std::vector<std::string> &debug) {
  zFnGraph(boundary3D).clear();
  zFnGraph(boundary2D).clear();
  boundaryWidthIds.clear();
  if (positions3D.size() != positions2D.size()) {
    debug.push_back("[createIsoBoundaries] Boundary construction failed: surface and unfolded vertex arrays do not correspond.");
    return;
  }
  std::unordered_map<unsigned long long, int> edgeUse;
  size_t offset = 0;
  for (int count : counts) {
    for (int i = 0; i < count; ++i)
      ++edgeUse[makeEdgeKey(connects[offset + i], connects[offset + (i + 1) % count])];
    offset += count;
  }
  zPointArray points3D, points2D;
  zIntArray edges;
  std::unordered_map<int, int> localVertices;
  const auto addVertex = [&](int id) {
    const auto found = localVertices.find(id);
    if (found != localVertices.end()) return found->second;
    const int localId = static_cast<int>(points3D.size());
    localVertices.emplace(id, localId);
    points3D.push_back(positions3D[id]);
    points2D.push_back(positions2D[id]);
    return localId;
  };
  std::vector<int> meshBoundaryDegree(positions3D.size(), 0);
  std::vector<int> sideDegree(positions3D.size(), 0);
  std::vector<char> touchesNonManifold(positions3D.size(), 0);
  bool nonManifold = false, rejectedBoundary = false, sideOnInterior = false;
  offset = 0;
  for (int count : counts) {
    for (int i = 0; i < count; ++i) {
      const int a = connects[offset + i], b = connects[offset + (i + 1) % count];
      const unsigned long long key = makeEdgeKey(a, b);
      const auto found = widthIds.find(key);
      // Only singleton mesh edges whose source lies on an input side face.
      // Never close gaps by connecting unrelated ends or inventing width IDs.
      if (edgeUse[key] > 2) {
        nonManifold = true;
        touchesNonManifold[a] = touchesNonManifold[b] = 1;
      }
      if (edgeUse[key] != 1) {
        if (found != widthIds.end()) sideOnInterior = true;
        continue;
      }
      ++meshBoundaryDegree[a];
      ++meshBoundaryDegree[b];
      if (found == widthIds.end()) {
        rejectedBoundary = true;
        continue;
      }
      ++sideDegree[a];
      ++sideDegree[b];
      edges.push_back(addVertex(a));
      edges.push_back(addVertex(b));
      boundaryWidthIds.push_back(found->second);
    }
    offset += count;
  }
  bool filterBreak = false, meshBoundaryIrregular = false, sideBranch = false;
  bool endpointOnNonManifold = false;
  zIntArray endpoints;
  for (size_t v = 0; v < sideDegree.size(); ++v) {
    if (sideDegree[v] == 1) {
      endpoints.push_back(static_cast<int>(v));
      filterBreak |= meshBoundaryDegree[v] > sideDegree[v];
      meshBoundaryIrregular |= meshBoundaryDegree[v] != 2;
      endpointOnNonManifold |= touchesNonManifold[v] != 0;
    }
    sideBranch |= sideDegree[v] > 2;
  }
  if (filterBreak)
    debug.push_back("[createIsoBoundaries] SideFilterGap: an open endpoint touches a mesh boundary edge excluded because it has no side-face width mapping. Check source-side classification; this may also be an intentionally excluded hole edge.");
  if (meshBoundaryIrregular)
    debug.push_back("[createIsoBoundaries] IrregularMeshBoundary: an open side endpoint already has abnormal connectivity in the unfiltered mesh boundary. Inspect missing or duplicate faces and non-manifold vertex connections.");
  if (nonManifold)
    debug.push_back("[createIsoBoundaries] NonManifoldMeshEdge: a surface edge is shared by too many faces; it cannot be treated as a simple boundary edge.");
  if (endpointOnNonManifold)
    debug.push_back("[createIsoBoundaries] EndpointOnNonManifoldEdge: an open side endpoint also touches a non-manifold surface edge; check overlapping or duplicate extracted faces.");
  if (sideBranch)
    debug.push_back("[createIsoBoundaries] BranchedSideBoundary: several retained side boundary edges meet at a vertex; check touching components or non-manifold surface vertices.");
  if (!endpoints.empty()) {
    double scale = 0.0;
    for (const auto &point : positions3D)
      scale = (std::max)(scale, std::sqrt(vecDistSq(point, positions3D.front())));
    const double tolerance = (std::max)(1e-9, scale * 1e-7);
    bool nearDistinctEndpoints = false;
    for (size_t i = 0; i < endpoints.size() && !nearDistinctEndpoints; ++i)
      for (size_t j = i + 1; j < endpoints.size(); ++j)
        if (vecDistSq(positions3D[endpoints[i]], positions3D[endpoints[j]]) <= tolerance * tolerance) {
          nearDistinctEndpoints = true;
          break;
        }
    if (nearDistinctEndpoints)
      debug.push_back("[createIsoBoundaries] NearlyCoincidentEndpoints: distinct vertex indices have nearly coincident positions. Unmerged isosurface intersections or touching components may explain the gap; this is not proof that welding is valid.");
    if (sideOnInterior)
      debug.push_back("[createIsoBoundaries] SideMappedInteriorEdge: some side-mapped edges are shared by multiple faces and excluded from the boundary. Check cap or coplanar-face connectivity near the gap; this alone does not prove a classification error.");
  }
  if (edges.empty()) {
    debug.push_back(rejectedBoundary
        ? "[createIsoBoundaries] NoSideBoundary: mesh boundary edges exist but none passed the side-face mapping filter."
        : "[createIsoBoundaries] NoMeshBoundary: no singleton surface edges are available; check whether the surface is empty, closed or contains duplicate faces.");
    return;
  }
  zFnGraph(boundary3D).create(points3D, edges);
  zFnGraph(boundary2D).create(points2D, edges);
}

ZSPACE_INLINE bool HarmonicSlicer::computeContours(bool useParallel) {
  clearContours();
  if (std::find(_Tet.UVWValid.begin(), _Tet.UVWValid.end(), true) == _Tet.UVWValid.end()) {
    _Status = "computeContours failed: needs computeVolumetricParameterization() first.";
    printDebug("[computeContours] " + _Status);
    printError("[HarmonicSlicer::computeContours] Failed");
    return false;
  }

  double shortest = 0.0, longest = 0.0;
  const int interiorCount = columnLengths(shortest, longest)
      ? std::max(0, layersFromSpacing(shortest, i_Settings.LayerSpacing) - 2) : 0;
  const int numLayers = _RequestedSectionCount>0?_RequestedSectionCount:interiorCount+2;
  _RequestedSectionCount=0;
  if(i_Settings.HPrintHeight.Enabled&&i_Settings._HIsoRelaxIterations==0){
    _Status="computeContours failed: height optimization requires enabled iso relaxation.";return false;
  }
  zFloatArray vertexHeights;
  if (!computeVertexLayerHeights(numLayers, vertexHeights)) {
    _Status = "computeContours failed: harmonic gradients cannot provide finite positive layer heights.";
    printDebug("[computeContours] " + _Status);
    printError("[HarmonicSlicer::computeContours] Failed");
    return false;
  }
  o_IsoMeshes.assign(numLayers, zObjMesh());
  o_IsoMeshUV.assign(numLayers, zPointArray());
  o_VertexLayerHeights.assign(numLayers, zFloatArray());
  _HarmonicIsoEdgeWidthIds.assign(numLayers, EdgeWidthIds());
  o_IsoMeshes2D.assign(numLayers, zObjMesh());
  o_IsoMeshBoundary.assign(numLayers, zObjGraph());
  o_IsoMeshBoundary2D.assign(numLayers, zObjGraph());
  o_IsoMeshBoundaryWidthId.assign(numLayers, zIntArray());
  std::vector<char> generated(numLayers, 0);
  std::vector<std::vector<std::string>> layerDebug(numLayers);
  _IsoValues.resize(numLayers);for(int k=0;k<numLayers;++k)_IsoValues[k]=double(k)/(numLayers-1);

  WallWidthLookup lookup;
  // buildTetVolume preserves input side-face width IDs on the refined surface;
  // caps carry -1. Use tet IDs here, including surface Steiner vertices.
  for (size_t faceId = 0; faceId < _Tet.Faces.size(); ++faceId) {
    const auto &face = _Tet.Faces[faceId];
    if (_Tet.FaceTags[faceId] != HarmonicField::TetFaceTag::Side || _Tet.FaceWidthIds[faceId] < 0) continue;
    const auto &tri = face.VertexIds;
    auto key = tri;
    std::sort(key.begin(), key.end());
    const int widthId = _Tet.FaceWidthIds[faceId];
    lookup.Triangles[key] = widthId;
    for (int e = 0; e < 3; ++e)
      lookup.Edges[makeEdgeKey(tri[e], tri[(e + 1) % 3])] = widthId;
  }

  zTsDetail::ParallelErrors parallelErrors1(numLayers);
  #pragma omp parallel for schedule(dynamic) if(useParallel)
  for (int k = 0; k < numLayers; ++k) {
    try {
      if (k == 0 || k == numLayers - 1) {
        extractCapMesh(k == 0 ? VertexTag::Bottom : VertexTag::Top,
                       o_IsoMeshes[k], o_IsoMeshUV[k], vertexHeights, o_VertexLayerHeights[k],
                       lookup, _HarmonicIsoEdgeWidthIds[k], layerDebug[k]);
      } else {
        const float iso = static_cast<float>(_IsoValues[k]);
        extractIsoMesh(iso, o_IsoMeshes[k], o_IsoMeshUV[k], vertexHeights, o_VertexLayerHeights[k],
                       lookup, _HarmonicIsoEdgeWidthIds[k], layerDebug[k]);
      }
    } catch (...) { parallelErrors1.capture(k); }
  }
  parallelErrors1.rethrow();
  if(i_Settings.HPrintHeight.Enabled&&i_Settings.HPrintHeight.RedistributeIsoValues){
    std::string error;
    if(!zTsRelax::redistributeIsoValues(relaxationLayers(),i_Settings.HPrintHeight,_IsoValues,useParallel,error)){
      _Status="computeContours failed: height iso schedule: "+error;return false;
    }
    zTsDetail::ParallelErrors scheduleErrors(numLayers);
    #pragma omp parallel for schedule(dynamic) if(useParallel)
    for(int k=1;k<numLayers-1;++k){try{
      _HarmonicIsoEdgeWidthIds[k].clear();layerDebug[k].clear();
      extractIsoMesh(float(_IsoValues[k]),o_IsoMeshes[k],o_IsoMeshUV[k],vertexHeights,o_VertexLayerHeights[k],lookup,_HarmonicIsoEdgeWidthIds[k],layerDebug[k]);
    }catch(...){scheduleErrors.capture(k);}}
    scheduleErrors.rethrow();
  }
  std::string relaxationError;
  bool relaxed=relaxContours(relaxationError,useParallel);
  if(!relaxed&&i_Settings.HPrintHeight.Enabled&&i_Settings.HPrintHeight.RedistributeIsoValues){
    printDebug("[heightOptimizer] Scheduled seed rejected ("+relaxationError+"); restoring uniform iso values.");
    for(int k=1;k<numLayers-1;++k){_IsoValues[k]=double(k)/(numLayers-1);_HarmonicIsoEdgeWidthIds[k].clear();layerDebug[k].clear();
      extractIsoMesh(float(_IsoValues[k]),o_IsoMeshes[k],o_IsoMeshUV[k],vertexHeights,o_VertexLayerHeights[k],lookup,_HarmonicIsoEdgeWidthIds[k],layerDebug[k]);}
    relaxationError.clear();relaxed=relaxContours(relaxationError,useParallel);
  }
  if (!relaxed) {
    clearContours();
    _Status = "computeContours failed: iso relaxation: " + relaxationError;
    printDebug("[computeContours] " + _Status);
    return false;
  }
  if(i_Settings.HPrintHeight.Enabled&&i_Settings.HPrintHeight.OptimizeSurfaces){
    std::string error;if(!optimizeSurfaceHeightDomain(nullptr,error,useParallel)){
      _Status="computeContours failed: surface height optimization: "+error;return false;}
  }
  zTsDetail::ParallelErrors parallelUnroll(numLayers);
  #pragma omp parallel for schedule(dynamic) if(useParallel)
  for (int k = 0; k < numLayers; ++k) {
    try {
      zPointArray positions3D, positions2D;
      if (unrollIsoMesh(o_IsoMeshes[k], o_IsoMeshUV[k], o_IsoMeshes2D[k],
                         positions3D, positions2D)) {
        generated[k] = 1;
        zIntArray connects, counts;
        zFnMesh(o_IsoMeshes[k]).getPolygonData(connects, counts);
        createIsoBoundaries(positions3D, positions2D, connects, counts, _HarmonicIsoEdgeWidthIds[k],
                             o_IsoMeshBoundary[k], o_IsoMeshBoundary2D[k], o_IsoMeshBoundaryWidthId[k], layerDebug[k]);
      }

      // Both boundary graphs share connectivity. Count connected components
      // independently so a branched component cannot hide other closed loops.
      zFnGraph boundary(o_IsoMeshBoundary[k]);
      zIntArray edges;
      boundary.getEdgeData(edges);
      std::vector<std::vector<int>> adjacency(boundary.numVertices());
      for (size_t e = 0; e + 1 < edges.size(); e += 2) {
        adjacency[edges[e]].push_back(edges[e + 1]);
        adjacency[edges[e + 1]].push_back(edges[e]);
      }
      int openEnds = 0, branchVertices = 0, closedLoops = 0, usablePaths = 0;
      std::vector<char> visited(adjacency.size(), 0);
      for (size_t start = 0; start < adjacency.size(); ++start) {
        if (visited[start] || adjacency[start].empty()) continue;
        std::vector<int> pending(1, static_cast<int>(start));
        visited[start] = 1;
        int vertices = 0, ends = 0, branches = 0;
        while (!pending.empty()) {
          const int v = pending.back();
          pending.pop_back();
          ++vertices;
          if (adjacency[v].size() == 1) ++ends;
          if (adjacency[v].size() > 2) ++branches;
          for (int next : adjacency[v]) {
            if (visited[next]) continue;
            visited[next] = 1;
            pending.push_back(next);
          }
        }
        openEnds += ends;
        branchVertices += branches;
        const bool closed = branches == 0 && ends == 0 && vertices >= 3;
        if (closed) ++closedLoops;
        if (closed || (branches == 0 && ends == 2)) ++usablePaths;
      }
      const std::string tag = "[computeContours] ";
      if (openEnds > 0 || branchVertices > 0)
        layerDebug[k].push_back(tag + "OpenOrBranchedBoundary: side edges contain open ends or branches; a simple closed boundary is required for slot placement.");
      if (closedLoops > 1)
        layerDebug[k].push_back(tag + "MultipleBoundaryLoops: disjoint closed boundaries exist; slot placement selects only one loop.");
      if (usablePaths == 0)
        layerDebug[k].push_back(tag + " NoUsableBoundaryPath: no simple closed loop or open chain" +
                     (generated[k] ? std::string(".") : std::string("; no usable unfolded mesh.")));
    } catch (...) { parallelUnroll.capture(k); }
  }
  parallelUnroll.rethrow();

  // Workers only append to their own layer. Print serially in zero-based layer order.
  for (int k = 0; k < numLayers; ++k) {
    for (const std::string &message : layerDebug[k]) {
      const size_t prefixEnd = message.find(']');
      const std::string context = " layer=" + std::to_string(k);
      if (prefixEnd == std::string::npos) printDebug("[computeContours]" + context + " " + message);
      else printDebug(message.substr(0, prefixEnd + 1) + context + message.substr(prefixEnd + 1));
    }
  }

  const int generatedLayers = static_cast<int>(std::count(generated.begin(), generated.end(), char(1)));
  const bool unrolled = generatedLayers > 0;
  _Status = "computeContours: generated " + std::to_string(generatedLayers) + " layers.";
  if (!unrolled) printDebug("[computeContours] No usable unfolded contour mesh was generated.");

  if (unrolled) printSuccess("[HarmonicSlicer::computeContours] Success");
  else printError("[HarmonicSlicer::computeContours] Failed");
  return unrolled;
}

} // namespace zSpace
