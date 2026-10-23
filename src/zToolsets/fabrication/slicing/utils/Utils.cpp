// This file is part of zspace, a simple C++ collection of geometry
// data-structures & algorithms, data analysis & visualization framework.
//
// Copyright (C) 2019 ZSPACE
//
// This Source Code Form is subject to the terms of the MIT License
// If a copy of the MIT License was not distributed with this file, You can
// obtain one at https://opensource.org/licenses/MIT.
//
// Author : Weijie Tao <weijie.tao@zspace.com>
//

#ifndef EIGEN_DONT_PARALLELIZE
#define EIGEN_DONT_PARALLELIZE
#endif

#include <zspace/zToolsets/fabrication/slicing/utils/Parallel.h>
#include <zspace/zToolsets/fabrication/slicing/utils/Utils.h>

#ifndef ZSPACE_INLINE
#define ZSPACE_INLINE ZSPACE_TOOLSETS_INLINE
#endif

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <mutex>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace zSpace {
namespace zTsUtils {

//---- Thread-Safe Console Logging

namespace {
// Private: shared mutex guarding console output.
std::mutex &_getConsoleLogMutex() {
  static std::mutex mtx;
  return mtx;
}
} // anonymous namespace

//--------------------------
//---- PUBLIC FUNCTIONS
//--------------------------

ZSPACE_INLINE void printSuccess(const std::string &msg) {
  std::lock_guard<std::mutex> lock(_getConsoleLogMutex());
  std::cout << "\033[1;32m" << msg << "\033[0m\n" << std::flush;
}

ZSPACE_INLINE void printDebug(const std::string &msg) {
  std::lock_guard<std::mutex> lock(_getConsoleLogMutex());
  std::cerr << "\033[1;33m" << msg << "\033[0m\n" << std::flush;
}

ZSPACE_INLINE void printError(const std::string &msg) {
  std::lock_guard<std::mutex> lock(_getConsoleLogMutex());
  std::cerr << "\033[1;31m" << msg << "\033[0m\n" << std::flush;
}

//---- Basic Vector Utilities

ZSPACE_INLINE zVector vecSub(const zPoint &a, const zPoint &b) {
  return zVector(a.x - b.x, a.y - b.y, a.z - b.z);
}

ZSPACE_INLINE double vecDistSq(const zPoint &a, const zPoint &b) {
  const double dx = a.x - b.x;
  const double dy = a.y - b.y;
  const double dz = a.z - b.z;
  return dx * dx + dy * dy + dz * dz;
}

ZSPACE_INLINE double vecDot(const zVector &a, const zVector &b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

ZSPACE_INLINE zVector vecCross(const zVector &a, const zVector &b) {
  return zVector(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z,
                 a.x * b.y - a.y * b.x);
}

ZSPACE_INLINE unsigned long long makeEdgeKey(int u, int v) {
  const int lo = (std::min)(u, v);
  const int hi = (std::max)(u, v);
  return (static_cast<unsigned long long>(lo) << 32) |
         static_cast<unsigned int>(hi);
}

ZSPACE_INLINE Eigen::Vector3d toEigenVec(const zPoint &p) { return Eigen::Vector3d(p.x, p.y, p.z); }

ZSPACE_INLINE double triangleArea2(const zPoint &a, const zPoint &b, const zPoint &c) {
  return vecCross(vecSub(b, a), vecSub(c, a)).length();
}

ZSPACE_INLINE double meanEdgeLength(const zPointArray &positions, const std::vector<std::array<int, 3>> &triangles) {
  double lengthSum = 0.0;
  int edgeCount = 0;
  for (const auto &tri : triangles) {
    const zPoint &a = positions[tri[0]], &b = positions[tri[1]], &c = positions[tri[2]];
    if (triangleArea2(a, b, c) < 1e-8) continue; // skip degenerate triangles
    lengthSum += std::sqrt(vecDistSq(a, b)) + std::sqrt(vecDistSq(b, c)) + std::sqrt(vecDistSq(c, a));
    edgeCount += 3;
  }
  return (edgeCount > 0) ? (lengthSum / static_cast<double>(edgeCount)) : 0.0;
}

ZSPACE_INLINE bool turnAngleAt(const zPoint &prev, const zPoint &cur,
                               const zPoint &next, double &outAngleRad) {
  zVector v0 = vecSub(prev, cur);
  zVector v1 = vecSub(next, cur);
  const double len0 = v0.length();
  const double len1 = v1.length();
  if (len0 < 1e-3 || len1 < 1e-3)
    return false;

  double cosTheta = vecDot(v0, v1) / (len0 * len1);
  cosTheta = (std::max)(-1.0, (std::min)(1.0, cosTheta));
  outAngleRad = std::acos(cosTheta);
  return true;
}

ZSPACE_INLINE double signedArea2D(const zPointArray &positions) {
  double area = 0.0;
  const int n = static_cast<int>(positions.size());
  for (int i = 0; i < n; i++) {
    const zPoint &p0 = positions[i];
    const zPoint &p1 = positions[(i + 1) % n];
    area += (p0.x * p1.y - p1.x * p0.y);
  }
  return 0.5 * area;
}

ZSPACE_INLINE int outerBoundaryLoop(const zIntArray &connects, const zIntArray &counts,
                                   const zPointArray &positions, zIntArray &outLoop) {
  outLoop.clear();
  std::unordered_map<unsigned long long, int> edgeUse;
  size_t offset = 0;
  for (int count : counts) {
    for (int i = 0; i < count; ++i)
      ++edgeUse[makeEdgeKey(connects[offset + i], connects[offset + (i + 1) % count])];
    offset += count;
  }
  std::unordered_map<int, int> nextVertex;
  offset = 0;
  for (int count : counts) {
    for (int i = 0; i < count; ++i) {
      const int a = connects[offset + i], b = connects[offset + (i + 1) % count];
      if (edgeUse[makeEdgeKey(a, b)] == 1) nextVertex[a] = b;
    }
    offset += count;
  }

  int loopCount = 0;
  double bestLength = -1.0;
  std::unordered_set<int> visited;
  for (const auto &entry : nextVertex) {
    if (visited.count(entry.first)) continue;
    zIntArray loop;
    double length = 0.0;
    int current = entry.first;
    do {
      if (!visited.insert(current).second) break; // pinched / malformed loop guard
      loop.push_back(current);
      const auto found = nextVertex.find(current);
      if (found == nextVertex.end()) break;
      length += std::sqrt(vecDistSq(positions[current], positions[found->second]));
      current = found->second;
    } while (current != entry.first);
    if (loop.size() < 3) continue;
    ++loopCount;
    if (length > bestLength) {
      bestLength = length;
      outLoop = loop;
    }
  }
  return loopCount;
}

//---- HalfEdge Mesh Topology Walk Utilities

ZSPACE_INLINE zItMeshHalfEdge nextEdgeRing(zItMeshHalfEdge he) {
  return he.getNext().getNext().getSym();
}

ZSPACE_INLINE zItMeshHalfEdge nextEdgeLoop(zItMeshHalfEdge he) {
  return he.getNext().getSym().getNext();
}

ZSPACE_INLINE zItMeshHalfEdge zUpEdge(zItMeshVertex v) {
  zItMeshHalfEdgeArray hEdges;
  v.getConnectedHalfEdges(hEdges);
  if (hEdges.empty())
    return zItMeshHalfEdge();

  return *std::max_element(hEdges.begin(), hEdges.end(),
                           [](zItMeshHalfEdge a, zItMeshHalfEdge b) {
                             return a.getVector().z < b.getVector().z;
                           });
}

ZSPACE_INLINE void getFaceVerticesFromHalfedge(zItMeshHalfEdge &heStart,
                                               bool forward,
                                               zIntArray &faceVerts) {
  faceVerts.clear();
  zItMeshHalfEdge he = heStart;
  do {
    faceVerts.push_back(forward ? he.getVertex().getId()
                                : he.getStartVertex().getId());
    he = forward ? he.getNext() : he.getPrev();
  } while (he != heStart);
}

//---- Combine Multiple Graph Objects

ZSPACE_INLINE void combineGraphObjects(const zObjGraphArray &sources,
                                       zObjGraph &outGraph) {
  zPointArray positions;
  zIntArray edgeConnects;
  for (const zObjGraph &source : sources) {
    zFnGraph fnSource(const_cast<zObjGraph &>(source));
    zPointArray srcPos;
    zIntArray srcEdges;
    fnSource.getVertexPositions(srcPos);
    fnSource.getEdgeData(srcEdges);
    const int baseIndex = static_cast<int>(positions.size());
    positions.insert(positions.end(), srcPos.begin(), srcPos.end());
    for (int id : srcEdges)
      edgeConnects.push_back(baseIndex + id);
  }
  zFnGraph fnOut(outGraph);
  fnOut.clear();
  if (!positions.empty())
    fnOut.create(positions, edgeConnects);
}

//---- Input Tags

ZSPACE_INLINE bool vertexTagsFromFaceTags(zObjMesh &mesh,
                                          const zIntArray &faceTags,
                                          zIntArray &outVertexTags,
                                          std::string *outError) {
  outVertexTags.clear();
  const auto fail = [&](const std::string &reason) {
    if (outError)
      *outError = reason;
    return false;
  };
  zFnMesh fn(mesh);
  zIntArray connects, counts;
  fn.getPolygonData(connects, counts);
  if (faceTags.size() != counts.size())
    return fail(std::to_string(faceTags.size()) + " face tags for " +
                std::to_string(counts.size()) + " faces");

  const int kBottom = static_cast<int>(VertexTag::Bottom);
  const int kTop = static_cast<int>(VertexTag::Top);
  const int kSide = static_cast<int>(VertexTag::Side);
  std::vector<char> onBottom(fn.numVertices(), 0), onTop(fn.numVertices(), 0);
  int bottomFaces = 0, topFaces = 0;
  size_t offset = 0;
  for (size_t f = 0; f < counts.size(); f++) {
    const int tag = faceTags[f];
    if (tag != kBottom && tag != kTop && tag != kSide)
      return fail("face " + std::to_string(f) + " has tag " +
                  std::to_string(tag) + " (expected Bottom/Top/Side)");
    if (tag == kBottom)
      bottomFaces++;
    if (tag == kTop)
      topFaces++;
    for (int k = 0; k < counts[f]; k++) {
      const int v = connects[offset + k];
      if (tag == kBottom)
        onBottom[v] = 1;
      else if (tag == kTop)
        onTop[v] = 1;
    }
    offset += counts[f];
  }
  if (bottomFaces == 0 || topFaces == 0)
    return fail("needs at least one bottom and one top face (bottom=" +
                std::to_string(bottomFaces) + ", top=" +
                std::to_string(topFaces) + ")");

  outVertexTags.assign(fn.numVertices(), kSide);
  int conflicts = 0, firstConflict = -1;
  for (int v = 0; v < fn.numVertices(); v++) {
    if (onBottom[v] && onTop[v]) {
      if (firstConflict < 0)
        firstConflict = v;
      conflicts++;
    } else if (onBottom[v]) {
      outVertexTags[v] = kBottom;
    } else if (onTop[v]) {
      outVertexTags[v] = kTop;
    }
  }
  if (conflicts > 0) {
    outVertexTags.clear();
    return fail(std::to_string(conflicts) +
                " vertices lie on both a bottom and a top face (first: vertex " +
                std::to_string(firstConflict) +
                ") - add at least one side face between floor and ceiling");
  }
  return true;
}

ZSPACE_INLINE int widthGroupCount(const zIntArray &faceTags,
                                   const zIntArray &faceWidthIds) {
  int count = 0;
  for (size_t f = 0; f < faceTags.size(); f++) {
    if (faceTags[f] != static_cast<int>(VertexTag::Side))
      continue;
    const int id = (f < faceWidthIds.size()) ? faceWidthIds[f] : 0;
    count = (std::max)(count, id + 1);
  }
  return count;
}

//---- Layering Shared by the Slicer Operators

ZSPACE_INLINE double shortestBottomToTopLength(const zPointArray &positions,
                                               const zIntArray &edgePairs,
                                               const zIntArray &vertexTags) {
  const int n = static_cast<int>(positions.size());
  if (n == 0 || static_cast<int>(vertexTags.size()) != n)
    return -1.0;

  std::vector<std::vector<std::pair<int, double>>> adjacency(n);
  for (size_t e = 0; e + 1 < edgePairs.size(); e += 2) {
    const int a = edgePairs[e], b = edgePairs[e + 1];
    if (a < 0 || b < 0 || a >= n || b >= n || a == b)
      continue;
    const double length = std::sqrt(vecDistSq(positions[a], positions[b]));
    adjacency[a].push_back({b, length});
    adjacency[b].push_back({a, length});
  }

  typedef std::pair<double, int> QueueItem; // (distance, vertex)
  std::priority_queue<QueueItem, std::vector<QueueItem>,
                      std::greater<QueueItem>>
      queue;
  std::vector<double> distance(n, (std::numeric_limits<double>::max)());
  for (int v = 0; v < n; v++) {
    if (vertexTags[v] == static_cast<int>(VertexTag::Bottom)) {
      distance[v] = 0.0;
      queue.push({0.0, v});
    }
  }
  while (!queue.empty()) {
    const QueueItem item = queue.top();
    queue.pop();
    const double d = item.first;
    const int v = item.second;
    if (d > distance[v])
      continue;
    if (vertexTags[v] == static_cast<int>(VertexTag::Top))
      return d; // first Top vertex settled = shortest Bottom->Top path
    for (const std::pair<int, double> &next : adjacency[v]) {
      const double nd = d + next.second;
      if (nd < distance[next.first]) {
        distance[next.first] = nd;
        queue.push({nd, next.first});
      }
    }
  }
  return -1.0;
}

ZSPACE_INLINE double shortestBottomToTopLength(zObjMesh &mesh,
                                               const zIntArray &vertexTags) {
  zFnMesh fn(mesh);
  zPointArray positions;
  fn.getVertexPositions(positions);
  zIntArray edgePairs;
  for (zItMeshHalfEdge he(mesh); !he.end(); he++) {
    edgePairs.push_back(he.getStartVertex().getId());
    edgePairs.push_back(he.getVertex().getId());
  }
  return shortestBottomToTopLength(positions, edgePairs, vertexTags);
}

ZSPACE_INLINE int layersFromSpacing(double length, double spacing) {
  if (!(length > 1e-6) || !(spacing > 0.0))
    return 0;
  return (std::max)(1, static_cast<int>(std::ceil(length / spacing))) + 1;
}

//---- Single-Chain Graph Closed/Open Detection

ZSPACE_INLINE bool isClosedGraph(const zObjGraph &graph) {
  zFnGraph fn(const_cast<zObjGraph &>(graph));
  const int N = fn.numVertices();
  if (N < 3)
    return false;
  zIntArray edges;
  fn.getEdgeData(edges);
  for (size_t e = 0; e + 1 < edges.size(); e += 2) {
    if ((edges[e] == 0 && edges[e + 1] == N - 1) ||
        (edges[e + 1] == 0 && edges[e] == N - 1))
      return true;
  }
  return false;
}

//---- Universal Polyline Sampling by Arc-Length with Optional Tangent Vector

ZSPACE_INLINE bool graphPolylines(const zObjGraph &graph,
                                 std::vector<GraphPolyline> &outCurves) {
  outCurves.clear();
  zFnGraph fn(const_cast<zObjGraph &>(graph));
  const int n = fn.numVertices();
  zIntArray edges;
  fn.getEdgeData(edges);
  if (n == 0 || edges.empty() || edges.size() % 2) return false;
  std::vector<zIntArray> adjacency(n);
  for (int e = 0; e < static_cast<int>(edges.size() / 2); ++e) {
    const int a = edges[2 * e], b = edges[2 * e + 1];
    if (a < 0 || b < 0 || a >= n || b >= n || a == b) return false;
    adjacency[a].push_back(e);
    adjacency[b].push_back(e);
    if (adjacency[a].size() > 2 || adjacency[b].size() > 2) return false;
  }
  std::vector<char> used(edges.size() / 2, 0);
  const auto walk = [&](int start) {
    GraphPolyline curve;
    int current = start;
    do {
      curve.Vertices.push_back(current);
      int edge = -1;
      for (int e : adjacency[current]) if (!used[e]) { edge = e; break; }
      if (edge < 0) break;
      used[edge] = 1;
      curve.Edges.push_back(edge);
      current = edges[2 * edge] == current ? edges[2 * edge + 1] : edges[2 * edge];
      if (current == start) { curve.Closed = true; break; }
    } while (true);
    if (!curve.Edges.empty()) outCurves.push_back(std::move(curve));
  };
  for (int v = 0; v < n; ++v)
    if (adjacency[v].size() == 1 && !used[adjacency[v][0]]) walk(v);
  for (int e = 0; e < static_cast<int>(used.size()); ++e)
    if (!used[e]) walk(edges[2 * e]);
  return !outCurves.empty();
}

ZSPACE_INLINE bool graphtopolyline(const zObjGraph &graph,
                                   zPointArray &outPoints,
                                   zScalarArray &outArcLengths) {
  outPoints.clear();
  outArcLengths.clear();
  zFnGraph fn(const_cast<zObjGraph &>(graph));
  zPointArray positions;
  zIntArray edges;
  fn.getVertexPositions(positions);
  fn.getEdgeData(edges);
  const int n = static_cast<int>(positions.size());
  if (n < 2 || edges.empty() || edges.size() % 2 != 0) return false;
  for (const zPoint &p : positions)
    if (!std::isfinite(p.x) || !std::isfinite(p.y) || !std::isfinite(p.z))
      return false;

  struct Neighbor { int vertex; double length; };
  std::vector<std::array<Neighbor, 2>> adjacency(n);
  std::vector<int> degree(n, 0);
  for (size_t e = 0; e < edges.size(); e += 2) {
    const int a = edges[e], b = edges[e + 1];
    if (a < 0 || b < 0 || a >= n || b >= n || a == b ||
        degree[a] == 2 || degree[b] == 2) return false;
    for (int j = 0; j < degree[a]; ++j)
      if (adjacency[a][j].vertex == b) return false;
    const double length = std::hypot(positions[b].x - positions[a].x,
                                     positions[b].y - positions[a].y,
                                     positions[b].z - positions[a].z);
    if (!std::isfinite(length) || length <= 0.0) return false;
    adjacency[a][degree[a]++] = {b, length};
    adjacency[b][degree[b]++] = {a, length};
  }

  int endpoints = 0, start = -1;
  for (int v = 0; v < n; ++v) {
    if (degree[v] == 0) return false;
    if (degree[v] == 1) {
      ++endpoints;
      if (start < 0) start = v;
    }
  }
  if (endpoints != 0 && endpoints != 2) return false;
  const bool closed = endpoints == 0;
  if (closed) start = 0;
  const size_t edgeCount = edges.size() / 2;
  if (edgeCount != static_cast<size_t>(closed ? n : n - 1)) return false;

  zPointArray points;
  zScalarArray arc;
  points.reserve(edgeCount + 1);
  arc.reserve(edgeCount + 1);
  points.push_back(positions[start]);
  arc.push_back(0.0);
  std::vector<unsigned char> visited(n, 0);
  visited[start] = 1;
  int previous = -1, current = start, visitedCount = 1;
  double total = 0.0;
  for (size_t e = 0; e < edgeCount; ++e) {
    int next = -1;
    double length = 0.0;
    for (int j = 0; j < degree[current]; ++j) {
      const Neighbor &neighbor = adjacency[current][j];
      if (neighbor.vertex != previous) {
        next = neighbor.vertex;
        length = neighbor.length;
        break;
      }
    }
    if (next < 0) return false;
    if (visited[next]) {
      if (!closed || next != start || e + 1 != edgeCount) return false;
    } else {
      visited[next] = 1;
      ++visitedCount;
    }
    total += length;
    const zScalar s = static_cast<zScalar>(total);
    if (!std::isfinite(s) || s <= arc.back()) return false;
    points.push_back(positions[next]);
    arc.push_back(s);
    previous = current;
    current = next;
  }
  if (visitedCount != n || (closed ? current != start : degree[current] != 1))
    return false;
  outPoints.swap(points);
  outArcLengths.swap(arc);
  return true;
}

ZSPACE_INLINE bool samplePolyline(const zPoint *positions,
                                  const zScalar *scalars, size_t count,
                                  double t, zPoint &outPoint,
                                  zVector *outTangent, double minSegmentLength) {
  if (!positions || !scalars || count < 1)
    return false;

  if (count == 1 || t <= 0.0) {
    outPoint = positions[0];
    if (outTangent && count > 1) {
      *outTangent = vecSub(positions[1], positions[0]);
      if (outTangent->length() > 1e-3)
        outTangent->normalize();
    }
    return true;
  }
  if (t >= 1.0) {
    outPoint = positions[count - 1];
    if (outTangent && count > 1) {
      *outTangent = vecSub(positions[count - 1], positions[count - 2]);
      if (outTangent->length() > 1e-3)
        outTangent->normalize();
    }
    return true;
  }

  const double targetDist = t * scalars[count - 1];
  auto it = std::upper_bound(scalars, scalars + count, targetDist);
  size_t s = (it == scalars) ? 0 : std::distance(scalars, it) - 1;
  if (s + 1 >= count) {
    outPoint = positions[count - 1];
    return true;
  }

  const double segLen = scalars[s + 1] - scalars[s];
  const double localT =
      (segLen > minSegmentLength) ? (targetDist - scalars[s]) / segLen : 0.0;
  const zPoint &a = positions[s];
  const zPoint &b = positions[s + 1];
  outPoint.x = a.x + (b.x - a.x) * localT;
  outPoint.y = a.y + (b.y - a.y) * localT;
  outPoint.z = a.z + (b.z - a.z) * localT;

  if (outTangent) {
    *outTangent = vecSub(b, a);
    if (outTangent->length() > 1e-3)
      outTangent->normalize();
  }
  return true;
}

ZSPACE_INLINE bool samplePolyline(const zObjGraph &polylineGraph,
                                  const zScalarArray &scalars, double t,
                                  zPoint &outPoint, zVector *outTangent) {
  zFnGraph fnGraph(const_cast<zObjGraph &>(polylineGraph));
  zPointArray positions;
  fnGraph.getVertexPositions(positions);
  return samplePolyline(positions.data(), scalars.data(),
                        (std::min)(positions.size(), scalars.size()), t, outPoint, outTangent);
}

//---- Universal Mesh Barycentric Parameterizer (2D unrolled and 3D space meshes)

ZSPACE_INLINE bool parameterizeGraph(const zObjGraph &inGraph, zObjMesh &inMesh,
                                     ParameterGraph &outParams,
                                     const zPointArray *positionsOverride,
                                     bool allowNearest, std::string *outError,
                                     const zPointArray *diagnosticTargetPositions) {
  if (outError) outError->clear();
  const auto fail = [&](const std::string &reason) {
    if (outError) *outError = reason;
    outParams.clear();
    return false;
  };
  outParams.clear();
  zFnGraph fnGraph(const_cast<zObjGraph &>(inGraph));
  zPointArray graphPositions;
  fnGraph.getVertexPositions(graphPositions);
  fnGraph.getEdgeData(outParams.edgeConnects);
  if (graphPositions.empty() || outParams.edgeConnects.empty())
    return fail("EmptyGraph");

  zFnMesh fnMesh(inMesh);
  zPointArray meshPositions;
  if (positionsOverride) {
    if (int(positionsOverride->size()) != fnMesh.numVertices())
      return fail("PositionCountMismatch");
    meshPositions = *positionsOverride;
  } else {
    fnMesh.getVertexPositions(meshPositions);
  }
  if (meshPositions.empty() || fnMesh.numPolygons() == 0)
    return fail("EmptyChart");

  struct TriCache {
    int faceId, fanId;
    int v0, v1, v2;
    zPoint p0, p1, p2;
    Eigen::Vector3d e0, e1, normal;
    double denom;
    bool valid;
  };

  std::vector<TriCache> allTris;
  allTris.reserve(fnMesh.numPolygons() * 2);

  for (zItMeshFace f(inMesh); !f.end(); f++) {
    zIntArray faceVIds;
    f.getVertices(faceVIds);
    if (faceVIds.size() < 3)
      continue;

    for (size_t tri = 1; tri + 1 < faceVIds.size(); tri++) {
      TriCache t;
      t.faceId = f.getId();
      t.fanId = int(tri - 1);
      t.v0 = faceVIds[0];
      t.v1 = faceVIds[tri];
      t.v2 = faceVIds[tri + 1];
      t.p0 = meshPositions[t.v0];
      t.p1 = meshPositions[t.v1];
      t.p2 = meshPositions[t.v2];

      // Cross products avoid cancellation in the Gram determinant for thin triangles.
      // Compute in double precision; size alone must not invalidate a triangle.
      t.e0 = toEigenVec(t.p1) - toEigenVec(t.p0);
      t.e1 = toEigenVec(t.p2) - toEigenVec(t.p0);
      t.normal = t.e0.cross(t.e1);
      t.denom = t.normal.squaredNorm();
      t.valid = t.e0.allFinite() && t.e1.allFinite() &&
                std::isfinite(t.denom) && t.denom > 0.0;
      allTris.push_back(t);
    }
  }

  if (allTris.empty())
    return fail("NoTriangles");
  if (std::none_of(allTris.begin(), allTris.end(), [](const TriCache &t) { return t.valid; }))
    return fail("NoValidTriangles: all triangles are degenerate or non-finite");

  // Allow 0.01% of the chart diagonal, with an absolute floor of 0.1 mm.
  Eigen::Vector3d minCorner = toEigenVec(meshPositions.front());
  Eigen::Vector3d maxCorner = minCorner;
  for (const auto &position : meshPositions) {
    const Eigen::Vector3d value = toEigenVec(position);
    if (!value.allFinite()) continue;
    minCorner = minCorner.cwiseMin(value);
    maxCorner = maxCorner.cwiseMax(value);
  }
  const double chartScale = (maxCorner - minCorner).norm();
  const double projectionDistanceTolerance = (std::max)(0.1,
      std::isfinite(chartScale) ? 1e-4 * chartScale : 0.0);
  const double projectionDistanceSqTolerance = projectionDistanceTolerance * projectionDistanceTolerance;

  const size_t numGraphVerts = graphPositions.size();
  outParams.vertices.resize(numGraphVerts);

  for (size_t i = 0; i < numGraphVerts; i++) {
    const zPoint &P = graphPositions[i];
    const auto pointFailure = [&](const std::string &reason) {
      return fail("vertex " + std::to_string(i) + " (" + std::to_string(P.x) +
          ", " + std::to_string(P.y) + ", " + std::to_string(P.z) + "): " + reason);
    };
    if (!std::isfinite(P.x) || !std::isfinite(P.y) || !std::isfinite(P.z)) {
      return pointFailure("NonFinitePoint");
    }
    double bestDistSq = (std::numeric_limits<double>::max)();
    int bestV0 = allTris[0].v0, bestV1 = allTris[0].v1, bestV2 = allTris[0].v2;
    double bestU = 1.0, bestV = 0.0, bestW = 0.0;
    const TriCache *bestTriangle = nullptr;
    double bestRaw[3] = {0.0, 0.0, 0.0};
    bool foundInside = false;

    // Strict mapping searches actual containment first. Only fall back to
    // tolerance-expanded triangles if no valid containing triangle exists.
    // This avoids rejecting two tolerance candidates before a later exact hit.
    const int passCount = allowNearest ? 1 : 2;
    for (int pass = 0; pass < passCount; ++pass) {
      const double containmentTolerance = (!allowNearest && pass == 0) ? 0.0 : 1e-4;
      for (const auto &t : allTris) {
        if (!t.valid)
          continue;

        const Eigen::Vector3d v2 = toEigenVec(P) - toEigenVec(t.p0);
        const double v = v2.cross(t.e1).dot(t.normal) / t.denom;
        const double w = t.e0.cross(v2).dot(t.normal) / t.denom;
        if (!std::isfinite(v) || !std::isfinite(w)) continue;
        const double u = 1.0 - v - w;

        if (u >= -containmentTolerance && v >= -containmentTolerance &&
            w >= -containmentTolerance) {
          const double cu = (std::max)(0.0, u);
          const double cv = (std::max)(0.0, v);
          const double cw = (std::max)(0.0, w);
          const double csum = cu + cv + cw;
          const double nu = cu / csum, nv = cv / csum, nw = cw / csum;

          zPoint proj(t.p0.x * nu + t.p1.x * nv + t.p2.x * nw,
                      t.p0.y * nu + t.p1.y * nv + t.p2.y * nw,
                      t.p0.z * nu + t.p1.z * nv + t.p2.z * nw);
          const double distSq = vecDistSq(P, proj);
          if (!allowNearest && (!std::isfinite(distSq) || distSq > projectionDistanceSqTolerance)) continue;
          if (!allowNearest && distSq <= projectionDistanceSqTolerance && foundInside && bestDistSq <= projectionDistanceSqTolerance) {
            // Adjacent triangles agree on the shared edge/vertex. Different
            // positive-weight supports mean the chart overlaps at this point.
            const int ids[3] = {t.v0, t.v1, t.v2};
            const double weights[3] = {nu, nv, nw};
            const int bestIds[3] = {bestV0, bestV1, bestV2};
            const double bestWeights[3] = {bestU, bestV, bestW};
            for (int j = 0; j < 3; ++j) {
              if (weights[j] <= 1e-4) continue;
              bool shared = false;
              for (int k = 0; k < 3; ++k)
                shared |= ids[j] == bestIds[k] && bestWeights[k] > 1e-4;
              if (!shared) {
                const double minA = (std::min)(bestRaw[0], (std::min)(bestRaw[1], bestRaw[2]));
                const double minB = (std::min)(u, (std::min)(v, w));
                return pointFailure(minA > 1e-4 && minB > 1e-4
                    ? "AmbiguousChartOverlap: the point lies inside distinct triangles; the chart overlaps and the mapping is not unique."
                    : "AmbiguousChartOverlap: candidates near an edge have incompatible interpolation supports; boundary tolerance or overlapping triangles prevent a unique mapping.");
              }
            }
          }
          if (!foundInside || distSq < bestDistSq) {
            bestDistSq = distSq;
            bestV0 = t.v0;
            bestV1 = t.v1;
            bestV2 = t.v2;
            bestU = nu;
            bestV = nv;
            bestW = nw;
            bestTriangle = &t;
            bestRaw[0] = u; bestRaw[1] = v; bestRaw[2] = w;
            foundInside = true;
            if (allowNearest && distSq < 1e-2)
              break;
          }
        } else if (allowNearest && !foundInside) {
          const double cu = (std::max)(0.0, (std::min)(1.0, u));
          const double cv = (std::max)(0.0, (std::min)(1.0, v));
          const double cw = (std::max)(0.0, (std::min)(1.0, w));
          const double csum = cu + cv + cw;
          if (csum > 1e-6) {
            const double nu = cu / csum, nv = cv / csum, nw = cw / csum;
            zPoint proj(t.p0.x * nu + t.p1.x * nv + t.p2.x * nw,
                        t.p0.y * nu + t.p1.y * nv + t.p2.y * nw,
                        t.p0.z * nu + t.p1.z * nv + t.p2.z * nw);
            const double distSq = vecDistSq(P, proj);
            if (distSq < bestDistSq) {
              bestDistSq = distSq;
              bestV0 = t.v0;
              bestV1 = t.v1;
              bestV2 = t.v2;
              bestU = nu;
              bestV = nv;
              bestW = nw;
            }
          }
        }
      }

    if (foundInside) break;
    }

    if (!std::isfinite(bestDistSq) || bestDistSq == (std::numeric_limits<double>::max)() ||
        (!allowNearest && (!foundInside || bestDistSq > projectionDistanceSqTolerance))) {
      return pointFailure(!foundInside ? "NoContainingTriangle: no valid chart triangle contains the point within the allowed containment and projection tolerances" :
          "ProjectionDistanceExceeded: the point is too far from its candidate triangle.");
    }
    outParams.vertices[i] = {bestV0,
                             bestV1,
                             bestV2,
                             static_cast<float>(bestU),
                             static_cast<float>(bestV),
                             static_cast<float>(bestW)};
  }

  return true;
}

ZSPACE_INLINE bool interpolateParameterGraph(const ParameterGraph &inParams,
                                             const zPointArray &vertexValues,
                                             zPointArray &outValues) {
  outValues.assign(inParams.vertices.size(), zPoint());
  const int numValues = int(vertexValues.size());
  for (size_t i = 0; i < inParams.vertices.size(); i++) {
    const auto &bv = inParams.vertices[i];
    if (bv.v0 < 0 || bv.v0 >= numValues || bv.v1 < 0 || bv.v1 >= numValues ||
        bv.v2 < 0 || bv.v2 >= numValues)
      return false;
    const zPoint &p0 = vertexValues[bv.v0];
    const zPoint &p1 = vertexValues[bv.v1];
    const zPoint &p2 = vertexValues[bv.v2];
    outValues[i] = zPoint(p0.x * bv.u + p1.x * bv.v + p2.x * bv.w,
                          p0.y * bv.u + p1.y * bv.v + p2.y * bv.w,
                          p0.z * bv.u + p1.z * bv.v + p2.z * bv.w);
  }
  return true;
}

//---- Universal Barycentric Mesh Graph Reconstructor

ZSPACE_INLINE bool realizeParameterGraph(const ParameterGraph &inParams,
                                         zObjMesh &targetMesh,
                                         zObjGraph &outGraph,
                                         zVectorArray *outNormals,
                                         const zFloatArray *vertexHeights,
                                         zFloatArray *outHeights,
                                         const zVectorArray *vertexNormals) {
  if (inParams.empty())
    return false;

  zFnMesh fnTarget(targetMesh);
  zPointArray targetPositions;
  fnTarget.getVertexPositions(targetPositions);
  if (targetPositions.empty())
    return false;

  const size_t numVerts = inParams.vertices.size();
  zPointArray newPositions(numVerts);
  if (outNormals)
    outNormals->assign(numVerts, zVector(0, 0, 1));
  if (outHeights)
    outHeights->assign(numVerts, 0.0f);

  const size_t numMeshVerts = targetPositions.size();
  const bool smoothNormals = vertexNormals && !vertexNormals->empty();
  if (smoothNormals && vertexNormals->size() != numMeshVerts) return false;
  const bool wantHeights = vertexHeights && outHeights;

  for (size_t i = 0; i < numVerts; i++) {
    const auto &bv = inParams.vertices[i];
    if (bv.v0 < 0 || bv.v0 >= (int)numMeshVerts || bv.v1 < 0 ||
        bv.v1 >= (int)numMeshVerts || bv.v2 < 0 || bv.v2 >= (int)numMeshVerts) {
      return false;
    }

    const zPoint &p0 = targetPositions[bv.v0];
    const zPoint &p1 = targetPositions[bv.v1];
    const zPoint &p2 = targetPositions[bv.v2];

    newPositions[i] = zPoint(p0.x * bv.u + p1.x * bv.v + p2.x * bv.w,
                             p0.y * bv.u + p1.y * bv.v + p2.y * bv.w,
                             p0.z * bv.u + p1.z * bv.v + p2.z * bv.w);
    if (!std::isfinite(newPositions[i].x) || !std::isfinite(newPositions[i].y) ||
        !std::isfinite(newPositions[i].z)) return false;

    if (wantHeights) {
      float h0 =
          (bv.v0 < (int)vertexHeights->size()) ? (*vertexHeights)[bv.v0] : 0.0f;
      float h1 =
          (bv.v1 < (int)vertexHeights->size()) ? (*vertexHeights)[bv.v1] : 0.0f;
      float h2 =
          (bv.v2 < (int)vertexHeights->size()) ? (*vertexHeights)[bv.v2] : 0.0f;
      (*outHeights)[i] = bv.u * h0 + bv.v * h1 + bv.w * h2;
    }

    if (outNormals) {
      zVector e0 = vecSub(p1, p0);
      zVector e1 = vecSub(p2, p0);
      zVector norm = vecCross(e0, e1);
      if (smoothNormals) {
        const zVector &n0 = (*vertexNormals)[bv.v0], &n1 = (*vertexNormals)[bv.v1],
                      &n2 = (*vertexNormals)[bv.v2];
        zVector interpolated(n0.x*bv.u+n1.x*bv.v+n2.x*bv.w,
                                  n0.y*bv.u+n1.y*bv.v+n2.y*bv.w,
                                  n0.z*bv.u+n1.z*bv.v+n2.z*bv.w);
        if (!std::isfinite(interpolated.x) || !std::isfinite(interpolated.y) ||
            !std::isfinite(interpolated.z) || interpolated.length() <= 1e-6) return false;
        norm = interpolated;
      }
      if (norm.length() > 1e-6)
        norm.normalize();
      else
        norm = zVector(0, 0, 1);
      (*outNormals)[i] = norm;
    }
  }

  zFnGraph fnOut(outGraph);
  fnOut.clear();
  fnOut.create(newPositions, const_cast<zIntArray &>(inParams.edgeConnects));
  return true;
}

//---- Barycentric Mesh-to-Mesh Graph Projection

ZSPACE_INLINE bool remapGraph(zObjGraph &graph, zObjMesh &inMesh2D,
                              zObjMesh &projectionMesh3D,
                              zVectorArray *mappedNormals,
                              const zFloatArray *vertexHeights,
                              zFloatArray *mappedHeights, std::string *outError,
                              const zVectorArray *vertexNormals) {
  if (outError) outError->clear();
  ParameterGraph tempParams;
  zIntArray sourceConnects, sourceCounts, targetConnects, targetCounts;
  zFnMesh(inMesh2D).getPolygonData(sourceConnects, sourceCounts);
  zFnMesh(projectionMesh3D).getPolygonData(targetConnects, targetCounts);
  if (zFnMesh(inMesh2D).numVertices() != zFnMesh(projectionMesh3D).numVertices() ||
      sourceConnects != targetConnects || sourceCounts != targetCounts) {
    if (outError) *outError = "TopologyMismatch: 2D/3D vertex counts or face connectivity differ";
    return false;
  }
  zPointArray diagnosticPositions;
  if (outError) zFnMesh(projectionMesh3D).getVertexPositions(diagnosticPositions);
  if (!parameterizeGraph(graph, inMesh2D, tempParams, nullptr, false, outError,
                          outError ? &diagnosticPositions : nullptr))
    return false;
  if (!realizeParameterGraph(tempParams, projectionMesh3D, graph,
                               mappedNormals, vertexHeights, mappedHeights, vertexNormals)) {
    if (outError) *outError = "RealizationFailed: empty target, invalid vertex index or non-finite 3D position";
    return false;
  }
  return true;
}

//---- Barycentric Point Transfer

ZSPACE_INLINE bool mapPointToMeshes(const zPoint &point, zObjMesh &fromMesh,
                                    const std::vector<const zObjMesh *> &toMeshes,
                                    zPointArray &outPoints,
                                    std::vector<char> &outValid,
                                    double *outDistance, bool useParallel) {
  outPoints.assign(toMeshes.size(), zPoint());
  outValid.assign(toMeshes.size(), 0);

  // parameterizeGraph() takes a graph: a single point as a zero-length edge.
  zPointArray pts = {point, point};
  zIntArray edges = {0, 1};
  zObjGraph pointGraph;
  zFnGraph(pointGraph).create(pts, edges);
  ParameterGraph params;
  if (!parameterizeGraph(pointGraph, fromMesh, params))
    return false;

  if (outDistance) {
    zObjGraph located;
    if (!realizeParameterGraph(params, fromMesh, located))
      return false;
    zPointArray locatedPts;
    zFnGraph(located).getVertexPositions(locatedPts);
    *outDistance = locatedPts.empty()
                       ? (std::numeric_limits<double>::max)()
                       : std::sqrt(vecDistSq(point, locatedPts[0]));
  }

  // Repeated target pointers are allowed. Resolve core's lazy storage once
  // on the calling thread, then workers only read plain coordinate arrays.
  std::vector<zPointArray> targetPositions(toMeshes.size());
  for (size_t i = 0; i < toMeshes.size(); ++i)
    if (toMeshes[i]) zFnMesh(const_cast<zObjMesh &>(*toMeshes[i])).getVertexPositions(targetPositions[i]);
  zTsDetail::ParallelErrors parallelErrors1(int(toMeshes.size()));
#pragma omp parallel for schedule(dynamic) if(useParallel)
  for (int i = 0; i < int(toMeshes.size()); i++) {
    try {
      if (!toMeshes[i])
        continue;
      zPointArray realizedPts;
      if (!interpolateParameterGraph(params, targetPositions[i], realizedPts))
        continue;
      if (realizedPts.empty())
        continue;
      if (!std::all_of(realizedPts.begin(), realizedPts.end(), [](const zPoint &p) {
            return std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z);
          })) continue;
      outPoints[i] = realizedPts[0];
      outValid[i] = 1;
    } catch (...) { parallelErrors1.capture(i); }
  }
  parallelErrors1.rethrow();
  return true;
}

ZSPACE_INLINE bool mapPoint(const zPoint &point, zObjMesh &fromMesh,
                            zObjMesh &toMesh, zPoint &outPoint,
                            double *outDistance) {
  zPointArray points;
  std::vector<char> valid;
  if (!mapPointToMeshes(point, fromMesh, {&toMesh}, points, valid,
                        outDistance, false) ||
      !valid[0])
    return false;
  outPoint = points[0];
  return true;
}

//---- Toolpath Pipe Mesh Construction

ZSPACE_INLINE void buildPathMesh(const PrintPath &path,
                                 float defaultLayerHeight,
                                 zObjMesh &outMesh) {
  zFnMesh fnOut(outMesh);
  fnOut.clear();

  zFnGraph fnContour(const_cast<zObjGraph &>(path.graph));
  zPointArray points;
  fnContour.getVertexPositions(points);
  const int N = int(points.size());
  if (N < 2)
    return;

  std::vector<GraphPolyline> curves;
  if (!graphPolylines(path.graph, curves)) return;
  zIntArray previous(N, -1), nextVertex(N, -1), pathEdges;
  for (const GraphPolyline &curve : curves) {
    const int count = static_cast<int>(curve.Vertices.size());
    for (int j = 0; j < count; ++j) {
      const int v = curve.Vertices[j];
      if (j > 0 || curve.Closed) previous[v] = curve.Vertices[(j + count - 1) % count];
      if (j + 1 < count || curve.Closed) {
        nextVertex[v] = curve.Vertices[(j + 1) % count];
        pathEdges.push_back(v);
        pathEdges.push_back(nextVertex[v]);
      }
    }
  }

  const zFloatArray *layerWeights = &path.widths;
  const zFloatArray *layerHeights = &path.heights;
  const zVectorArray *layerNormals = &path.normals;

  zPointArray meshPositions;
  zIntArray polyCounts;
  zIntArray polyConnects;
  meshPositions.reserve(N * 4);

  for (int i = 0; i < N; i++) {
    zPoint p = points[i];

    // 1. Tangent direction
    zVector pathDir;
    const int prev = previous[i] >= 0 ? previous[i] : i;
    const int next = nextVertex[i] >= 0 ? nextVertex[i] : i;
    pathDir = vecSub(points[next], points[prev]);
    if (pathDir.length() < 1e-3)
      pathDir = zVector(1, 0, 0);
    pathDir.normalize();

    // 2. 3D surface normal
    zVector norm = (layerNormals && i < int(layerNormals->size()))
                       ? (*layerNormals)[i]
                       : zVector(0, 0, 1);
    if (norm.length() < 1e-6)
      norm = zVector(0, 0, 1);
    norm.normalize();

    // 3. In-plane binormal perpendicular to tangent
    zVector binormal = norm ^ pathDir;
    if (binormal.length() < 1e-6) {
      zVector temp(0, 0, 1);
      if (std::abs(pathDir.z) > 0.9)
        temp = zVector(0, 1, 0);
      binormal = temp ^ pathDir;
    }
    binormal.normalize();

    // 4. Per-vertex width and height
    float w = (layerWeights && i < int(layerWeights->size()))
                  ? (*layerWeights)[i]
                  : 0.20f;
    if (w <= 0.1f)
      w = 0.20f;
    float h = (layerHeights && i < int(layerHeights->size()))
                  ? (*layerHeights)[i]
                  : defaultLayerHeight;
    if (h <= 0.1f)
      h = defaultLayerHeight;

    // 5. Generate 4 corner vertices for rectangular cross section
    const float halfW = w * 0.5f;
    zPoint c0(p.x + binormal.x * halfW, p.y + binormal.y * halfW,
              p.z + binormal.z * halfW); // Top-Right
    zPoint c1(p.x - binormal.x * halfW, p.y - binormal.y * halfW,
              p.z - binormal.z * halfW); // Top-Left
    zPoint c2(c1.x - norm.x * h, c1.y - norm.y * h,
              c1.z - norm.z * h); // Bottom-Left
    zPoint c3(c0.x - norm.x * h, c0.y - norm.y * h,
              c0.z - norm.z * h); // Bottom-Right

    meshPositions.push_back(c0);
    meshPositions.push_back(c1);
    meshPositions.push_back(c2);
    meshPositions.push_back(c3);
  }

  // 6. Construct 4 quads per segment along the toolpath tube
  for (size_t e = 0; e + 1 < pathEdges.size(); e += 2) {
    const int i = pathEdges[e], next = pathEdges[e + 1];
    // Top face
    polyCounts.push_back(4);
    polyConnects.push_back(i * 4 + 0);
    polyConnects.push_back(next * 4 + 0);
    polyConnects.push_back(next * 4 + 1);
    polyConnects.push_back(i * 4 + 1);

    // Left face
    polyCounts.push_back(4);
    polyConnects.push_back(i * 4 + 1);
    polyConnects.push_back(next * 4 + 1);
    polyConnects.push_back(next * 4 + 2);
    polyConnects.push_back(i * 4 + 2);

    // Bottom face
    polyCounts.push_back(4);
    polyConnects.push_back(i * 4 + 2);
    polyConnects.push_back(next * 4 + 2);
    polyConnects.push_back(next * 4 + 3);
    polyConnects.push_back(i * 4 + 3);

    // Right face
    polyCounts.push_back(4);
    polyConnects.push_back(i * 4 + 3);
    polyConnects.push_back(next * 4 + 3);
    polyConnects.push_back(next * 4 + 0);
    polyConnects.push_back(i * 4 + 0);
  }

  // Cap ends for open curves
  for (const GraphPolyline &curve : curves) {
    if (curve.Closed) continue;
    const int first = curve.Vertices.front() * 4, last = curve.Vertices.back() * 4;
    polyCounts.push_back(4);
    polyConnects.insert(polyConnects.end(), {first, first + 3, first + 2, first + 1});
    polyCounts.push_back(4);
    polyConnects.insert(polyConnects.end(), {last, last + 1, last + 2, last + 3});
  }

  fnOut.create(meshPositions, polyCounts, polyConnects);
}

//---- Generic Mesh/Graph Transform & Segmentation

ZSPACE_INLINE void translateMesh(zObjMesh &mesh, const zVector &offset) {
  zFnMesh fn(mesh);
  if (fn.numVertices() <= 0)
    return;
  zPointArray pts;
  fn.getVertexPositions(pts);
  for (auto &p : pts) {
    p.x += offset.x;
    p.y += offset.y;
    p.z += offset.z;
  }
  fn.setVertexPositions(pts);
}

ZSPACE_INLINE void translateGraph(zObjGraph &graph, const zVector &offset) {
  zFnGraph fn(graph);
  if (fn.numVertices() <= 0)
    return;
  zPointArray pts;
  fn.getVertexPositions(pts);
  for (auto &p : pts) {
    p.x += offset.x;
    p.y += offset.y;
    p.z += offset.z;
  }
  fn.setVertexPositions(pts);
}

ZSPACE_INLINE bool segmentGraph(const zObjGraph &sourceBoundaryGraph,
                                const std::vector<int> &splitIndices,
                                std::vector<zObjGraph> &outBoundarySegments) {
  outBoundarySegments.clear();
  zFnGraph fnSource(const_cast<zObjGraph &>(sourceBoundaryGraph));
  zPointArray positions;
  fnSource.getVertexPositions(positions);

  const int N = int(positions.size());
  if (N < 2)
    return false;

  const int numSegs = int(splitIndices.size());
  if (numSegs == 0) {
    outBoundarySegments.push_back(sourceBoundaryGraph);
    return true;
  }

  outBoundarySegments.reserve(numSegs);

  for (int seg = 0; seg < numSegs; seg++) {
    const int startIdx = splitIndices[seg];
    const int endIdx =
        (seg + 1 < numSegs) ? splitIndices[seg + 1] : splitIndices[0];

    zPointArray outPositions;
    zIntArray outEdges;

    if (startIdx <= endIdx) {
      for (int k = startIdx; k <= endIdx; k++) {
        outPositions.push_back(positions[k]);
      }
    } else {
      for (int k = startIdx; k < N; k++) {
        outPositions.push_back(positions[k]);
      }
      for (int k = 0; k <= endIdx; k++) {
        outPositions.push_back(positions[k]);
      }
    }

    if (outPositions.size() >= 2) {
      for (size_t e = 0; e + 1 < outPositions.size(); e++) {
        outEdges.push_back(int(e));
        outEdges.push_back(int(e + 1));
      }

      zObjGraph segGraph;
      zFnGraph fnSeg(segGraph);
      fnSeg.create(outPositions, outEdges);
      outBoundarySegments.push_back(segGraph);
    }
  }

  return !outBoundarySegments.empty();
}

} // namespace zTsUtils
} // namespace zSpace
