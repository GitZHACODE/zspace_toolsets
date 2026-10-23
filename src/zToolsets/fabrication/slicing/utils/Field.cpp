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

#include <zspace/zToolsets/fabrication/slicing/utils/Utils.h>

#ifndef ZSPACE_INLINE
#define ZSPACE_INLINE ZSPACE_TOOLSETS_INLINE
#endif

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace zSpace {
namespace zTsField {
ZSPACE_INLINE bool appendBracingSegments(const zObjGraph &graph, const SlicerSettings &settings,
    zScalar trim, std::vector<FieldSegment> &segments, int layer) {
  std::vector<zTsUtils::GraphPolyline> curves;
  if (!zTsUtils::graphPolylines(graph, curves) || (!settings.BracingCurveWidths.empty() && curves.size() != settings.BracingCurveWidths.size())) {
    zTsUtils::printDebug("[appendBracingSegments] layer=" + std::to_string(layer) +
        " Bracing topology does not match per-curve widths; field generation stopped for this layer.");
    return false;
  }
  std::sort(curves.begin(), curves.end(), [](const auto &a, const auto &b) {
    return *std::min_element(a.Vertices.begin(), a.Vertices.end()) <
           *std::min_element(b.Vertices.begin(), b.Vertices.end());
  });
  zPointArray points;
  zFnGraph(const_cast<zObjGraph &>(graph)).getVertexPositions(points);
  for (size_t c = 0; c < curves.size(); ++c) {
    zPointArray local; zIntArray edges;
    for (int id : curves[c].Vertices) local.push_back(points[id]);
    for (int i = 1; i < static_cast<int>(local.size()); ++i) { edges.push_back(i-1); edges.push_back(i); }
    if (curves[c].Closed) { edges.push_back(static_cast<int>(local.size())-1); edges.push_back(0); }
    zObjGraph rib; zFnGraph(rib).create(local, edges);
    const double width = settings.BracingCurveWidths.empty() ? settings.BracingWidth : settings.BracingCurveWidths[c];
    double length = 0.0;
    for (size_t e = 0; e < edges.size(); e += 2) {
      const auto &a = local[edges[e]], &b = local[edges[e + 1]];
      length += std::hypot(double(b.x) - a.x, double(b.y) - a.y, double(b.z) - a.z);
    }
    if (!std::isfinite(trim) || trim < 0 || trim > 1 ||
        !std::isfinite(settings.SlotOffset) || settings.SlotOffset < 0 ||
        !std::isfinite(width) || width <= 0 || !std::isfinite(length) || length <= 0) {
      zTsUtils::printDebug("[appendBracingSegments] layer=" + std::to_string(layer) +
          " Invalid BracingSlot, SlotOffset, width or curve length.");
      return false;
    }
    double ratio = double(trim) + (layer % 2 == 0 ? 1.0 : -1.0) * settings.SlotOffset * width / length;
    ratio = curves[c].Closed ? ratio - std::floor(ratio) : std::clamp(ratio, 0.0, 1.0);
    segments.emplace_back(rib, 1, static_cast<zScalar>(width), 0.0, static_cast<zScalar>(ratio));
  }
  return true;
}


using namespace zTsUtils;

//--------------------------
//---- PRIVATE HELPERS
//--------------------------

namespace {

// ---- extractPrintPaths helpers ----

// Raw marching-triangles isocontour extraction: zTsField::extractRawIsoContour.

// Spur-pruning + gap-welding + chain/loop extraction: splits the raw
// iso-contour, walks it into a single oriented, inner-before-outer closed
// ring (plus bracing open chains), with a separate PrintPath per curve.
void extractTrackPaths(const zPointArray &rawPositions,
                       const zIntArray &rawEdges, const zFloatArray &rawWidths,
                       const zFloatArray &rawTrackIds, double weldTolerance,
                       double spurLengthThreshold,
                       std::vector<PrintPath> &outPaths) {
  outPaths.clear();
  const int numV = static_cast<int>(rawPositions.size());
  if (numV < 2 || rawEdges.size() < 2)
    return;

  const bool hasTrackIds = rawTrackIds.size() == static_cast<size_t>(numV);
  auto isOuterVertex = [&](int id) -> bool {
    return hasTrackIds && id < int(rawTrackIds.size()) &&
           rawTrackIds[id] >= 0.5f;
  };

  // Contract short edges before building adjacency: every incident edge must
  // use the same representative, rather than leaving a gap by dropping an edge.
  std::vector<int> parent(numV);
  for (int i = 0; i < numV; ++i) parent[i] = i;
  const auto root = [&](int id) {
    while (parent[id] != id) {
      parent[id] = parent[parent[id]];
      id = parent[id];
    }
    return id;
  };
  const double mergeToleranceSq = weldTolerance * weldTolerance;
  for (size_t e = 0; e + 1 < rawEdges.size(); e += 2) {
    const int a = rawEdges[e], b = rawEdges[e + 1];
    if (a < 0 || b < 0 || a >= numV || b >= numV) continue;
    if (vecDistSq(rawPositions[a], rawPositions[b]) <= mergeToleranceSq) {
      const int ra = root(a), rb = root(b);
      // Deterministic representative retains its position, width and track ID.
      if (ra != rb) parent[(std::max)(ra, rb)] = (std::min)(ra, rb);
    }
  }

  std::vector<std::vector<int>> adj(numV);
  std::unordered_set<unsigned long long> edgeSet;
  for (size_t e = 0; e + 1 < rawEdges.size(); e += 2) {
    int a = rawEdges[e], b = rawEdges[e + 1];
    if (a < 0 || b < 0 || a >= numV || b >= numV || a == b)
      continue;
    a = root(a);
    b = root(b);
    if (a == b) continue;
    const unsigned long long key = makeEdgeKey(a, b);
    if (!edgeSet.insert(key).second)
      continue;
    adj[a].push_back(b);
    adj[b].push_back(a);
  }

  std::vector<std::unordered_set<int>> activeAdj(numV);
  for (int i = 0; i < numV; i++)
    for (int nb : adj[i])
      activeAdj[i].insert(nb);

  bool prunedAny = true;
  while (prunedAny) {
    prunedAny = false;
    for (int i = 0; i < numV; i++) {
      if (activeAdj[i].size() >= 3) {
        std::vector<int> neighbors(activeAdj[i].begin(), activeAdj[i].end());
        for (int startNb : neighbors) {
          std::vector<int> branchPath = {i, startNb};
          int prev = i;
          int curr = startNb;
          bool isDeadEnd = false;
          double branchLen =
              vecSub(rawPositions[curr], rawPositions[prev]).length();

          while (true) {
            if (activeAdj[curr].size() == 1) {
              isDeadEnd = true;
              break;
            }
            if (activeAdj[curr].size() != 2)
              break;
            int next = -1;
            for (int nb : activeAdj[curr]) {
              if (nb != prev) {
                next = nb;
                break;
              }
            }
            if (next < 0 || next == i)
              break;
            branchLen +=
                vecSub(rawPositions[next], rawPositions[curr]).length();
            prev = curr;
            curr = next;
            branchPath.push_back(curr);
          }

          if (isDeadEnd && branchLen < spurLengthThreshold) {
            for (size_t k = 0; k + 1 < branchPath.size(); k++) {
              int u = branchPath[k];
              int v = branchPath[k + 1];
              activeAdj[u].erase(v);
              activeAdj[v].erase(u);
            }
            prunedAny = true;
          }
        }
      }
    }
  }

  std::vector<std::vector<int>> openChains;
  std::vector<std::vector<int>> closedLoops;
  std::unordered_set<unsigned long long> visitedEdges;

  auto walkChain = [&](int start,
                       int next) -> std::pair<std::vector<int>, bool> {
    std::vector<int> seq = {start};
    int prev = start, curr = next;
    visitedEdges.insert(makeEdgeKey(prev, curr));

    while (true) {
      if (curr == start)
        return {seq, true};
      seq.push_back(curr);
      if (activeAdj[curr].size() != 2)
        break;

      int nxt = -1;
      for (int nb : activeAdj[curr]) {
        if (nb != prev) {
          nxt = nb;
          break;
        }
      }
      if (nxt < 0 || visitedEdges.count(makeEdgeKey(curr, nxt)))
        break;
      visitedEdges.insert(makeEdgeKey(curr, nxt));
      prev = curr;
      curr = nxt;
    }
    return {seq, false};
  };

  for (int i = 0; i < numV; i++) {
    if (activeAdj[i].size() == 1) {
      int next = *activeAdj[i].begin();
      if (visitedEdges.count(makeEdgeKey(i, next)))
        continue;
      openChains.push_back(walkChain(i, next).first);
    }
  }

  for (int i = 0; i < numV; i++) {
    if (activeAdj[i].size() == 2) {
      for (int nb : activeAdj[i]) {
        if (visitedEdges.count(makeEdgeKey(i, nb)))
          continue;

        auto [loop, closed] = walkChain(i, nb);
        if (closed && loop.size() >= 3) {
          closedLoops.push_back(std::move(loop));
        } else if (!loop.empty()) {
          openChains.push_back(std::move(loop));
        }
      }
    }
  }

  // 4. Self-close any open chain whose two ends are within weldTolerance -
  // safe unconditionally now that nothing is track-filtered, so the
  // boundary ring can never be mistakenly split by a filtering artifact
  // that this step would then need to leave alone.
  const double weldTolSq = weldTolerance * weldTolerance;
  for (auto it = openChains.begin(); it != openChains.end();) {
    if (it->size() >= 3 && vecDistSq(rawPositions[it->front()],
                                     rawPositions[it->back()]) <= weldTolSq) {
      closedLoops.push_back(*it);
      it = openChains.erase(it);
      continue;
    }
    ++it;
  }

  // 5. Weld remaining open-chain fragments (e.g. bracing ribs split at
  // triangle boundaries) end-to-end by nearest-endpoint matching.
  bool mergedAny = true;
  while (mergedAny && openChains.size() > 1) {
    mergedAny = false;
    size_t bestI = 0, bestJ = 0;
    int bestMode = -1;
    double bestDistSq = weldTolSq;

    for (size_t i = 0; i < openChains.size(); i++) {
      for (size_t j = i + 1; j < openChains.size(); j++) {
        const auto &c1 = openChains[i];
        const auto &c2 = openChains[j];

        double d0 =
            vecDistSq(rawPositions[c1.back()], rawPositions[c2.front()]);
        if (d0 < bestDistSq) {
          bestDistSq = d0;
          bestI = i;
          bestJ = j;
          bestMode = 0;
        }
        double d1 = vecDistSq(rawPositions[c1.back()], rawPositions[c2.back()]);
        if (d1 < bestDistSq) {
          bestDistSq = d1;
          bestI = i;
          bestJ = j;
          bestMode = 1;
        }
        double d2 =
            vecDistSq(rawPositions[c1.front()], rawPositions[c2.back()]);
        if (d2 < bestDistSq) {
          bestDistSq = d2;
          bestI = i;
          bestJ = j;
          bestMode = 2;
        }
        double d3 =
            vecDistSq(rawPositions[c1.front()], rawPositions[c2.front()]);
        if (d3 < bestDistSq) {
          bestDistSq = d3;
          bestI = i;
          bestJ = j;
          bestMode = 3;
        }
      }
    }

    if (bestMode >= 0) {
      std::vector<int> &c1 = openChains[bestI];
      std::vector<int> c2 = openChains[bestJ];

      if (bestMode == 0) {
        c1.insert(c1.end(), c2.begin(), c2.end());
      } else if (bestMode == 1) {
        std::reverse(c2.begin(), c2.end());
        c1.insert(c1.end(), c2.begin(), c2.end());
      } else if (bestMode == 2) {
        c2.insert(c2.end(), c1.begin(), c1.end());
        c1 = c2;
      } else if (bestMode == 3) {
        std::reverse(c2.begin(), c2.end());
        c2.insert(c2.end(), c1.begin(), c1.end());
        c1 = c2;
      }

      openChains.erase(openChains.begin() + bestJ);
      mergedAny = true;

      if (c1.size() >= 3 && vecDistSq(rawPositions[c1.front()],
                                      rawPositions[c1.back()]) <= weldTolSq) {
        closedLoops.push_back(c1);
        openChains.erase(openChains.begin() + bestI);
      }
    }
  }

  if (closedLoops.empty() && openChains.empty())
    return;

  // 6. Orient the longest closed loop (the boundary ring) CCW
  if (!closedLoops.empty()) {
    size_t mainIdx = 0;
    for (size_t i = 1; i < closedLoops.size(); i++) {
      if (closedLoops[i].size() > closedLoops[mainIdx].size())
        mainIdx = i;
    }
    std::vector<int> &mainRing = closedLoops[mainIdx];
    const int m = static_cast<int>(mainRing.size());

    zPointArray ringPositions(m);
    for (int i = 0; i < m; i++)
      ringPositions[i] = rawPositions[mainRing[i]];
    if (signedArea2D(ringPositions) < 0.0)
      std::reverse(mainRing.begin(), mainRing.end());

    if (hasTrackIds) {
      bool hasInner = false, hasOuter = false;
      for (int id : mainRing) {
        if (isOuterVertex(id))
          hasOuter = true;
        else
          hasInner = true;
      }
      if (hasInner && hasOuter) {
        for (int i = 0; i < m; i++) {
          const int prev = (i - 1 + m) % m;
          if (!isOuterVertex(mainRing[i]) && isOuterVertex(mainRing[prev])) {
            std::rotate(mainRing.begin(), mainRing.begin() + i, mainRing.end());
            break;
          }
        }
      }
    }
    if (mainIdx != 0) std::swap(closedLoops[0], closedLoops[mainIdx]);
  }

  // Each connected curve is a separate PrintPath. Array adjacency must not
  // create a bridge between distinct loops or chains during resampling.
  auto appendSeq = [&](const std::vector<int> &seq, bool isLoop) {
    if (seq.size() < 2) return;
    zPointArray positions;
    zIntArray edges;
    PrintPath path;
    for (int id : seq) {
      positions.push_back(rawPositions[id]);
      path.widths.push_back(id < int(rawWidths.size()) ? rawWidths[id]
                                                        : 0.04f);
      path.trackIds.push_back(isOuterVertex(id) ? 1.0f : 0.0f);
    }
    for (size_t i = 0; i + 1 < seq.size(); i++) {
      edges.push_back(int(i));
      edges.push_back(int(i) + 1);
    }
    if (isLoop && seq.size() >= 3) {
      edges.push_back(int(seq.size()) - 1);
      edges.push_back(0);
    }
    zFnGraph(path.graph).create(positions, edges);
    outPaths.push_back(std::move(path));
  };

  for (auto &loop : closedLoops)
    appendSeq(loop, true);
  for (auto &chain : openChains)
    appendSeq(chain, false);

}


struct FastSegment {
  float x0, y0;
  float x1, y1; // Preserve the directed endpoints for winding tests.
  float dx, dy;
  float lenSqInv;
  float minX, maxX;
  float minY, maxY;
};

// Shared directed-edge preparation for distance and winding queries.
// Short edges remain available to winding; a zero inverse marks an edge
// ignored by distance queries, preserving their existing length threshold.
inline void computeFastSegments(const zObjGraph &graph,
                                std::vector<FastSegment> &outSegments) {
  outSegments.clear();
  zFnGraph fn(const_cast<zObjGraph &>(graph));
  zPointArray positions;
  zIntArray edges;
  fn.getVertexPositions(positions);
  fn.getEdgeData(edges);
  if (positions.empty() || edges.empty())
    return;

  outSegments.reserve(edges.size() / 2);
  for (size_t e = 0; e + 1 < edges.size(); e += 2) {
    int v0 = edges[e];
    int v1 = edges[e + 1];
    if (v0 < 0 || v0 >= (int)positions.size() || v1 < 0 ||
        v1 >= (int)positions.size())
      continue;

    float x0 = static_cast<float>(positions[v0].x);
    float y0 = static_cast<float>(positions[v0].y);
    float x1 = static_cast<float>(positions[v1].x);
    float y1 = static_cast<float>(positions[v1].y);

    float dx = x1 - x0;
    float dy = y1 - y0;
    float lenSq = dx * dx + dy * dy;
    FastSegment seg;
    seg.x0 = x0;
    seg.y0 = y0;
    seg.x1 = x1;
    seg.y1 = y1;
    seg.dx = dx;
    seg.dy = dy;
    seg.lenSqInv = (lenSq >= 1e-6f) ? 1.0f / lenSq : 0.0f;
    seg.minX = (std::min)(x0, x1);
    seg.maxX = (std::max)(x0, x1);
    seg.minY = (std::min)(y0, y1);
    seg.maxY = (std::max)(y0, y1);
    outSegments.push_back(seg);
  }
}

// Computes per-point unsigned distance to the nearest FastSegment, minus
// offset.
inline void computeDF(const zPoint *points, int numPoints,
                                     const std::vector<FastSegment> &segments,
                                     zScalarArray &outDistances,
                                     float offset = 0.0f) {
  outDistances.resize(numPoints);
  if (std::none_of(segments.begin(), segments.end(),
                   [](const FastSegment &s) { return s.lenSqInv > 0.0f; })) {
    std::fill(outDistances.begin(), outDistances.end(),
              (std::numeric_limits<zScalar>::max)());
    return;
  }

  const int numSegs = static_cast<int>(segments.size());

  for (int i = 0; i < numPoints; i++) {
    const float px = static_cast<float>(points[i].x);
    const float py = static_cast<float>(points[i].y);
    float minD2 = (std::numeric_limits<float>::max)();

    for (int k = 0; k < numSegs; k++) {
      const FastSegment &seg = segments[k];
      if (seg.lenSqInv == 0.0f) continue;

      float dboxX = (std::max)(0.0f, (std::max)(seg.minX - px, px - seg.maxX));
      float dboxY = (std::max)(0.0f, (std::max)(seg.minY - py, py - seg.maxY));
      if (dboxX * dboxX + dboxY * dboxY >= minD2)
        continue;

      float t =
          ((px - seg.x0) * seg.dx + (py - seg.y0) * seg.dy) * seg.lenSqInv;
      t = (std::max)(0.0f, (std::min)(1.0f, t));
      float projX = seg.x0 + t * seg.dx;
      float projY = seg.y0 + t * seg.dy;
      float dist2 = (px - projX) * (px - projX) + (py - projY) * (py - projY);
      if (dist2 < minD2)
        minD2 = dist2;
    }

    outDistances[i] = static_cast<zScalar>(std::sqrt(minD2) - offset);
  }
}

// Computes per-point inside/outside sign (-1/1) via winding number against
// prepared directed edges. Distance pruning must not be used for winding.
inline void computeSDF(const zPoint *points, int numPoints,
                               const std::vector<FastSegment> &polyEdges,
                               zScalarArray &outPolyField) {
  outPolyField.resize(numPoints);
  if (numPoints == 0)
    return;

  if (polyEdges.empty()) {
    std::fill(outPolyField.begin(), outPolyField.end(), 1.0);
    return;
  }

  const int numEdges = static_cast<int>(polyEdges.size());

  for (int i = 0; i < numPoints; i++) {
    const float px = static_cast<float>(points[i].x);
    const float py = static_cast<float>(points[i].y);

    int windingNumber = 0;

    for (int e = 0; e < numEdges; e++) {
      const FastSegment &edge = polyEdges[e];

      if (py < edge.minY || py > edge.maxY)
        continue;

      if (edge.y0 <= py) {
        if (edge.y1 > py) {
          float cross = edge.dx * (py - edge.y0) -
                        (px - edge.x0) * edge.dy;
          if (cross > 0.0f)
            windingNumber++;
        }
      } else {
        if (edge.y1 <= py) {
          float cross = edge.dx * (py - edge.y0) -
                        (px - edge.x0) * edge.dy;
          if (cross < 0.0f)
            windingNumber--;
        }
      }
    }

    outPolyField[i] = (windingNumber != 0) ? -1.0 : 1.0;
  }
}

// Shared by ratio-based and anchored slots; tangent need not be normalized.
inline bool createSlotSegmentXY(const zPoint &center, const zVector &tangent,
                                zScalar halfLength, zObjGraph &outSlot) {
  const double length = std::hypot(tangent.x, tangent.y);
  if (!std::isfinite(length) || length <= 1e-9 ||
      !std::isfinite(halfLength) || halfLength <= 0.0) return false;
  const double px = -tangent.y / length * halfLength;
  const double py = tangent.x / length * halfLength;
  zPointArray positions = {
      zPoint(center.x + px, center.y + py, center.z),
      zPoint(center.x - px, center.y - py, center.z)};
  zIntArray edges = {0, 1};
  zFnGraph(outSlot).create(positions, edges);
  return true;
}

} // anonymous namespace

//--------------------------
//---- PUBLIC FUNCTIONS
//--------------------------

ZSPACE_INLINE bool createSlot(const zObjGraph &sourceGraph,
                                           zObjGraph &outGraph, zScalar t,
                                           zScalar trimLength) {
  // Extract before clearing the output so sourceGraph may alias outGraph.
  std::vector<GraphPolyline> curves;
  const bool valid = graphPolylines(sourceGraph, curves);
  zPointArray sourcePoints;
  zFnGraph(const_cast<zObjGraph &>(sourceGraph)).getVertexPositions(sourcePoints);
  zFnGraph(outGraph).clear();
  if (!valid || !std::isfinite(t) || !std::isfinite(trimLength) ||
      trimLength <= 0.0) return false;

  const double ratio = std::clamp(static_cast<double>(t), 0.0, 1.0);
  zPointArray slotPoints;
  zIntArray slotEdges;
  for (const GraphPolyline &curve : curves) {
    zPointArray polyline;
    for (int v : curve.Vertices) polyline.push_back(sourcePoints[v]);
    if (curve.Closed) polyline.push_back(polyline.front());
    zScalarArray arc(polyline.size(), 0.0);
    for (size_t i = 1; i < polyline.size(); ++i)
      arc[i] = static_cast<zScalar>(arc[i - 1] + std::sqrt(vecDistSq(polyline[i - 1], polyline[i])));
    if (arc.empty() || !std::isfinite(arc.back()) || arc.back() <= 1e-3) continue;
    zPoint center;
    zVector tangent;
    if (!samplePolyline(polyline.data(), arc.data(), polyline.size(), ratio, center, &tangent)) continue;
    zObjGraph slot;
    if (!createSlotSegmentXY(center, tangent, trimLength, slot)) continue;
    zPointArray points;
    zFnGraph(slot).getVertexPositions(points);
    const int base = static_cast<int>(slotPoints.size());
    slotPoints.insert(slotPoints.end(), points.begin(), points.end());
    slotEdges.push_back(base);
    slotEdges.push_back(base + 1);
  }
  if (slotEdges.empty()) return false;
  zFnGraph(outGraph).create(slotPoints, slotEdges);
  return true;
}

ZSPACE_INLINE bool validateFieldSettings(const SlicerSettings &s) {
  const auto valid = [](double v) { return std::isfinite(v) && v >= MinPhysicalParameter; };
  if (!valid(s.FieldCellSize) || s.BoundaryWidth.empty()) {
    zTsUtils::printDebug("[validateFieldSettings] FieldCellSize must be finite and at least 1e-4 mm; BoundaryWidth must not be empty.");
    return false;
  }
  const auto validWidth = [&](double w, bool boundary) {
    return valid(w) && valid(boundary && !s.DoubleBoundary ? w : 0.5 * w + 0.5 * s.OverlapWidth);
  };
  for (float w : s.BoundaryWidth) if (!validWidth(w, true)) {
    zTsUtils::printDebug("[validateFieldSettings] Boundary width or effective bead width is below 1e-4 mm or non-finite.");
    return false;
  }
  if (s.BracingCurveWidths.empty()) {
    if (validWidth(s.BracingWidth, false)) return true;
  } else {
    bool validAll = true;
    for (float w : s.BracingCurveWidths) validAll = validAll && validWidth(w, false);
    if (validAll) return true;
  }
  zTsUtils::printDebug("[validateFieldSettings] Bracing width or effective bead width is below 1e-4 mm or non-finite.");
  return false;
}

ZSPACE_INLINE float widthOfGroup(const std::vector<float> &widths, int id) {
  if (widths.empty()) {
    zTsUtils::printDebug("[widthOfGroup] Empty width table.");
    return std::numeric_limits<float>::quiet_NaN();
  }
  const int index = (std::max)(0, (std::min)(id, int(widths.size()) - 1));
  if (!std::isfinite(widths[index]) || widths[index] < MinPhysicalParameter) {
    zTsUtils::printDebug("[widthOfGroup] Width must be finite and at least 1e-4 mm.");
    return std::numeric_limits<float>::quiet_NaN();
  }
  return widths[index];
}

ZSPACE_INLINE zScalar printedWidth(zScalar width, bool isBoundary,
                                   bool doubleBoundary, zScalar overlapWidth) {
  const zScalar result = isBoundary && !doubleBoundary ? width : zScalar(0.5 * width + 0.5 * overlapWidth);
  if (!std::isfinite(result) || result < MinPhysicalParameter) {
    zTsUtils::printDebug("[printedWidth] Effective width must be finite and at least 1e-4 mm.");
    return std::numeric_limits<zScalar>::quiet_NaN();
  }
  return result;
}

ZSPACE_INLINE int loopWidthSegments(const zObjGraph &loop,
                                    const zIntArray &edgeWidthIds,
                                    const std::vector<float> &widths,
                                    std::vector<FieldSegment> &outSegments) {
  zFnGraph fn(const_cast<zObjGraph &>(loop));
  zPointArray positions;
  zIntArray edges;
  fn.getVertexPositions(positions);
  fn.getEdgeData(edges);
  const size_t firstSegment = outSegments.size();
  // Group actual directed edges by width. Their distance field is the union
  // of the pieces; no connections are added between disjoint pieces.
  std::unordered_map<int, zIntArray> groups;
  for (int e = 0; e < static_cast<int>(edges.size() / 2); ++e) {
    const int a = edges[2 * e], b = edges[2 * e + 1];
    if (a < 0 || b < 0 || a >= static_cast<int>(positions.size()) ||
        b >= static_cast<int>(positions.size()) || a == b) continue;
    const int id = e < static_cast<int>(edgeWidthIds.size()) ? edgeWidthIds[e] : 0;
    groups[id].push_back(a);
    groups[id].push_back(b);
  }
  zIntArray groupIds;
  for (const auto &group : groups) groupIds.push_back(group.first);
  std::sort(groupIds.begin(), groupIds.end());
  for (int id : groupIds) {
    zPointArray points;
    zIntArray connects;
    std::unordered_map<int, int> local;
    for (int v : groups[id]) {
      auto found = local.find(v);
      if (found == local.end()) {
        const int index = static_cast<int>(points.size());
        local.emplace(v, index);
        points.push_back(positions[v]);
        connects.push_back(index);
      } else connects.push_back(found->second);
    }
    zObjGraph segment;
    zFnGraph(segment).create(points, connects);
    const float w = widthOfGroup(widths, id);
    outSegments.emplace_back(segment, 0, w, 0.5 * w);
  }
  return static_cast<int>(outSegments.size() - firstSegment);
}

ZSPACE_INLINE void computeMeshDomainSigns(const zObjMeshScalarField &field,
    const zObjMesh &mesh, zScalarArray &outSigns) {
  zFnMesh fnMesh(const_cast<zObjMesh &>(mesh));
  zPointArray positions;
  zIntArray connects, counts;
  fnMesh.getVertexPositions(positions);
  fnMesh.getPolygonData(connects, counts);
  std::unordered_map<unsigned long long, int> edgeUse;
  size_t offset = 0;
  for (int count : counts) {
    for (int i = 0; i < count; ++i)
      ++edgeUse[makeEdgeKey(connects[offset + i], connects[offset + (i + 1) % count])];
    offset += count;
  }
  zIntArray edges;
  offset = 0;
  for (int count : counts) {
    for (int i = 0; i < count; ++i) {
      const int a = connects[offset + i], b = connects[offset + (i + 1) % count];
      if (edgeUse[makeEdgeKey(a, b)] != 1) continue;
      edges.push_back(a);
      edges.push_back(b);
    }
    offset += count;
  }
  zFnMesh fnGrid(const_cast<zObjMeshScalarField &>(field));
  outSigns.assign(fnGrid.numVertices(), 1.0);
  if (edges.empty() || positions.empty() || outSigns.empty()) return;
  zObjGraph domainBoundary;
  zFnGraph(domainBoundary).create(positions, edges);
  std::vector<FastSegment> segments;
  computeFastSegments(domainBoundary, segments);
  zPointArray gridPositions;
  fnGrid.getVertexPositions(gridPositions);
  computeSDF(gridPositions.data(), static_cast<int>(gridPositions.size()), segments, outSigns);
}

ZSPACE_INLINE bool locateGraphLoopAnchor(const zObjGraph &graph, const zPoint &point,
    zPointArray &outLoop, zIntArray &outEdgeIds, LoopAnchor &outAnchor) {
  outLoop.clear();
  outEdgeIds.clear();
  outAnchor = LoopAnchor();
  std::vector<GraphPolyline> curves;
  if (!graphPolylines(graph, curves)) return false;
  zPointArray positions;
  zFnGraph(const_cast<zObjGraph &>(graph)).getVertexPositions(positions);
  double bestDistance = (std::numeric_limits<double>::max)();
  for (const GraphPolyline &curve : curves) {
    if (!curve.Closed || curve.Vertices.size() < 3) continue;
    zPointArray points;
    for (int v : curve.Vertices) points.push_back(positions[v]);
    LoopAnchor candidate;
    if (!locateLoopAnchor(points, point, candidate) || candidate.distance >= bestDistance) continue;
    bestDistance = candidate.distance;
    outLoop.swap(points);
    outEdgeIds = curve.Edges;
    outAnchor = candidate;
  }
  return outAnchor.valid();
}

ZSPACE_INLINE bool locateLoopAnchor(const zPointArray &loop,
                                    const zPoint &anchor, LoopAnchor &out,
                                    const zPointArray *senseLoop3D) {
  out = LoopAnchor();
  const zPointArray &loop3D = loop;
  const int n = static_cast<int>(loop3D.size());
  if (n < 3)
    return false;

  double bestDistSq = (std::numeric_limits<double>::max)();
  for (int i = 0; i < n; i++) {
    const zPoint &a = loop3D[i];
    const zPoint &b = loop3D[(i + 1) % n];
    const double dx = b.x - a.x, dy = b.y - a.y, dz = b.z - a.z;
    const double lenSq = dx * dx + dy * dy + dz * dz;
    double t = 0.0;
    if (lenSq > 1e-18)
      t = ((anchor.x - a.x) * dx + (anchor.y - a.y) * dy +
           (anchor.z - a.z) * dz) /
          lenSq;
    t = (std::max)(0.0, (std::min)(1.0, t));
    const double px = a.x + t * dx - anchor.x;
    const double py = a.y + t * dy - anchor.y;
    const double pz = a.z + t * dz - anchor.z;
    const double distSq = px * px + py * py + pz * pz;
    if (distSq < bestDistSq) {
      bestDistSq = distSq;
      out.edge = i;
      out.t = t;
      out.point = zPoint(a.x + t * dx, a.y + t * dy, a.z + t * dz);
    }
  }
  // An unrolled/chart image may be mirrored: take the sense from the world
  // loop when given (same vertex order).
  const bool useSense = senseLoop3D && int(senseLoop3D->size()) == n;
  out.ccw = signedArea2D(useSense ? *senseLoop3D : loop3D) >= 0.0;
  out.distance = std::sqrt(bestDistSq);
  return out.valid();
}

ZSPACE_INLINE bool createAnchoredSlot(const zPointArray &loop2D,
                                      const LoopAnchor &anchor, double offset,
                                      zScalar halfLength, zObjGraph &outSlot,
                                      zPoint *outCenter) {
  zFnGraph(outSlot).clear();
  const int n = static_cast<int>(loop2D.size());
  if (!anchor.valid() || n < 3 || anchor.edge >= n ||
      !std::isfinite(anchor.t) || !std::isfinite(offset) ||
      !std::isfinite(halfLength) || halfLength <= 0.0)
    return false;

  // Cumulative XY arc length around the closed loop (the fields live in XY).
  zPointArray points = loop2D;
  points.push_back(loop2D.front());
  zScalarArray arc(n + 1, 0.0);
  for (int i = 0; i < n; i++) {
    const zPoint &a = loop2D[i];
    const zPoint &b = loop2D[(i + 1) % n];
    arc[i + 1] = arc[i] + std::sqrt((b.x - a.x) * (b.x - a.x) +
                                    (b.y - a.y) * (b.y - a.y));
  }
  const double total = arc[n];
  if (!std::isfinite(total) || total <= 1e-9)
    return false;

  // Increasing index is counter-clockwise exactly when the loop is.
  const double start =
      arc[anchor.edge] + anchor.t * (arc[anchor.edge + 1] - arc[anchor.edge]);
  double s = std::fmod(start + (anchor.ccw ? offset : -offset), total);
  if (s < 0.0)
    s += total;

  zPoint center;
  zVector tangent;
  if (!std::isfinite(s) ||
      !samplePolyline(points.data(), arc.data(), points.size(), s / total,
                      center, &tangent, 1e-9) ||
      !createSlotSegmentXY(center, tangent, halfLength, outSlot)) return false;
  if (outCenter)
    *outCenter = center;
  return true;
}

ZSPACE_INLINE bool locateSlotAnchor(const zPointArray &floorLoop3D,
                                    const zPoint &anchor, bool hasAnchor,
                                    const std::string &source,
                                    LoopAnchor &out) {
  const std::string tag = "[locateSlotAnchor] ";
  if (!hasAnchor) {
    out = LoopAnchor();
    printDebug(tag + "slot point is required; no default position is used.");
    return false;
  }
  if (!locateLoopAnchor(floorLoop3D, anchor, out)) {
    printDebug(tag + "could not locate the slot anchor on the floor "
                     "boundary - no boundary slot.");
    return false;
  }
  std::ostringstream msg;
  msg << "boundary slot anchored on floor edge " << out.edge
      << " (t=" << out.t << ", anchor " << out.distance << " away)";

  return true;
}

ZSPACE_INLINE bool createParitySlot(const zPointArray &loop2D,
                                    const LoopAnchor &anchor, bool isEven,
                                    double offsetLength, zScalar w,
                                    zObjGraph &outSlot) {
  return createAnchoredSlot(loop2D, anchor,
                            isEven ? offsetLength : -offsetLength, w,
                            outSlot);
}

ZSPACE_INLINE bool computeGridBounds(const zPoint &minBB, const zPoint &maxBB,
                                     double scaleFactor, double cellSize,
                                     zDomain<zPoint> &outBounds, int &outResX,
                                     int &outResY) {
  if (minBB.x >= maxBB.x || minBB.y >= maxBB.y)
    return false;

  const double centerX = (minBB.x + maxBB.x) * 0.5;
  const double centerY = (minBB.y + maxBB.y) * 0.5;
  const double rawWidth = (maxBB.x - minBB.x) * scaleFactor;
  const double rawHeight = (maxBB.y - minBB.y) * scaleFactor;

  if (!std::isfinite(cellSize) || cellSize < MinPhysicalParameter) {
    zTsUtils::printDebug("[computeGridBounds] Cell size must be finite and at least 1e-4 mm.");
    return false;
  }
  const double cell = cellSize;
  const int resX = (std::max)(1, static_cast<int>(std::ceil(rawWidth / cell)));
  const int resY =
      (std::max)(1, static_cast<int>(std::ceil(rawHeight / cell)));

  const double halfW = (resX * cell) * 0.5;
  const double halfH = (resY * cell) * 0.5;

  outBounds = zDomain<zPoint>(zPoint(centerX - halfW, centerY - halfH, 0.0),
                              zPoint(centerX + halfW, centerY + halfH, 0.0));
  outResX = resX;
  outResY = resY;
  return true;
}

ZSPACE_INLINE void
computeSegmentDistances(const zObjMeshScalarField &field,
                        const std::vector<FieldSegment> &segments,
                        std::vector<zScalarArray> &outSegDistances,
                        zScalarArray *outBoundarySigns) {
  if (outBoundarySigns) outBoundarySigns->clear();
  zFnMesh fnMesh(const_cast<zObjMeshScalarField &>(field));
  zPointArray gridPositions;
  fnMesh.getVertexPositions(gridPositions);
  const zPoint *points = gridPositions.data();
  const int n = fnMesh.numVertices();
  const int numSegs = static_cast<int>(segments.size());
  outSegDistances.clear();
  outSegDistances.resize(numSegs);
  if (!points || n <= 0)
    return;

  std::vector<FastSegment> fastEdges;
  std::vector<FastSegment> boundaryEdges;
  for (int k = 0; k < numSegs; k++) {
    computeFastSegments(segments[k].graph, fastEdges);
    if (outBoundarySigns && segments[k].type == 0)
      boundaryEdges.insert(boundaryEdges.end(), fastEdges.begin(), fastEdges.end());
    computeDF(points, n, fastEdges, outSegDistances[k]);
  }
  if (outBoundarySigns)
    computeSDF(points, n, boundaryEdges, *outBoundarySigns);
}

ZSPACE_INLINE void computeWidthField(
    const zObjMeshScalarField &field,
    const std::vector<zScalarArray> &segDistances,
    const std::vector<FieldSegment> &segments, bool doubleBoundary,
    zScalar overlapWidth, zScalarArray *outWidthField,
    zScalarArray *outIdField) {
  if (!outWidthField && !outIdField) return;
  const int n = zFnMeshScalarField(const_cast<zObjMeshScalarField &>(field)).numFieldValues();
  if (outWidthField) outWidthField->assign(n, 0.0);
  if (outIdField) outIdField->assign(n, 0.0);
  if (segments.empty() || segDistances.size() != segments.size()) return;

  zScalarArray minMetrics(n, (std::numeric_limits<zScalar>::max)());
  for (size_t k = 0; k < segments.size(); ++k) {
    if (segDistances[k].size() < static_cast<size_t>(n)) continue;
    const FieldSegment &segment = segments[k];
    const zScalar width = segment.width;
    if (!std::isfinite(width) || width < MinPhysicalParameter) {
      zTsUtils::printDebug("[computeWidthField] Segment width must be finite and at least 1e-4 mm.");
      if (outWidthField) outWidthField->clear();
      if (outIdField) outIdField->clear();
      return;
    }
    const zScalar beadWidth = printedWidth(segment.width, segment.type == 0,
                                           doubleBoundary, overlapWidth);
    for (int i = 0; i < n; ++i) {
      // Keep division and strict comparison to preserve the winning-group ties.
      const zScalar metric = (segDistances[k][i] - segment.offset) / width;
      if (metric < minMetrics[i]) {
        minMetrics[i] = metric;
        if (outWidthField) (*outWidthField)[i] = beadWidth;
        if (outIdField) (*outIdField)[i] = static_cast<zScalar>(k);
      }
    }
  }
}

ZSPACE_INLINE void computePathField(
    const zObjMeshScalarField &field, const zObjGraph &boundaryPolygon,
    const std::vector<FieldSegment> &segments,
    const std::vector<zScalarArray> &segDistances, bool doubleBoundary,
    zScalar overlapWidth, zObjGraph *outTrimSlotsGraph,
    zScalarArray *outPathField, zScalarArray *outTrackIdField,
    const zObjGraph *boundarySlot, zScalar boundarySlotWidth,
    const zScalarArray *boundarySigns) {
  if (outTrimSlotsGraph) zFnGraph(*outTrimSlotsGraph).clear();
  zObjMeshScalarField &grid = const_cast<zObjMeshScalarField &>(field);
  zFnMeshScalarField fnField(grid);
  zPointArray gridPositions;
  zFnMesh(grid).getVertexPositions(gridPositions);
  const zPoint *points = gridPositions.data();
  const int n = fnField.numFieldValues();
  if (outPathField) outPathField->assign(n, 1.0); // No zero-isocontour on failure.
  if (outTrackIdField) outTrackIdField->assign(n, 0.0);
  if (!points || n <= 0 || segments.empty() ||
      segments.size() != segDistances.size()) return;
  for (const zScalarArray &distances : segDistances)
    if (distances.size() < static_cast<size_t>(n)) return;
  if (!outPathField && !outTrackIdField && !outTrimSlotsGraph) return;

  struct BoundaryOffsets { size_t index; zScalar outer, inner; };
  std::vector<BoundaryOffsets> boundaries;
  std::vector<size_t> bracing;
  for (size_t k = 0; k < segments.size(); ++k) {
    const FieldSegment &segment = segments[k];
    if (segment.type == 0) {
      const zScalar bead = printedWidth(segment.width, true, doubleBoundary, overlapWidth);
      const zScalar halfBead = static_cast<zScalar>(0.5) * bead;
      boundaries.push_back({k, halfBead, segment.width - halfBead});
    } else {
      bracing.push_back(k);
    }
  }

  zScalarArray computedSigns;
  const zScalarArray *signs = boundarySigns;
  if (!boundaries.empty() && (!signs || signs->size() != static_cast<size_t>(n))) {
    std::vector<FastSegment> edges;
    computeFastSegments(boundaryPolygon, edges);
    computeSDF(points, n, edges, computedSigns);
    signs = &computedSigns;
  }

  // Only bracing requires these full-grid buffers. Reuse them across ribs.
  zScalarArray bracingField, ribField, trimmedRib, slotDistances;
  if (!bracing.empty()) bracingField.assign(n, (std::numeric_limits<zScalar>::max)());
  std::vector<FastSegment> slotEdges;
  std::vector<zObjGraph> slotGraphs;
  for (size_t k : bracing) {
    const FieldSegment &segment = segments[k];
    const zScalar bead = printedWidth(segment.width, false, doubleBoundary,
                                      doubleBoundary ? overlapWidth : zScalar(0.0));
    const zScalar offset = 0.5 * segment.width - 0.5 * bead;
    ribField.resize(n);
    for (int i = 0; i < n; ++i) ribField[i] = segDistances[k][i] - offset;
    const zScalarArray *rib = &ribField;
    if (segment.trim >= 0.0 && segment.trim <= 1.0) {
      zObjGraph slot;
      if (createSlot(segment.graph, slot, segment.trim, segment.width * 1.5)) {
        if (outTrimSlotsGraph) slotGraphs.push_back(slot);
        computeFastSegments(slot, slotEdges);
        computeDF(points, n, slotEdges, slotDistances, static_cast<float>(0.5 * bead));
        fnField.boolean_subtract(ribField, slotDistances, trimmedRib, false);
        rib = &trimmedRib;
      }
    }
    for (int i = 0; i < n; ++i)
      if ((*rib)[i] < bracingField[i]) bracingField[i] = (*rib)[i];
  }

  zScalarArray discardedPath;
  zScalarArray &path = outPathField ? *outPathField : discardedPath;
  if (boundaries.empty()) {
    if (!bracing.empty()) path.swap(bracingField);
  } else {
    zScalarArray outer(n), inner, scratch;
    if (doubleBoundary) inner.resize(n);
    for (int i = 0; i < n; ++i) {
      const bool inside = (*signs)[i] <= 0.0;
      double maxOuter = -(std::numeric_limits<double>::max)();
      double maxInner = -(std::numeric_limits<double>::max)();
      for (const BoundaryOffsets &boundary : boundaries) {
        const zScalar d = segDistances[boundary.index][i];
        const zScalar signedDistance = inside ? -d : d;
        const zScalar outerValue = signedDistance + boundary.outer;
        if (outerValue > maxOuter) maxOuter = outerValue;
        if (doubleBoundary) {
          const zScalar innerValue = signedDistance + boundary.inner;
          if (innerValue > maxInner) maxInner = innerValue;
        }
      }
      outer[i] = maxOuter;
      if (doubleBoundary) inner[i] = maxInner;
    }

    if (doubleBoundary) {
      if (!bracing.empty()) {
        fnField.boolean_subtract(inner, bracingField, scratch, false);
        inner.swap(scratch);
      }
      if (outTrackIdField)
        for (int i = 0; i < n; ++i)
          (*outTrackIdField)[i] = outer[i] >= -inner[i] ? 1.0 : 0.0;
      fnField.boolean_subtract(outer, inner, scratch, false);
      outer.swap(scratch); // Reuse outer's storage for the combined band.
      if (boundarySlot && boundarySlotWidth > 0.0 &&
          zFnGraph(const_cast<zObjGraph &>(*boundarySlot)).numEdges() > 0) {
        const zScalar bead = printedWidth(boundarySlotWidth, true, true, overlapWidth);
        if (outTrimSlotsGraph) slotGraphs.push_back(*boundarySlot);
        computeFastSegments(*boundarySlot, slotEdges);
        computeDF(points, n, slotEdges, slotDistances, static_cast<float>(0.5 * bead));
        if (outTrackIdField)
          for (int i = 0; i < n; ++i)
            if (outer[i] < -slotDistances[i]) (*outTrackIdField)[i] = 1.0;
        fnField.boolean_subtract(outer, slotDistances, path, false);
      } else {
        path.swap(outer);
      }
    } else if (!bracing.empty()) {
      fnField.boolean_subtract(outer, bracingField, path, false);
    } else {
      path.swap(outer);
    }
  }
  if (outTrimSlotsGraph && !slotGraphs.empty())
    combineGraphObjects(slotGraphs, *outTrimSlotsGraph);
}

ZSPACE_INLINE void
extractRawIsoContour(const zObjMeshScalarField &field,
                     const zScalarArray &pathVals,
                     const zScalarArray &widthVals,
                     const zScalarArray &trackIdVals, zPointArray &outPositions,
                     zIntArray &outEdgeConnects, zFloatArray &outVertexWidths,
                     zFloatArray &outVertexTrackIds) {
  outPositions.clear();
  outEdgeConnects.clear();
  outVertexWidths.clear();
  outVertexTrackIds.clear();

  zFnMesh fnMesh(const_cast<zObjMeshScalarField &>(field));
  // Workers read owned snapshots; acquiring a mutable v3 raw view would
  // invalidate the shared grid topology even though extraction only reads it.
  zPointArray gridPositions;
  fnMesh.getVertexPositions(gridPositions);
  const zPoint *rawPositions = gridPositions.data();
  const int numVerts = fnMesh.numVertices();
  const int numFaces = fnMesh.numPolygons();
  const bool hasTrackIds = trackIdVals.size() == static_cast<size_t>(numVerts);
  if (!rawPositions || numVerts == 0 || numFaces == 0 ||
      pathVals.size() != static_cast<size_t>(numVerts) ||
      widthVals.size() != static_cast<size_t>(numVerts)) {
    return;
  }

  std::unordered_map<unsigned long long, int> edgeCrossMap;
  edgeCrossMap.reserve(numFaces / 4);

  auto getEdgeIntersection = [&](int v0, int v1) -> int {
    const double s0 = pathVals[v0];
    const double s1 = pathVals[v1];
    const double denom = s1 - s0;
    const double t = (std::abs(denom) > 1e-9)
                         ? (std::max)(0.0, (std::min)(1.0, -s0 / denom))
                         : 0.5;
    const unsigned long long edgeKey = t == 0.0 ? makeEdgeKey(v0, v0)
        : t == 1.0 ? makeEdgeKey(v1, v1) : makeEdgeKey(v0, v1);
    auto it = edgeCrossMap.find(edgeKey);
    if (it != edgeCrossMap.end()) return it->second;

    const zPoint &p0 = rawPositions[v0], &p1 = rawPositions[v1];
    outPositions.push_back(zPoint(p0.x + t * (p1.x - p0.x),
                                  p0.y + t * (p1.y - p0.y),
                                  p0.z + t * (p1.z - p0.z)));
    outVertexWidths.push_back(static_cast<float>(
        widthVals[v0] + t * (widthVals[v1] - widthVals[v0])));

    if (hasTrackIds) {
      const double id0 = trackIdVals[v0];
      const double id1 = trackIdVals[v1];
      const double idInterp = id0 + t * (id1 - id0);
      outVertexTrackIds.push_back(idInterp >= 0.5 ? 1.0f : 0.0f);
    } else {
      outVertexTrackIds.push_back(0.0f);
    }

    const int newIdx = static_cast<int>(outPositions.size()) - 1;
    edgeCrossMap[edgeKey] = newIdx;
    return newIdx;
  };

  zIntArray polyConnects, polyCounts;
  fnMesh.getPolygonData(polyConnects, polyCounts);

  size_t cOffset = 0;
  for (int count : polyCounts) {
    for (int tri = 1; tri + 1 < count; tri++) {
      const int v0 = polyConnects[cOffset];
      const int v1 = polyConnects[cOffset + tri];
      const int v2 = polyConnects[cOffset + tri + 1];

      const double s0 = pathVals[v0], s1 = pathVals[v1], s2 = pathVals[v2];
      const int mask = (s0 > 0.0) | ((s1 > 0.0) << 1) | ((s2 > 0.0) << 2);
      if (mask == 0 || mask == 7)
        continue;

      int p0 = -1, p1 = -1;
      switch (mask) {
      case 1:
      case 6:
        p0 = getEdgeIntersection(v0, v1);
        p1 = getEdgeIntersection(v2, v0);
        break;
      case 2:
      case 5:
        p0 = getEdgeIntersection(v0, v1);
        p1 = getEdgeIntersection(v1, v2);
        break;
      case 3:
      case 4:
        p0 = getEdgeIntersection(v1, v2);
        p1 = getEdgeIntersection(v2, v0);
        break;
      }

      if (p0 != p1 && p0 >= 0 && p1 >= 0) {
        outEdgeConnects.push_back(p0);
        outEdgeConnects.push_back(p1);
      }
    }
    cOffset += count;
  }
}

ZSPACE_INLINE bool extractPrintPaths(const zObjMeshScalarField &fieldGrid,
    const zScalarArray &pathField, const zScalarArray &widthField,
    const zScalarArray &trackField, std::vector<PrintPath> &outPaths, int layer, std::vector<std::string> *debug) {
  const auto report = [&](const std::string &reason) {
    if (debug) debug->push_back("[extractPrintPaths] " + reason);
    else printDebug("[extractPrintPaths] " + (layer >= 0 ? "layer=" + std::to_string(layer) + " " : std::string()) + reason);
  };
  const auto fail = [&](bool invalidInput = false) {
    report("openCurves=0 closedCurves=0 invalidGraphs=" + std::to_string(invalidInput ? 1 : 0));
    return false;
  };
  outPaths.clear();
  const int n = zFnMesh(const_cast<zObjMeshScalarField &>(fieldGrid)).numVertices();
  if (n == 0 || pathField.size() != size_t(n) || widthField.size() != size_t(n) ||
      trackField.size() != size_t(n))
    return fail(true);
  zFnMesh gridMesh(const_cast<zObjMeshScalarField &>(fieldGrid));
  zPointArray gridPoints;
  zIntArray gridConnects, gridCounts;
  gridMesh.getVertexPositions(gridPoints);
  gridMesh.getPolygonData(gridConnects, gridCounts);
  double cellLength = (std::numeric_limits<double>::max)();
  size_t gridOffset = 0;
  for (int count : gridCounts) {
    for (int i = 0; i < count; ++i) {
      const double length = std::sqrt(vecDistSq(gridPoints[gridConnects[gridOffset + i]],
          gridPoints[gridConnects[gridOffset + (i + 1) % count]]));
      if (length > 0.0) cellLength = (std::min)(cellLength, length);
    }
    gridOffset += count;
  }
  if (cellLength == (std::numeric_limits<double>::max)()) {
    return fail(true);
  }
  const double weldTolerance = (std::max)(1e-3, 1e-4 * cellLength);
  const double spurLength = 2.0 * cellLength;

  zPointArray positions;
  zIntArray edges;
  zFloatArray widths, tracks;
  extractRawIsoContour(fieldGrid, pathField, widthField, trackField,
                      positions, edges, widths, tracks);
  if (edges.empty()) return fail();
  extractTrackPaths(positions, edges, widths, tracks, weldTolerance, spurLength, outPaths);
  if (outPaths.empty()) return fail();
  size_t openCount = 0, closedCount = 0, invalidGraphCount = 0;
  size_t longestIndex = outPaths.size();
  double longestLength = 0.0;
  for (size_t index = 0; index < outPaths.size(); ++index) {
    const PrintPath &path = outPaths[index];
    std::vector<GraphPolyline> curves;
    if (!graphPolylines(path.graph, curves) || curves.size() != 1) {
      ++invalidGraphCount;
      continue;
    }
    const GraphPolyline &curve = curves.front();
    if (!curve.Closed) {
      ++openCount;
      continue;
    }
    zPointArray points;
    zFnGraph(const_cast<zObjGraph &>(path.graph)).getVertexPositions(points);
    auto ids = curve.Vertices;
    if (ids.size() > 1 && ids.front() == ids.back()) ids.pop_back();
    bool valid = ids.size() >= 3 && ids.size() == points.size() &&
        path.widths.size() == points.size() &&
        (path.trackIds.empty() || path.trackIds.size() == points.size());
    double length = 0.0;
    for (size_t i = 0; valid && i < ids.size(); ++i) {
      const int a = ids[i], b = ids[(i + 1) % ids.size()];
      if (a < 0 || b < 0 || size_t(a) >= points.size() || size_t(b) >= points.size()) {
        valid = false;
        break;
      }
      length += std::sqrt(vecDistSq(points[a], points[b]));
    }
    if (!valid || !std::isfinite(length) || length <= 0.0) {
      ++invalidGraphCount;
      continue;
    }
    ++closedCount;
    if (length > longestLength) {
      longestLength = length;
      longestIndex = index;
    }
  }
  if (outPaths.size() != 1 || longestIndex == outPaths.size())
    report("openCurves=" + std::to_string(openCount) +
        " closedCurves=" + std::to_string(closedCount) +
        " invalidGraphs=" + std::to_string(invalidGraphCount) +
        " kept=" + std::to_string(longestIndex < outPaths.size() ? 1 : 0) +
        " discarded=" + std::to_string(outPaths.size() - (longestIndex < outPaths.size() ? 1 : 0)));
  if (longestIndex == outPaths.size()) {
    outPaths.clear();
    return false;
  }
  if (longestIndex != 0) std::swap(outPaths.front(), outPaths[longestIndex]);
  outPaths.resize(1);
  return true;
}

ZSPACE_INLINE bool mapPrintPaths(const std::vector<PrintPath> &paths2D,
    zObjMesh &chart2D, zObjMesh &mesh3D, const zFloatArray *vertexHeights,
    std::vector<PrintPath> &outPaths3D, int layer, const zVectorArray *vertexNormals) {
  outPaths3D = paths2D;
  const zFloatArray emptyHeights;
  if (!vertexHeights) vertexHeights = &emptyHeights;
  bool success = true;
  for (size_t i = 0; i < outPaths3D.size(); ++i) {
    auto &path = outPaths3D[i];
    std::string reason;
    if (!remapGraph(path.graph, chart2D, mesh3D, &path.normals,
                    vertexHeights, &path.heights, &reason, vertexNormals)) {
      printDebug("[mapPrintPaths] " + (layer >= 0 ? "layer=" + std::to_string(layer) + " " : std::string()) + reason);
      path = PrintPath();
      success = false;
    }
  }
  return success;
}

} // namespace zTsField
} // namespace zSpace
