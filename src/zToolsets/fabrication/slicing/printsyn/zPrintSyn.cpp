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
#include <zspace/zToolsets/fabrication/slicing/printsyn/zPrintSyn.h>

#ifndef ZSPACE_INLINE
#define ZSPACE_INLINE ZSPACE_TOOLSETS_INLINE
#endif

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <limits>
#include <map>
#include <memory>
#include "../utils/SurfaceQuery.h"
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#if defined(_OPENMP)
#include <omp.h>
#endif

namespace zSpace {

using namespace zTsUtils;

//--------------------------
//---- PRIVATE HELPERS (grouped by the step below that uses them; a free
//---- function must be fully defined before its first use here, since none
//---- of these are declared in zPrintSyn.h)
//--------------------------

namespace {

// ---- resamplePrintPath helpers ----

// Flags original vertex IDs whose turn angle is below the feature threshold.
void detectFeatureAngles(const zPointArray &positions,
                         const std::vector<int> &sequence,
                         float angleThresholdDeg, bool closed,
                         std::vector<int> &outFeatureFlags) {
  const int N = static_cast<int>(sequence.size());
  outFeatureFlags.assign(positions.size(), 0);
  if (N < 3)
    return;

  const double thresholdRad =
      static_cast<double>(angleThresholdDeg) * (3.14159265358979323846 / 180.0);

  for (int i = 0; i < N; i++) {
    int prevId = sequence[(i - 1 + N) % N];
    int curId = sequence[i];
    int nextId = sequence[(i + 1) % N];
    if (!closed && (i == 0 || i == N - 1)) {
      outFeatureFlags[curId] = 1;
      continue;
    }

    // Features on short edges are still features. turnAngleAt's 1e-3
    // physical cutoff is unsuitable here; reject only coincident edges.
    const Eigen::Vector3d a = toEigenVec(positions[prevId]) - toEigenVec(positions[curId]);
    const Eigen::Vector3d b = toEigenVec(positions[nextId]) - toEigenVec(positions[curId]);
    const double product = a.norm() * b.norm();
    if (product > 0.0 && std::isfinite(product) &&
        std::acos((std::max)(-1.0, (std::min)(1.0, a.dot(b) / product))) < thresholdRad) {
      outFeatureFlags[curId] = 1;
    }
  }
}

// Fixed-spacing arc-length resample of a single ioPath.graph, carrying
// widths/heights/normals along.
bool resampleSinglePath(PrintPath &ioPath, float resampleLength,
                        float featureAngleTol) {
  zFnGraph fnGraph(ioPath.graph);
  const int N = fnGraph.numVertices();
  if (!std::isfinite(resampleLength) || resampleLength < zTsField::MinPhysicalParameter) {
    zTsUtils::printDebug("[resampleSinglePath] ResampleLength must be finite and at least 1e-4 mm.");
    return false;
  }
  const float sampleLength = resampleLength;
  if (N < 2)
    return false;
  if ((!ioPath.widths.empty() && ioPath.widths.size() != size_t(N)) ||
      (!ioPath.heights.empty() && ioPath.heights.size() != size_t(N)) ||
      (!ioPath.normals.empty() && ioPath.normals.size() != size_t(N)) ||
      (!ioPath.trackIds.empty() && ioPath.trackIds.size() != size_t(N)) ||
      (!ioPath.featureFlags.empty() && ioPath.featureFlags.size() != size_t(N))) return false;

  std::vector<GraphPolyline> curves;
  if (!graphPolylines(ioPath.graph, curves) || curves.size() != 1 ||
      curves[0].Vertices.size() != static_cast<size_t>(N)) return false;
  const bool closed = curves[0].Closed;

  zPointArray positions;
  fnGraph.getVertexPositions(positions);
  for (const auto &p : positions)
    if (!std::isfinite(p.x) || !std::isfinite(p.y) || !std::isfinite(p.z)) return false;

  const zIntArray &sequence = curves[0].Vertices;

  std::vector<int> featureFlags;
  detectFeatureAngles(positions, sequence, featureAngleTol, closed, featureFlags);
  for (int i = 0; i < N; ++i) {
    // Open endpoints are mandatory, but are not angle features.
    if (!closed && (i == sequence.front() || i == sequence.back())) featureFlags[i] = 0;
    if (!ioPath.featureFlags.empty()) featureFlags[i] |= ioPath.featureFlags[i];
  }

  const zFloatArray &weights = ioPath.widths;
  const zFloatArray &heights = ioPath.heights;
  const zVectorArray &normals = ioPath.normals;
  const zFloatArray &trackIds = ioPath.trackIds;
  const bool hasHeights = !heights.empty();
  const bool hasNormals = !normals.empty();
  const bool hasTrackIds = !trackIds.empty();

  zPointArray outPoints;
  zFloatArray outWeights;
  zFloatArray outHeights;
  zVectorArray outNormals;
  zFloatArray outTrackIds;
  zIntArray outEdges;
  zIntArray outFeatures;
  std::vector<char> outRequired;

  const auto popPoint = [&]() {
    outPoints.pop_back(); outWeights.pop_back();
    if (hasHeights) outHeights.pop_back();
    if (hasNormals) outNormals.pop_back();
    if (hasTrackIds) outTrackIds.pop_back();
    outFeatures.pop_back(); outRequired.pop_back();
  };

  auto appendPoint = [&](const zPoint &p, float w, float h, const zVector &nrm,
                          float trk, int features, bool required) {
    // A nearby ordinary sample yields to an exact feature/endpoint. Distinct
    // required points never weld, however short the edge between them is.
    const double clearance = (std::min)(0.01, double(sampleLength) * 0.1);
    while (required && !outPoints.empty() && !outRequired.back() &&
           vecDistSq(p, outPoints.back()) < clearance * clearance) popPoint();
    if (!outPoints.empty()) {
      if (vecDistSq(p, outPoints.back()) == 0.0) {
        outFeatures.back() |= features;
        outRequired.back() |= required;
        return;
      }
    }
    outPoints.push_back(p);
    outWeights.push_back(w);
    if (hasHeights)
      outHeights.push_back(h);
    if (hasNormals)
      outNormals.push_back(nrm);
    if (hasTrackIds)
      outTrackIds.push_back(trk);
    outFeatures.push_back(features);
    outRequired.push_back(required);
  };

  auto weightAt = [&](int id) {
    return (id < static_cast<int>(weights.size())) ? weights[id] : 0.04f;
  };
  auto heightAt = [&](int id) {
    return (hasHeights && id < static_cast<int>(heights.size())) ? heights[id]
                                                                 : 0.0f;
  };
  auto normalAt = [&](int id) {
    return (hasNormals && id < static_cast<int>(normals.size()))
               ? normals[id]
               : zVector(0, 0, 1);
  };
  auto trackAt = [&](int id) {
    return (hasTrackIds && id < static_cast<int>(trackIds.size()))
               ? trackIds[id]
               : 0.0f;
  };

  const int numSegments = closed ? N : (N - 1);

  // Emit the path's start vertex and reset the sampling phase there.
  const int firstId = sequence[0];
  appendPoint(positions[firstId], weightAt(firstId), heightAt(firstId),
              normalAt(firstId), trackAt(firstId), featureFlags[firstId], true);
  double distSinceSample = 0.0;

  for (int seg = 0; seg < numSegments; seg++) {
    int i0 = sequence[seg];
    int i1 = sequence[(seg + 1) % N];
    const zPoint &p0 = positions[i0];
    const zPoint &p1 = positions[i1];
    float w0 = weightAt(i0), w1 = weightAt(i1);
    float h0 = heightAt(i0), h1 = heightAt(i1);
    zVector n0 = normalAt(i0), n1 = normalAt(i1);
    float trk0 = trackAt(i0), trk1 = trackAt(i1);

    zVector segVec = vecSub(p1, p0);
    const double segLen = segVec.length();

    if (segLen > 0.0) {
      // March along this segment, continuing the running phase carried
      // over from previous segments, placing one sample every
      // sampleLength regardless of how long/short this segment is.
      double remaining = sampleLength - distSinceSample;
      while (remaining < segLen) {
        const double t = remaining / segLen;
        zPoint interpP(p0.x + t * segVec.x, p0.y + t * segVec.y,
                       p0.z + t * segVec.z);
        float interpW = static_cast<float>(w0 + t * (w1 - w0));
        float interpH = static_cast<float>(h0 + t * (h1 - h0));
        zVector interpN(n0.x + t * (n1.x - n0.x), n0.y + t * (n1.y - n0.y),
                        n0.z + t * (n1.z - n0.z));
        if (interpN.length() > 1e-6)
          interpN.normalize();
        else
          interpN = zVector(0, 0, 1);
        const float interpTrk =
            (trk0 + t * (trk1 - trk0)) >= 0.5f ? 1.0f : 0.0f;
        appendPoint(interpP, interpW, interpH, interpN, interpTrk, 0, false);
        remaining += sampleLength;
      }
      distSinceSample = segLen - (remaining - sampleLength);
    }

    // Force-emit i1 (and reset the phase there) if it's the open path's
    // true endpoint or a preserved feature/corner vertex.
    const bool isLastOpenVertex = (!closed && seg == numSegments - 1);
    const bool isFeature =
        (i1 < int(featureFlags.size())) && featureFlags[i1] != 0;
    if ((isLastOpenVertex || isFeature) && !(closed && i1 == firstId)) {
      appendPoint(p1, w1, h1, n1, trk1, featureFlags[i1], true);
      distSinceSample = 0.0;
    }
  }

  // Ordinary closing samples yield to the required start, including tiny
  // numerical residuals after a whole number of sample intervals. Required
  // closing features never participate in the old 0.1 mm closure weld.
  const double closingClearance = (std::min)(0.01, double(sampleLength) * 0.1);
  while (closed && outPoints.size() > 1 && !outRequired.back() &&
         vecDistSq(outPoints.front(), outPoints.back()) < closingClearance * closingClearance)
    popPoint();
  if (closed && outPoints.size() > 1 && vecDistSq(outPoints.front(), outPoints.back()) == 0.0) {
    outFeatures.front() |= outFeatures.back();
    popPoint();
  }
  const int numPts = static_cast<int>(outPoints.size());
  if (numPts >= 2) {
    for (int i = 0; i + 1 < numPts; i++) {
      outEdges.push_back(i);
      outEdges.push_back(i + 1);
    }
    if (closed && numPts >= 3) {
      const int finalN = static_cast<int>(outPoints.size());
      if (finalN >= 3) {
        outEdges.push_back(finalN - 1);
        outEdges.push_back(0);
      }
    }
  }

  if (outPoints.size() < (closed ? 3u : 2u) || outEdges.size() < (closed ? 6u : 2u))
    return false;

  fnGraph.clear();
  fnGraph.create(outPoints, outEdges);
  ioPath.widths = outWeights;
  ioPath.featureFlags = outFeatures;
  if (hasHeights)
    ioPath.heights = outHeights;
  if (hasNormals)
    ioPath.normals = outNormals;
  if (hasTrackIds)
    ioPath.trackIds = outTrackIds;
  return true;
}

// Smooths localized width glitches (isolated spikes/dips) in ioPath.widths.
void filterGlitches(PrintPath &ioPath) {
  const bool closed = isClosedGraph(ioPath.graph);
  zFloatArray &widths = ioPath.widths;
  const int n = static_cast<int>(widths.size());
  if (n <= 2)
    return;
  for (int i = 0; i < n; i++) {
    int prevId = (i == 0) ? (closed ? n - 1 : 0) : i - 1;
    int nextId = (i == n - 1) ? (closed ? 0 : n - 1) : i + 1;
    if (prevId == i || nextId == i)
      continue;
    float prevWidth = widths[prevId];
    float currentWidth = widths[i];
    float nextWidth = widths[nextId];
    if (prevWidth <= 0.0f || nextWidth <= 0.0f)
      continue;
    if (std::fabs(prevWidth - nextWidth) > 1.0f)
      continue;
    if (std::fabs(currentWidth - prevWidth) <= 1.0f)
      continue;
    widths[i] = prevWidth;
  }
}

void reversePrintPath(PrintPath &ioPath) {
  zFnGraph fnGraph(ioPath.graph);
  zPointArray positions;
  fnGraph.getVertexPositions(positions);
  const int numV = static_cast<int>(positions.size());
  if (numV < 2)
    return;

  const bool closed = isClosedGraph(ioPath.graph);

  std::reverse(positions.begin(), positions.end());
  if (!ioPath.widths.empty())
    std::reverse(ioPath.widths.begin(), ioPath.widths.end());
  if (!ioPath.heights.empty())
    std::reverse(ioPath.heights.begin(), ioPath.heights.end());
  if (!ioPath.normals.empty())
    std::reverse(ioPath.normals.begin(), ioPath.normals.end());
  if (!ioPath.trackIds.empty())
    std::reverse(ioPath.trackIds.begin(), ioPath.trackIds.end());
  if (!ioPath.featureFlags.empty())
    std::reverse(ioPath.featureFlags.begin(), ioPath.featureFlags.end());

  zIntArray newEdges;
  newEdges.reserve(numV * 2);
  for (int i = 0; i + 1 < numV; i++) {
    newEdges.push_back(i);
    newEdges.push_back(i + 1);
  }
  if (closed && numV >= 3) {
    newEdges.push_back(numV - 1);
    newEdges.push_back(0);
  }

  fnGraph.clear();
  fnGraph.create(positions, newEdges);
}

// Orients one layer's ring path in place - CCW if ccw is true, CW if false
// (extraction already produces at most one ring per layer, already
// inner-before-outer ordered).
bool orientLayer(std::vector<PrintPath> &ringPaths, bool ccw) {
  if (ringPaths.empty())
    return false;

  bool oriented = false;
  for (PrintPath &path : ringPaths) {
    if (!isClosedGraph(path.graph)) continue;
    zFnGraph fnGraph(path.graph);
    zPointArray positions;
    fnGraph.getVertexPositions(positions);
    if (positions.size() < 3)
      continue;
    const bool isCCW = signedArea2D(positions) >= 0.0;
    if (isCCW != ccw) {
      reversePrintPath(path);
    }
    oriented = true;
  }
  return oriented;
}

// ---- alignStartPoints helpers ----

std::vector<double> cumulativeArcLength(const zPointArray &positions) {
  const int n = static_cast<int>(positions.size());
  std::vector<double> cumLen(n + 1, 0.0);
  for (int i = 0; i < n; i++) {
    cumLen[i + 1] =
        cumLen[i] + vecSub(positions[(i + 1) % n], positions[i]).length();
  }
  return cumLen;
}

// Arc length of positions' inner-tagged prefix (index 0 is always inner
// and inner is a contiguous prefix, per extractTrackPaths), or the whole
// loop's length if trackIds carries no track distinction at all.
double innerSegmentLength(const zPointArray &positions,
                          const zFloatArray &trackIds,
                          const std::vector<double> &cumLen) {
  const int n = static_cast<int>(positions.size());
  int innerEnd = n;
  if (!trackIds.empty()) {
    if (std::none_of(trackIds.begin(), trackIds.end(), [](float id) { return id < 0.5f; }))
      return cumLen[n];
    innerEnd = 0;
    while (innerEnd < n && innerEnd < int(trackIds.size()) &&
           trackIds[innerEnd] < 0.5f)
      innerEnd++;
  }
  return cumLen[innerEnd];
}

// ---- failure diagnostics ----

// One-line topology report of a print path for failure messages: vertex /
// edge / component counts, whether the closing edge v0-v(N-1) exists, and -
// for an open chain - where its ends are (index, position, inner/outer
// track) plus the gap between the first two ends, with a hint at where the
// chain broke.
std::string describePathTopology(const PrintPath &path) {
  zFnGraph fn(const_cast<zObjGraph &>(path.graph));
  const int n = fn.numVertices();
  zIntArray edges;
  fn.getEdgeData(edges);
  zPointArray positions;
  fn.getVertexPositions(positions);

  std::vector<int> valence(n, 0), parent(n);
  for (int i = 0; i < n; i++)
    parent[i] = i;
  const auto findRoot = [&](int v) {
    while (parent[v] != v)
      v = parent[v] = parent[parent[v]];
    return v;
  };
  for (size_t e = 0; e + 1 < edges.size(); e += 2) {
    const int a = edges[e], b = edges[e + 1];
    if (a < 0 || b < 0 || a >= n || b >= n)
      continue;
    valence[a]++;
    valence[b]++;
    parent[findRoot(a)] = findRoot(b);
  }
  int components = 0, isolated = 0, branches = 0;
  std::vector<int> ends;
  for (int i = 0; i < n; i++) {
    if (findRoot(i) == i)
      components++;
    if (valence[i] == 0)
      isolated++;
    else if (valence[i] == 1)
      ends.push_back(i);
    else if (valence[i] > 2)
      branches++;
  }

  const auto trackOf = [&](int v) -> const char * {
    if (v >= int(path.trackIds.size()))
      return "no track";
    return path.trackIds[v] >= 0.5f ? "outer" : "inner";
  };

  std::ostringstream out;
  out << std::setprecision(6) << "vertices=" << n
      << " edges=" << edges.size() / 2 << " components=" << components
      << " open ends=" << ends.size() << " branch vertices=" << branches
      << " isolated vertices=" << isolated << ", closing edge v0-v"
      << (n - 1) << (isClosedGraph(path.graph) ? " present" : " missing");
  if (!path.trackIds.empty()) {
    int outer = 0;
    for (float t : path.trackIds)
      outer += (t >= 0.5f) ? 1 : 0;
    out << ", track vertices inner=" << int(path.trackIds.size()) - outer
        << " outer=" << outer;
  }
  for (size_t k = 0; k < ends.size() && k < 4; k++) {
    const zPoint &p = positions[ends[k]];
    out << " | end v" << ends[k] << " at (" << p.x << ", " << p.y << ", "
        << p.z << ", " << trackOf(ends[k]) << ")";
  }
  if (ends.size() > 4)
    out << " | ... " << ends.size() - 4 << " more ends";
  if (ends.size() >= 2) {
    const zPoint &a = positions[ends[0]];
    const zPoint &b = positions[ends[1]];
    const double gap = std::sqrt((a.x - b.x) * (a.x - b.x) +
                                 (a.y - b.y) * (a.y - b.y) +
                                 (a.z - b.z) * (a.z - b.z));
    out << " | gap between the first two ends=" << gap << " | hint: ";
    if (components > 1)
      out << "the path is split into " << components
          << " pieces - the zero-isocontour is broken (grid too coarse for "
             "the bead width, or cut by the field grid's edge)";
    else if (std::string(trackOf(ends[0])) != trackOf(ends[1]))
      out << "the chain breaks where the inner and outer tracks join "
             "(trim slot / double-boundary junction)";
    else
      out << "the chain breaks inside the " << trackOf(ends[0])
          << " track";
  }
  return out.str();
}

// Layers whose first path (the one alignStartPoints() rotates) has >= 3
// vertices but no closing edge.
std::vector<int> openLeadPathLayers(const PrintPathArray &layers) {
  std::vector<int> open;
  for (int l = 0; l < int(layers.size()); l++) {
    if (layers[l].empty())
      continue;
    zFnGraph fn(const_cast<zObjGraph &>(layers[l][0].graph));
    if (fn.numVertices() >= 3 && !isClosedGraph(layers[l][0].graph))
      open.push_back(l);
  }
  return open;
}

std::string layerList(const std::vector<int> &layers, size_t maxShown = 12) {
  std::ostringstream out;
  for (size_t k = 0; k < layers.size() && k < maxShown; k++)
    out << (k ? ", " : "") << layers[k];
  if (layers.size() > maxShown)
    out << ", ... (" << layers.size() - maxShown << " more)";
  return out.str();
}

bool insertAndRotateToStart(PrintPath &ioPath, double s,
                            std::string *outReason = nullptr,
                            double weldEps = 1e-3) {
  zFnGraph fnGraph(ioPath.graph);
  zPointArray positions;
  fnGraph.getVertexPositions(positions);
  const int n = static_cast<int>(positions.size());
  if (n < 3) {
    if (outReason)
      *outReason = "path has only " + std::to_string(n) + " vertices (needs >= 3)";
    return false;
  }
  if (!isClosedGraph(ioPath.graph)) {
    if (outReason)
      *outReason = "path is not a closed ring, so its start cannot be "
                   "rotated: " +
                   describePathTopology(ioPath);
    return false;
  }

  const std::vector<double> cumLen = cumulativeArcLength(positions);
  const double L = cumLen[n];
  if (L <= 1e-6) {
    if (outReason)
      *outReason = "path has zero total length (" + std::to_string(n) +
                   " coincident vertices)";
    return false;
  }
  s = std::fmod(s, L);
  if (s < 0.0)
    s += L;

  int seg = 0;
  for (; seg < n; seg++) {
    if (s <= cumLen[seg + 1] || seg == n - 1)
      break;
  }
  const int i0 = seg, i1 = (seg + 1) % n;
  const double segStart = cumLen[seg], segLen = cumLen[seg + 1] - segStart;
  const double t =
      (segLen > 1e-6)
          ? (std::max)(0.0, (std::min)(1.0, (s - segStart) / segLen))
          : 0.0;

  const bool hasHeights = !ioPath.heights.empty();
  const bool hasNormals = !ioPath.normals.empty();
  const bool hasTrackIds = !ioPath.trackIds.empty();

  zPointArray newPositions = positions;
  zFloatArray newWidths = ioPath.widths;
  zFloatArray newHeights = ioPath.heights;
  zVectorArray newNormals = ioPath.normals;
  zFloatArray newTrackIds = ioPath.trackIds;
  zIntArray newFeatures = ioPath.featureFlags;
  int startIdx;

  const double tEps = weldEps / (std::max)(segLen, 1e-6);
  if (t <= tEps) {
    startIdx = i0;
  } else if (t >= 1.0 - tEps) {
    startIdx = i1;
  } else {
    const zPoint &p0 = positions[i0];
    const zPoint &p1 = positions[i1];
    zPoint newP(p0.x + t * (p1.x - p0.x), p0.y + t * (p1.y - p0.y),
                p0.z + t * (p1.z - p0.z));

    const float w0 =
        (i0 < int(ioPath.widths.size())) ? ioPath.widths[i0] : 0.04f;
    const float w1 =
        (i1 < int(ioPath.widths.size())) ? ioPath.widths[i1] : 0.04f;
    const float newW = static_cast<float>(w0 + t * (w1 - w0));

    float newH = 0.0f;
    if (hasHeights) {
      const float h0 =
          (i0 < int(ioPath.heights.size())) ? ioPath.heights[i0] : 0.0f;
      const float h1 =
          (i1 < int(ioPath.heights.size())) ? ioPath.heights[i1] : 0.0f;
      newH = static_cast<float>(h0 + t * (h1 - h0));
    }
    zVector newN(0, 0, 1);
    if (hasNormals) {
      const zVector n0 = (i0 < int(ioPath.normals.size())) ? ioPath.normals[i0]
                                                           : zVector(0, 0, 1);
      const zVector n1 = (i1 < int(ioPath.normals.size())) ? ioPath.normals[i1]
                                                           : zVector(0, 0, 1);
      zVector n(n0.x + t * (n1.x - n0.x), n0.y + t * (n1.y - n0.y),
                n0.z + t * (n1.z - n0.z));
      if (n.length() > 1e-6)
        n.normalize();
      else
        n = zVector(0, 0, 1);
      newN = n;
    }
    float newTrk = 0.0f;
    if (hasTrackIds) {
      const float trk0 =
          (i0 < int(ioPath.trackIds.size())) ? ioPath.trackIds[i0] : 0.0f;
      const float trk1 =
          (i1 < int(ioPath.trackIds.size())) ? ioPath.trackIds[i1] : 0.0f;
      newTrk = (trk0 + t * (trk1 - trk0)) >= 0.5f ? 1.0f : 0.0f;
    }

    newPositions.insert(newPositions.begin() + i1, newP);
    newWidths.insert(newWidths.begin() + i1, newW);
    if (hasHeights)
      newHeights.insert(newHeights.begin() + i1, newH);
    if (hasNormals)
      newNormals.insert(newNormals.begin() + i1, newN);
    if (hasTrackIds)
      newTrackIds.insert(newTrackIds.begin() + i1, newTrk);
    if (!newFeatures.empty()) newFeatures.insert(newFeatures.begin() + i1, 0);
    startIdx = i1;
  }

  std::rotate(newPositions.begin(), newPositions.begin() + startIdx,
              newPositions.end());
  std::rotate(newWidths.begin(), newWidths.begin() + startIdx, newWidths.end());
  if (hasHeights)
    std::rotate(newHeights.begin(), newHeights.begin() + startIdx,
                newHeights.end());
  if (hasNormals)
    std::rotate(newNormals.begin(), newNormals.begin() + startIdx,
                newNormals.end());
  if (hasTrackIds)
    std::rotate(newTrackIds.begin(), newTrackIds.begin() + startIdx,
                newTrackIds.end());
  if (!newFeatures.empty())
    std::rotate(newFeatures.begin(), newFeatures.begin() + startIdx, newFeatures.end());

  const int newN = static_cast<int>(newPositions.size());
  zIntArray newEdges;
  newEdges.reserve(newN * 2);
  for (int i = 0; i + 1 < newN; i++) {
    newEdges.push_back(i);
    newEdges.push_back(i + 1);
  }
  newEdges.push_back(newN - 1);
  newEdges.push_back(0);

  fnGraph.clear();
  fnGraph.create(newPositions, newEdges);
  ioPath.widths = newWidths;
  ioPath.featureFlags = newFeatures;
  if (hasHeights)
    ioPath.heights = newHeights;
  if (hasNormals)
    ioPath.normals = newNormals;
  if (hasTrackIds)
    ioPath.trackIds = newTrackIds;
  return true;
}

// Finds the arc-length position of the point on path's closed loop nearest to
// target
bool nearestInnerArcPos(const PrintPath &path, const zPoint &target,
                        double &outArcPos) {
  zFnGraph fnGraph(const_cast<zObjGraph &>(path.graph));
  zPointArray positions;
  fnGraph.getVertexPositions(positions);
  const int n = static_cast<int>(positions.size());
  if (n < 3 || !isClosedGraph(path.graph))
    return false;

  const bool hasTrackIds = !path.trackIds.empty();
  const bool hasInnerTrack = std::any_of(path.trackIds.begin(), path.trackIds.end(),
                                       [](float id) { return id < 0.5f; });
  auto isOuter = [&](int id) {
    return hasTrackIds && id < int(path.trackIds.size()) &&
           path.trackIds[id] >= 0.5f;
  };

  const std::vector<double> cumLen = cumulativeArcLength(positions);

  double bestDistSq = (std::numeric_limits<double>::max)();
  double bestArcPos = 0.0;
  bool found = false;

  for (int i = 0; i < n; i++) {
    const int i0 = i, i1 = (i + 1) % n;
    if (hasInnerTrack && (isOuter(i0) || isOuter(i1)))
      continue; // keep the search strictly inside the inner segment

    const zPoint &p0 = positions[i0];
    const zPoint &p1 = positions[i1];
    zVector seg = vecSub(p1, p0);
    const double segLenSq = vecDot(seg, seg);
    double t = 0.0;
    if (segLenSq > 1e-6) {
      t = vecDot(vecSub(target, p0), seg) / segLenSq;
      t = (std::max)(0.0, (std::min)(1.0, t));
    }
    zPoint proj(p0.x + t * seg.x, p0.y + t * seg.y, p0.z + t * seg.z);
    const double distSq = vecDistSq(target, proj);
    if (distSq < bestDistSq) {
      bestDistSq = distSq;
      bestArcPos = cumLen[i] + t * (cumLen[i + 1] - cumLen[i]);
      found = true;
    }
  }

  if (found)
    outArcPos = bestArcPos;
  return found;
}

// Drops a closed loop's wraparound edge (last vertex -> first vertex), so
// the seam sits exactly at the loop's current start/end instead of being
// welded shut. Point order and count are unchanged - only the closing
// edge goes away. No-op if the path isn't closed.
void openLoopAtStart(PrintPath &ioPath) {
  if (!isClosedGraph(ioPath.graph))
    return;
  zFnGraph fnGraph(ioPath.graph);
  zPointArray positions;
  fnGraph.getVertexPositions(positions);
  const int n = static_cast<int>(positions.size());
  if (n < 2)
    return;

  zIntArray newEdges;
  newEdges.reserve((n - 1) * 2);
  for (int i = 0; i + 1 < n; i++) {
    newEdges.push_back(i);
    newEdges.push_back(i + 1);
  }

  fnGraph.clear();
  fnGraph.create(positions, newEdges);
}

// Cuts the closed ring open right before its (already finalized) start
// point
void cutGapBeforeStart(PrintPath &ioPath, double seamSize) {
  if (seamSize <= 1e-6) {
    openLoopAtStart(ioPath);
    return;
  }
  if (!isClosedGraph(ioPath.graph))
    return;

  zFnGraph fnGraph(ioPath.graph);
  zPointArray positions;
  fnGraph.getVertexPositions(positions);
  const int n = static_cast<int>(positions.size());
  if (n < 3)
    return;

  const std::vector<double> cumLen = cumulativeArcLength(positions);
  const double L = cumLen[n];
  if (L <= 1e-6)
    return;

  const double gap =
      (std::min)(seamSize, L * 0.5); // never eat more than half the loop
  const double targetArcPos = std::fmod(L - gap, L);

  int seg = n - 1;
  for (int i = 0; i < n; i++) {
    if (targetArcPos <= cumLen[i + 1]) {
      seg = i;
      break;
    }
  }
  const int i0 = seg, i1 = (seg + 1) % n;
  const double segStart = cumLen[seg], segLen = cumLen[seg + 1] - segStart;
  const double t =
      (segLen > 1e-6)
          ? (std::max)(0.0, (std::min)(1.0, (targetArcPos - segStart) / segLen))
          : 0.0;

  const bool hasHeights = !ioPath.heights.empty();
  const bool hasNormals = !ioPath.normals.empty();
  const bool hasTrackIds = !ioPath.trackIds.empty();

  zPointArray keptPositions;
  zFloatArray keptWidths;
  zFloatArray keptHeights;
  zVectorArray keptNormals;
  zFloatArray keptTrackIds;
  zIntArray keptFeatures;

  auto appendIdx = [&](int id) {
    keptPositions.push_back(positions[id]);
    keptWidths.push_back(id < int(ioPath.widths.size()) ? ioPath.widths[id]
                                                        : 0.04f);
    if (hasHeights)
      keptHeights.push_back(id < int(ioPath.heights.size()) ? ioPath.heights[id]
                                                            : 0.0f);
    if (hasNormals)
      keptNormals.push_back(id < int(ioPath.normals.size()) ? ioPath.normals[id]
                                                            : zVector(0, 0, 1));
    if (hasTrackIds)
      keptTrackIds.push_back(
          id < int(ioPath.trackIds.size()) ? ioPath.trackIds[id] : 0.0f);
    if (!ioPath.featureFlags.empty()) keptFeatures.push_back(ioPath.featureFlags[id]);
  };

  // Keep the surviving arc from the start (index 0) through i0, then
  // append the new end point (interpolated, or i1/i0 if it snaps onto an
  // existing vertex) - producing an open chain from start to this new
  // end, with the small [end, start) gap dropped.
  for (int id = 0; id <= i0; id++)
    appendIdx(id);

  const double tEps = 1e-3 / (std::max)(segLen, 1e-6);
  if (t >= 1.0 - tEps) {
    if (i1 != i0)
      appendIdx(i1);
  } else if (t > tEps) {
    const zPoint &p0 = positions[i0];
    const zPoint &p1 = positions[i1];
    zPoint newP(p0.x + t * (p1.x - p0.x), p0.y + t * (p1.y - p0.y),
                p0.z + t * (p1.z - p0.z));
    const float w0 =
        (i0 < int(ioPath.widths.size())) ? ioPath.widths[i0] : 0.04f;
    const float w1 =
        (i1 < int(ioPath.widths.size())) ? ioPath.widths[i1] : 0.04f;
    keptPositions.push_back(newP);
    keptWidths.push_back(static_cast<float>(w0 + t * (w1 - w0)));
    if (!ioPath.featureFlags.empty()) keptFeatures.push_back(0);
    if (hasHeights) {
      const float h0 =
          (i0 < int(ioPath.heights.size())) ? ioPath.heights[i0] : 0.0f;
      const float h1 =
          (i1 < int(ioPath.heights.size())) ? ioPath.heights[i1] : 0.0f;
      keptHeights.push_back(static_cast<float>(h0 + t * (h1 - h0)));
    }
    if (hasNormals) {
      const zVector n0 = (i0 < int(ioPath.normals.size())) ? ioPath.normals[i0]
                                                           : zVector(0, 0, 1);
      const zVector n1 = (i1 < int(ioPath.normals.size())) ? ioPath.normals[i1]
                                                           : zVector(0, 0, 1);
      zVector nrm(n0.x + t * (n1.x - n0.x), n0.y + t * (n1.y - n0.y),
                  n0.z + t * (n1.z - n0.z));
      if (nrm.length() > 1e-6)
        nrm.normalize();
      else
        nrm = zVector(0, 0, 1);
      keptNormals.push_back(nrm);
    }
    if (hasTrackIds) {
      const float trk0 =
          (i0 < int(ioPath.trackIds.size())) ? ioPath.trackIds[i0] : 0.0f;
      const float trk1 =
          (i1 < int(ioPath.trackIds.size())) ? ioPath.trackIds[i1] : 0.0f;
      keptTrackIds.push_back((trk0 + t * (trk1 - trk0)) >= 0.5f ? 1.0f : 0.0f);
    }
  }
  // t <= tEps: targetArcPos already snapped onto i0, which is already the
  // last appended point - nothing more to add.

  const int keptN = static_cast<int>(keptPositions.size());
  if (keptN < 2)
    return;

  zIntArray newEdges;
  newEdges.reserve((keptN - 1) * 2);
  for (int i = 0; i + 1 < keptN; i++) {
    newEdges.push_back(i);
    newEdges.push_back(i + 1);
  }

  fnGraph.clear();
  fnGraph.create(keptPositions, newEdges);
  ioPath.widths = keptWidths;
  const auto countFeatures = [](const zIntArray &flags) {
    return std::count_if(flags.begin(), flags.end(), [](int flag) { return flag != 0; });
  };
  const auto removedFeatures = countFeatures(ioPath.featureFlags) - countFeatures(keptFeatures);
  if (removedFeatures > 0)
    zTsUtils::printDebug("[cutGapBeforeStart] Explicit seam gap removed " +
        std::to_string(removedFeatures) + " feature vertices from the cut arc.");
  ioPath.featureFlags = keptFeatures;
  if (hasHeights)
    ioPath.heights = keptHeights;
  if (hasNormals)
    ioPath.normals = keptNormals;
  if (hasTrackIds)
    ioPath.trackIds = keptTrackIds;
}

} // namespace

ZSPACE_INLINE zIntArray zPrintSyn::cornerVertexIds(const PrintPath &path, float angleThreshold) {
  zPointArray points;
  zFnGraph(const_cast<zObjGraph &>(path.graph)).getVertexPositions(points);
  std::vector<GraphPolyline> curves;
  zIntArray result;
  if (!graphPolylines(path.graph, curves)) return result;
  for (const auto &curve : curves) {
    std::vector<int> flags;
    detectFeatureAngles(points, curve.Vertices, angleThreshold, curve.Closed, flags);
    for (size_t i = 0; i < curve.Vertices.size(); ++i) {
      if (!curve.Closed && (i == 0 || i + 1 == curve.Vertices.size())) continue;
      if (flags[curve.Vertices[i]]) result.push_back(curve.Vertices[i]);
    }
  }
  return result;
}

ZSPACE_INLINE bool zPrintSyn::reconstructCorners(PrintPath &path, float angleThreshold) {
  std::vector<GraphPolyline> curves;
  if (!graphPolylines(path.graph, curves) || curves.size() != 1) return false;
  zPointArray points;
  zFnGraph fn(path.graph);
  fn.getVertexPositions(points);
  auto sequence = curves[0].Vertices;
  const int n = static_cast<int>(sequence.size());
  const bool closed = curves[0].Closed;
  if (n != static_cast<int>(points.size())) return false;
  for (const auto &p : points)
    if (!std::isfinite(p.x) || !std::isfinite(p.y) || !std::isfinite(p.z)) return false;
  if ((!path.widths.empty() && path.widths.size() != points.size()) ||
      (!path.heights.empty() && path.heights.size() != points.size()) ||
      (!path.normals.empty() && path.normals.size() != points.size()) ||
      (!path.trackIds.empty() && path.trackIds.size() != points.size()) ||
      (!path.featureFlags.empty() && path.featureFlags.size() != points.size())) return false;
  std::vector<int> flags;
  detectFeatureAngles(points, sequence, angleThreshold, closed, flags);
  if (!closed) flags[sequence.front()] = flags[sequence.back()] = 0;
  const zIntArray previousFeatures = path.featureFlags;
  zIntArray recorded = previousFeatures;
  if (recorded.empty()) recorded.assign(points.size(), 0);
  for (int id = 0; id < n; ++id) {
    recorded[id] |= flags[id];
    // Features already reconstructed/recorded are anchors. A subsequent
    // resampling call must not reconstruct them into a different corner.
    if (!previousFeatures.empty() && previousFeatures[id]) flags[id] = 0;
  }
  path.featureFlags = recorded;
  if (n < 5) return true;
  auto merged = flags;
  for (int i = 0; i < n; ++i) {
    if (!closed && (i == 0 || i == n-1)) continue;
    if (flags[sequence[(i+n-1)%n]] && flags[sequence[(i+1)%n]] &&
        (previousFeatures.empty() || !previousFeatures[sequence[i]]))
      merged[sequence[i]] = 1;
  }
  if (closed) {
    // Cut the circular sequence between two support vertices, never in a group.
    int cut = -1;
    for (int i = 0; i < n; ++i)
      if (!merged[sequence[i]] && !merged[sequence[(i+n-1)%n]]) { cut = i; break; }
    if (cut < 0) return true;
    std::rotate(sequence.begin(), sequence.begin()+cut, sequence.end());
  }
  zPointArray output;
  zIntArray sources;
  zIntArray outputFeatures;
  for (int i = 0; i < n;) {
    if (!merged[sequence[i]]) {
      outputFeatures.push_back(recorded[sequence[i]]);
      sources.push_back(sequence[i]); output.push_back(points[sequence[i++]]); continue;
    }
    const int first = i;
    while (i < n && merged[sequence[i]]) ++i;
    const int last = i-1;
    bool usable = closed || (first >= 2 && last+2 < n);
    const int a = (first+n-2)%n, b = (first+n-1)%n;
    const int c = (last+1)%n, d = (last+2)%n;
    usable = usable && !merged[sequence[a]] && !merged[sequence[b]] &&
        !merged[sequence[c]] && !merged[sequence[d]];
    const auto &p = points[sequence[b]], &q = points[sequence[c]];
    const auto &pa = points[sequence[a]], &qd = points[sequence[d]];
    const double ux = double(p.x)-pa.x, uy = double(p.y)-pa.y;
    const double vx = double(qd.x)-q.x, vy = double(qd.y)-q.y;
    const double determinant = ux*vy-uy*vx;
    const double scale = std::sqrt((ux*ux+uy*uy)*(vx*vx+vy*vy));
    usable = usable && scale > 0 && std::abs(determinant) > 1e-8*scale;
    double x = 0, y = 0;
    if (usable) {
      const double t = ((double(q.x)-p.x)*vy-(double(q.y)-p.y)*vx)/determinant;
      const double u = ((double(q.x)-p.x)*uy-(double(q.y)-p.y)*ux)/determinant;
      x = p.x+t*ux; y = p.y+t*uy;
      // Incoming line extends forwards; outgoing line extends backwards.
      usable = std::isfinite(x) && std::isfinite(y) && t >= 0 && u <= 0;
    }
    if (usable) {
      const int representative = sequence[first];
      output.push_back(zPoint(x, y, points[representative].z));
      sources.push_back(representative);
      int features = PrintPath::IntersectionFeature;
      for (int k = first; k <= last; ++k) features |= recorded[sequence[k]];
      outputFeatures.push_back(features);
    } else {
      for (int k = first; k <= last; ++k) {
        sources.push_back(sequence[k]); output.push_back(points[sequence[k]]);
        outputFeatures.push_back(recorded[sequence[k]]);
      }
    }
  }
  if (output.size() < (closed ? 3u : 2u)) return true;
  // Keep the original starting vertex when it survives, or its merged corner.
  const auto start = std::find(sources.begin(), sources.end(), curves[0].Vertices.front());
  if (closed && start != sources.end()) {
    const auto offset = start-sources.begin();
    std::rotate(sources.begin(), start, sources.end());
    std::rotate(output.begin(), output.begin()+offset, output.end());
    std::rotate(outputFeatures.begin(), outputFeatures.begin()+offset, outputFeatures.end());
  }
  PrintPath updated;
  updated.featureFlags = outputFeatures;
  for (int id : sources) {
    if (!path.widths.empty()) updated.widths.push_back(path.widths[id]);
    if (!path.heights.empty()) updated.heights.push_back(path.heights[id]);
    if (!path.normals.empty()) updated.normals.push_back(path.normals[id]);
    if (!path.trackIds.empty()) updated.trackIds.push_back(path.trackIds[id]);
  }
  zIntArray edges;
  for (int i = 1; i < static_cast<int>(output.size()); ++i) { edges.push_back(i-1); edges.push_back(i); }
  if (closed) { edges.push_back(static_cast<int>(output.size())-1); edges.push_back(0); }
  zFnGraph(updated.graph).create(output, edges);
  path = std::move(updated);
  return true;
}


ZSPACE_INLINE zIntArray zPrintSyn::featureVertexIds(const PrintPath &path,float fallbackAngleThreshold) {
  if(path.featureFlags.empty())return cornerVertexIds(path,fallbackAngleThreshold);
  if(path.featureFlags.size()!=size_t(zFnGraph(const_cast<zObjGraph&>(path.graph)).numVertices()))return {};
  zIntArray ids;for(size_t i=0;i<path.featureFlags.size();++i)
    if(path.featureFlags[i]&(PrintPath::AngleFeature|PrintPath::IntersectionFeature))ids.push_back(int(i));
  return ids;
}


//--------------------------
//---- PRIVATE METHODS
//--------------------------

// The two align-chart overloads: alignCharts[l] is layer l's align chart
// (entries may repeat).
ZSPACE_INLINE bool zPrintSyn::alignStartPointsOnCharts(PrintPathArray &ioPaths,
                              const std::vector<const zObjMesh *> &alignCharts,
                              double seamOffset, double seamSize,
                              std::string &outStatus,
                              std::vector<int> *outBadLayers, bool useParallel) {
  if (outBadLayers)
    outBadLayers->clear();
  const int totalLayers = int(ioPaths.size());
  // Every layer's align chart shares one topology, so the seed layer's
  // barycentric start point can be realized on any other layer's chart.
  bool chartsValid =
      totalLayers > 0 && int(alignCharts.size()) >= totalLayers;
  for (int l = 0; chartsValid && l < totalLayers; l++)
    chartsValid = alignCharts[l] != nullptr;
  if (!chartsValid) {
    outStatus = "failed: align charts missing or size mismatch (totalLayers=" +
                std::to_string(totalLayers) + ", alignCharts2D=" +
                std::to_string(int(alignCharts.size())) + ").";
    zTsUtils::printDebug("[alignStartPoints] Alignment charts are missing or do not match the path layers.");
    return false;
  }

  // The seed's barycentric position on its own align chart, realized on
  // every layer's align chart (the same transfer both slicer operators use
  // for the boundary slot anchor).
  return alignStartPointsWith(
      ioPaths,
      [&, useParallel](int seedLayer, const zPoint &seed, zPointArray &outMapped,
          std::vector<char> &outValid) {
        return mapPointToMeshes(
            seed, const_cast<zObjMesh &>(*alignCharts[seedLayer]),
            alignCharts, outMapped, outValid, nullptr, useParallel);
      },
      seamOffset, seamSize, outStatus, outBadLayers, useParallel);
}

// Shared body of every alignStartPoints() overload; mapSeed carries the seed
// layer's start point to every layer.
ZSPACE_INLINE bool zPrintSyn::alignStartPointsWith(PrintPathArray &ioPaths,
                                                   const SeedMapper &mapSeed,
                                                   double seamOffset,
                                                   double seamSize,
                                                   std::string &outStatus,
                                                   std::vector<int> *outBadLayers, bool useParallel) {
  if (outBadLayers)
    outBadLayers->clear();
  const int totalLayers = int(ioPaths.size());
  for (int l = 0; l < totalLayers; ++l) for (const auto &path : ioPaths[l]) {
    if (!path.featureFlags.empty() && path.featureFlags.size() !=
        size_t(zFnGraph(const_cast<zObjGraph &>(path.graph)).numVertices())) {
      outStatus = "failed: feature count does not match path vertices on layer " + std::to_string(l);
      if (outBadLayers) outBadLayers->push_back(l);
      return false;
    }
  }
  int seedLayer = -1;
  for (int l = 0; l < totalLayers; l++) {
    if (!ioPaths[l].empty()) {
      zFnGraph fnGraph(ioPaths[l][0].graph);
      if (fnGraph.numVertices() >= 3) {
        seedLayer = l;
        break;
      }
    }
  }
  if (seedLayer < 0) {
    outStatus = "failed: no layer has a print path with >= 3 vertices (all " +
                std::to_string(totalLayers) +
                " layers are empty or degenerate).";
    zTsUtils::printDebug("[alignStartPoints] All paths are empty or degenerate; no usable seed ring exists.");
    return false;
  }

  // 1. Seed: walk seamOffset from the current start along the inner
  // segment (index 0 is always inner, as a contiguous prefix, per
  // extractTrackPaths), clamped so it can't cross into outer.
  PrintPath &seedPath = ioPaths[seedLayer][0];
  {
    zFnGraph fnGraph(seedPath.graph);
    zPointArray positions;
    fnGraph.getVertexPositions(positions);
    const std::vector<double> cumLen = cumulativeArcLength(positions);
    const double innerLen =
        innerSegmentLength(positions, seedPath.trackIds, cumLen);
    const double clampedWalk =
        (std::max)(0.0,
                   (std::min)(seamOffset, (std::max)(0.0, innerLen - 0.1)));
    std::string reason;
    if (!insertAndRotateToStart(seedPath, clampedWalk, &reason)) {
      outStatus = "failed: cannot place the seam start on seed layer " +
                  std::to_string(seedLayer) + " (the first layer with a "
                  ">= 3-vertex path; innerLen=" +
                  std::to_string(innerLen) + ", clampedWalk=" +
                  std::to_string(clampedWalk) + "): " + reason;
      zTsUtils::printDebug("[alignStartPoints] The seed seam point could not be inserted and ordered along the inner path.");
      if (outBadLayers)
        outBadLayers->assign(1, seedLayer);
      return false;
    }
  }

  zPoint seedPoint;
  {
    zFnGraph fnGraph(seedPath.graph);
    zPointArray positions;
    fnGraph.getVertexPositions(positions);
    if (positions.empty()) {
      outStatus = "failed: seed path has no vertices after "
                  "insertAndRotateToStart.";
      zTsUtils::printDebug("[alignStartPoints] The seed path has no vertices after inserting and rotating its start.");
      if (outBadLayers)
        outBadLayers->assign(1, seedLayer);
      return false;
    }
    seedPoint = positions[0];
  }

  // 2. The seed point carried to every layer (mapSeed).
  zPointArray mappedSeeds;
  std::vector<char> mappedValid;
  if (!mapSeed(seedLayer, seedPoint, mappedSeeds, mappedValid) ||
      int(mappedValid.size()) < totalLayers ||
      int(mappedSeeds.size()) < totalLayers) {
    outStatus = "failed: could not carry the seed point of layer " +
                std::to_string(seedLayer) + " to the other layers.";
    zTsUtils::printDebug("[alignStartPoints] The seed anchor could not be transferred to the other layer charts.");
    if (outBadLayers)
      outBadLayers->assign(1, seedLayer);
    return false;
  }

  // 3. Every layer: rotate its ring to start at the nearest inner-segment
  // point to its mapped seed. No sequential dependency between layers (the
  // target is the fixed seed, not the previous layer's result), so this is
  // safe to run in parallel.
  int alignedCount = 0;
  int skippedEmpty = 0, skippedDegenerate = 0, skippedRealize = 0,
      skippedNearest = 0, skippedRotate = 0;
  // Lowest failing layer per kind, for the status line.
  int firstRealizeLayer = -1, firstNearestLayer = -1, firstRotateLayer = -1;
  std::string firstRotateReason;
  std::vector<int> unaligned; // every layer keeping its unaligned start
  const auto noteFirst = [&](int &first, int l) {
    if (first < 0 || l < first)
      first = l;
    unaligned.push_back(l);
  };
  zTsDetail::ParallelErrors parallelErrors1(totalLayers);
#pragma omp parallel for schedule(dynamic) if(useParallel)
  for (int l = 0; l < totalLayers; l++) {
    try {
      if (ioPaths[l].empty()) {
#pragma omp atomic
        skippedEmpty++;
        continue;
      }
      for (PrintPath &path : ioPaths[l]) {
        zFnGraph fnGraph(path.graph);
        if (fnGraph.numVertices() < 3) {
  #pragma omp atomic
          skippedDegenerate++;
          continue;
        }
        if (!isClosedGraph(path.graph)) continue;

        if (!mappedValid[l]) {
  #pragma omp atomic
          skippedRealize++;
  #pragma omp critical(zPrintSynAlignFirst)
          noteFirst(firstRealizeLayer, l);
          continue;
        }

        double arcPos = 0.0;
        if (!nearestInnerArcPos(path, mappedSeeds[l], arcPos)) {
  #pragma omp atomic
          skippedNearest++;
  #pragma omp critical(zPrintSynAlignFirst)
          noteFirst(firstNearestLayer, l);
          continue;
        }

        std::string reason;
        if (insertAndRotateToStart(path, arcPos, &reason)) {
  #pragma omp atomic
          alignedCount++;
        } else {
  #pragma omp atomic
          skippedRotate++;
  #pragma omp critical(zPrintSynAlignFirst)
          {
            if (firstRotateLayer < 0 || l < firstRotateLayer) {
              firstRotateLayer = l;
              firstRotateReason = reason;
            }
            unaligned.push_back(l);
          }
        }
      }
    } catch (...) { parallelErrors1.capture(l); }
  }
  parallelErrors1.rethrow();

  zTsDetail::ParallelErrors parallelErrors2(totalLayers);
#pragma omp parallel for schedule(dynamic) if(useParallel)
  for (int l = 0; l < totalLayers; l++) {
    try {
      if (ioPaths[l].empty())
        continue;
      for (PrintPath &path : ioPaths[l])
        if (isClosedGraph(path.graph)) cutGapBeforeStart(path, seamSize);
    } catch (...) { parallelErrors2.capture(l); }
  }
  parallelErrors2.rethrow();

  const std::string counts =
      " layers (empty=" + std::to_string(skippedEmpty) +
      ", degenerate=" + std::to_string(skippedDegenerate) +
      ", realizeFailed=" + std::to_string(skippedRealize) +
      ", nearestFailed=" + std::to_string(skippedNearest) +
      ", rotateFailed=" + std::to_string(skippedRotate) + ")";
  std::sort(unaligned.begin(), unaligned.end());
  unaligned.erase(std::unique(unaligned.begin(), unaligned.end()), unaligned.end());
  if (outBadLayers)
    *outBadLayers = unaligned;
  std::string details;
  if (firstRealizeLayer >= 0)
    details += " | first realizeFailed layer " +
               std::to_string(firstRealizeLayer) +
               ": the seed point does not map onto this layer";
  if (firstNearestLayer >= 0)
    details += " | first nearestFailed layer " +
               std::to_string(firstNearestLayer) +
               ": no inner-track point near the mapped seed";
  if (firstRotateLayer >= 0)
    details += " | first rotateFailed layer " +
               std::to_string(firstRotateLayer) + ": " + firstRotateReason;
  if (alignedCount == 0) {
    outStatus = "failed: 0/" + std::to_string(totalLayers) + counts + "." +
                details;
    zTsUtils::printDebug("[alignStartPoints] No paths could be aligned; the path topology, seam cut or transferred anchor is unsuitable.");
    return false;
  }

  outStatus = "aligned " + std::to_string(alignedCount) + "/" +
              std::to_string(totalLayers) + counts + ".";

  if (!details.empty())
    zTsUtils::printDebug("[alignStartPoints] Some paths could not be aligned: the seed or transferred seam location is unavailable, or the path is unsuitable for seam placement.");
  return true;
}

//--------------------------
//---- PUBLIC METHODS (pipeline order)
//--------------------------

ZSPACE_INLINE bool zPrintSyn::resamplePrintPaths(PrintPathArray &ioPaths2D,
                                      const PrintSynSetting &settings,
                                      std::string &outStatus,
                                      std::vector<int> *outBadLayers, bool useParallel) {
  if (outBadLayers)
    outBadLayers->clear();
  if (!std::isfinite(settings.FeatureAngleTol) || settings.FeatureAngleTol < 0.0f ||
      settings.FeatureAngleTol > 180.0f || !std::isfinite(settings.ResampleLength) ||
      settings.ResampleLength < zTsField::MinPhysicalParameter) {
    outStatus = "failed: invalid feature angle or resample length.";
    return false;
  }
  const int totalLayers = int(ioPaths2D.size());
  if (totalLayers == 0) {
    outStatus = "failed: No print paths available.";
    zTsUtils::printDebug("[resamplePrintPaths] " + outStatus);
    return false;
  }

  const std::vector<int> openBefore = openLeadPathLayers(ioPaths2D);

  const auto wallStart = std::chrono::high_resolution_clock::now();
  int totalPaths = 0;
  std::vector<char> resampleFailed(totalLayers, 0);

  zTsDetail::ParallelErrors parallelErrors3(totalLayers);
#pragma omp parallel for schedule(dynamic) if(useParallel)
  for (int l = 0; l < totalLayers; l++) {
    try {
      for (auto &path : ioPaths2D[l]) {
        if (reconstructCorners(path, settings.FeatureAngleTol) &&
            resampleSinglePath(path, settings.ResampleLength,
                               settings.FeatureAngleTol)) {
          filterGlitches(path);
        } else resampleFailed[l] = 1;
      }
      orientLayer(ioPaths2D[l], true);
#pragma omp atomic
      totalPaths += int(ioPaths2D[l].size());
    } catch (...) { parallelErrors3.capture(l); }
  }
  parallelErrors3.rethrow();

  const auto wallEnd = std::chrono::high_resolution_clock::now();
  const double wallTotalMs =
      std::chrono::duration<double, std::milli>(wallEnd - wallStart).count();

  std::ostringstream out;
  out << "Resampled, glitch-filtered and oriented " << totalPaths
      << " print paths in " << std::fixed << std::setprecision(2) << wallTotalMs
      << " ms (sampleLen=" << settings.ResampleLength << "mm)";
  outStatus = out.str();


  // Lead paths closed before this step but open after it were opened by the
  // resample / glitch filter / orientation, not by the isocontour.
  std::vector<int> opened;
  for (int l : openLeadPathLayers(ioPaths2D))
    if (!std::binary_search(openBefore.begin(), openBefore.end(), l))
      opened.push_back(l);
  if (outBadLayers)
    *outBadLayers = opened;
  if (!opened.empty())
    zTsUtils::printDebug("[resamplePrintPaths] A previously closed path became open during resampling; corner filtering or path reconstruction did not preserve closure.");
  std::vector<int> failed = opened;
  for (int l = 0; l < totalLayers; ++l)
    if (resampleFailed[l]) failed.push_back(l);
  std::sort(failed.begin(), failed.end());
  failed.erase(std::unique(failed.begin(), failed.end()), failed.end());
  if (outBadLayers) *outBadLayers = failed;
  if (std::find(resampleFailed.begin(), resampleFailed.end(), char(1)) != resampleFailed.end()) {
    outStatus = "failed: invalid or disconnected print curves on layers " + layerList(failed);
    zTsUtils::printDebug("[resamplePrintPaths] Invalid or disconnected print curves could not be resampled while preserving their connectivity.");
    return false;
  }
  return true;
}

ZSPACE_INLINE bool zPrintSyn::alignStartPoints(PrintPathArray &ioPaths2D,
                                    const zObjMeshArray &alignCharts2D,
                                    double seamOffset, double seamSize,
                                    std::string &outStatus,
                                    std::vector<int> *outBadLayers, bool useParallel) {
  std::vector<const zObjMesh *> charts;
  charts.reserve(alignCharts2D.size());
  for (const zObjMesh &mesh : alignCharts2D)
    charts.push_back(&mesh);
  return alignStartPointsOnCharts(ioPaths2D, charts, seamOffset, seamSize,
                                  outStatus, outBadLayers, useParallel);
}

ZSPACE_INLINE bool zPrintSyn::alignStartPoints(PrintPathArray &ioPaths2D,
                                    const zObjMesh &alignChart2D,
                                    double seamOffset, double seamSize,
                                    std::string &outStatus,
                                    std::vector<int> *outBadLayers, bool useParallel) {
  const std::vector<const zObjMesh *> charts(ioPaths2D.size(), &alignChart2D);
  return alignStartPointsOnCharts(ioPaths2D, charts, seamOffset, seamSize,
                                  outStatus, outBadLayers, useParallel);
}

ZSPACE_INLINE bool zPrintSyn::alignStartPoints(PrintPathArray &ioPaths2D,
                                    const zObjMeshArray &charts2D,
                                    const zObjMeshArray &sections3D,
                                    double seamOffset, double seamSize,
                                    std::string &outStatus,
                                    std::vector<int> *outBadLayers, bool useParallel) {
  if (outBadLayers)
    outBadLayers->clear();
  const int totalLayers = int(ioPaths2D.size());
  if (totalLayers == 0 || int(charts2D.size()) < totalLayers ||
      int(sections3D.size()) < totalLayers) {
    outStatus = "failed: charts/sections missing or size mismatch "
                "(totalLayers=" +
                std::to_string(totalLayers) + ", charts2D=" +
                std::to_string(int(charts2D.size())) + ", sections3D=" +
                std::to_string(int(sections3D.size())) + ").";
    zTsUtils::printDebug("[alignStartPoints] The path layers, unfolded charts and surface sections do not correspond.");
    return false;
  }

  // Seed: chart -> its own section (world); every layer: nearest point on
  // its section -> its chart.
  return alignStartPointsWith(
      ioPaths2D,
      [&, useParallel](int seedLayer, const zPoint &seed, zPointArray &outMapped,
          std::vector<char> &outValid) {
        const bool parallelMapping = useParallel;
        zPoint seedWorld;
        if (!mapPoint(seed, const_cast<zObjMesh &>(charts2D[seedLayer]),
                      const_cast<zObjMesh &>(sections3D[seedLayer]),
                      seedWorld))
          return false;
        outMapped.assign(totalLayers, zPoint());
        outValid.assign(totalLayers, 0);
        zTsDetail::ParallelErrors parallelErrors4(totalLayers);
#pragma omp parallel for schedule(dynamic) if(parallelMapping)
        for (int l = 0; l < totalLayers; l++) {
          try {
            zObjMesh &section = const_cast<zObjMesh &>(sections3D[l]);
            if (zFnMesh(section).numVertices() == 0)
              continue;
            zPoint mapped;
            if (mapPoint(seedWorld, section,
                         const_cast<zObjMesh &>(charts2D[l]), mapped)) {
              outMapped[l] = mapped;
              outValid[l] = 1;
            }
          } catch (...) { parallelErrors4.capture(l); }
        }
        parallelErrors4.rethrow();
        return true;
      },
      seamOffset, seamSize, outStatus, outBadLayers, useParallel);
}

namespace {
// Fast exact/near-roundoff chart hits avoid scanning every triangle for every
// dense source vertex. Unsupported hits retain remapGraph's full diagnostics
// and tolerance handling; realization remains the shared barycentric method.
bool cachedRemap(zObjGraph &graph,zObjMesh &section,const zTsDetail::SurfaceQuery &query,
    double tolerance,zVectorArray *normals,const zFloatArray *heights,zFloatArray *mappedHeights,
    const zVectorArray *vertexNormals){
  zPointArray p;zFnGraph fn(graph);fn.getVertexPositions(p);ParameterGraph params;fn.getEdgeData(params.edgeConnects);
  params.vertices.reserve(p.size());
  for(const auto &point:p){if(!std::isfinite(point.x)||!std::isfinite(point.y)||!std::isfinite(point.z))return false;
    zTsDetail::SurfaceQuery::Hit hit;if(!query.uniqueChartHit(point,hit,tolerance))return false;
    const auto &ids=query.triangle(hit.Triangle).Ids;
    params.vertices.push_back({ids[0],ids[1],ids[2],float(hit.Weights[0]),float(hit.Weights[1]),float(hit.Weights[2])});}
  return realizeParameterGraph(params,section,graph,normals,heights,mappedHeights,vertexNormals);
}
}

ZSPACE_INLINE bool zPrintSyn::projectPrintPaths(const PrintPathArray &paths2D,
                                     const zObjMeshArray &charts2D,
                                     const zObjMeshArray &sections3D,
                                     const zFloatArray *vertexLayerHeights,
                                     PrintPathArray &outPaths3D,
                                     std::string &outStatus,
                                     std::vector<int> *outBadLayers, bool useParallel,
                                     const std::vector<zFloatArray> *perLayerVertexHeights,
                                     const std::vector<zVectorArray> *perLayerVertexNormals) {
  outPaths3D.clear();
  if (outBadLayers)
    outBadLayers->clear();
  const int totalLayers = int(paths2D.size());
  if (totalLayers == 0 || charts2D.empty() || sections3D.empty()) {
    outStatus =
        "failed: No 2D print paths or section meshes available for projection.";
    zTsUtils::printDebug("[projectPrintPaths] " + outStatus);
    return false;
  }
  static const zFloatArray emptyVertexLayerHeights;
  if (perLayerVertexNormals) {
    if (perLayerVertexNormals->size() != sections3D.size() || sections3D.size() < size_t(totalLayers)) {
      outStatus = "failed: normal tables and sections must have one entry per layer.";
      return false;
    }
    for (int l=0;l<totalLayers;++l) {
      const auto &normals=(*perLayerVertexNormals)[l];
      if ((!normals.empty() && normals.size()!=size_t(zFnMesh(const_cast<zObjMesh &>(sections3D[l])).numVertices())) ||
          std::any_of(normals.begin(),normals.end(),[](const zVector &n){return !std::isfinite(n.x)||!std::isfinite(n.y)||!std::isfinite(n.z)||double(n.x)*n.x+double(n.y)*n.y+double(n.z)*n.z<=1e-12;})) {
        outStatus="failed: invalid per-vertex normal table on layer "+std::to_string(l)+".";
        if(outBadLayers)outBadLayers->push_back(l);
        return false;
      }
    }
  }
  if (perLayerVertexHeights) {
    if (perLayerVertexHeights->size() != sections3D.size() || sections3D.size() < size_t(totalLayers)) {
      outStatus = "failed: height tables and sections must have one entry per layer.";
      return false;
    }
    for (int l = 0; l < totalLayers; ++l) {
      const auto &heights = (*perLayerVertexHeights)[l];
      if (heights.size() != size_t(zFnMesh(const_cast<zObjMesh &>(sections3D[l])).numVertices()) ||
          std::any_of(heights.begin(), heights.end(), [](float h) { return !std::isfinite(h) || h < 0.0f; })) {
        outStatus = "failed: invalid per-vertex height table on layer " + std::to_string(l) + ".";
        if (outBadLayers) outBadLayers->push_back(l);
        return false;
      }
    }
  }

  outPaths3D.assign(totalLayers, std::vector<PrintPath>());
  for (int l = 0; l < totalLayers; l++) {
    outPaths3D[l].assign(paths2D[l].size(), PrintPath());
  }

  int totalPaths = 0;
  std::vector<char> remapFailed(totalLayers, 0);
  std::vector<std::string> remapErrors(totalLayers);
  const auto wallStart = std::chrono::high_resolution_clock::now();

  zTsDetail::ParallelErrors parallelErrors5(totalLayers);
#pragma omp parallel for schedule(dynamic) if(useParallel)
  for (int l = 0; l < totalLayers; l++) {
    try {
      const zFloatArray &meshVertexHeights = perLayerVertexHeights
          ? (*perLayerVertexHeights)[l]
          : (vertexLayerHeights ? *vertexLayerHeights : emptyVertexLayerHeights);
      std::unique_ptr<zTsDetail::SurfaceQuery> chartQuery;double chartTolerance=0;
      if(l<int(charts2D.size())&&l<int(sections3D.size())){
        auto &chart=const_cast<zObjMesh&>(charts2D[l]);auto &section=const_cast<zObjMesh&>(sections3D[l]);
        zIntArray a,b,c,d;zFnMesh(chart).getPolygonData(a,b);zFnMesh(section).getPolygonData(c,d);
        zPointArray p;zFnMesh(chart).getVertexPositions(p);double scale=1;bool finite=true;
        for(const auto &v:p){finite&=std::isfinite(v.x)&&std::isfinite(v.y)&&std::isfinite(v.z);scale=(std::max)(scale,zTsUtils::toEigenVec(v).cwiseAbs().maxCoeff());}
        if(finite&&p.size()==size_t(zFnMesh(section).numVertices())&&a==c&&b==d){chartQuery=std::make_unique<zTsDetail::SurfaceQuery>(chart);chartTolerance=(std::max)(1e-7,8*std::numeric_limits<float>::epsilon()*scale);}
      }
      for (size_t r = 0; r < paths2D[l].size(); r++) {
        outPaths3D[l][r].graph = paths2D[l][r].graph;
        outPaths3D[l][r].widths = paths2D[l][r].widths;
        outPaths3D[l][r].trackIds = paths2D[l][r].trackIds;
        outPaths3D[l][r].featureFlags = paths2D[l][r].featureFlags;

        zFnGraph fnContour(outPaths3D[l][r].graph);
        if (fnContour.numVertices() > 0 && l < int(charts2D.size()) &&
            l < int(sections3D.size())) {
          std::string reason;
          const bool validFeatures = paths2D[l][r].featureFlags.empty() ||
              paths2D[l][r].featureFlags.size() == size_t(fnContour.numVertices());
          if (!validFeatures) reason = "FeatureCountMismatch";
          const bool cached=validFeatures&&chartQuery&&cachedRemap(outPaths3D[l][r].graph,
              const_cast<zObjMesh&>(sections3D[l]),*chartQuery,chartTolerance,&outPaths3D[l][r].normals,
              &meshVertexHeights,&outPaths3D[l][r].heights,perLayerVertexNormals?&(*perLayerVertexNormals)[l]:nullptr);
          if (!validFeatures || (!cached&&!remapGraph(outPaths3D[l][r].graph,
                          const_cast<zObjMesh &>(charts2D[l]),
                          const_cast<zObjMesh &>(sections3D[l]),
                          &outPaths3D[l][r].normals, &meshVertexHeights,
                          &outPaths3D[l][r].heights, &reason,
                          perLayerVertexNormals ? &(*perLayerVertexNormals)[l] : nullptr))) {
            remapFailed[l] = 1;
            if (remapErrors[l].empty()) remapErrors[l] = "layer " + std::to_string(l) +
                ", path " + std::to_string(r) + ": " + reason;
            outPaths3D[l][r] = PrintPath();
          }
        } else if (fnContour.numVertices() > 0) {
          remapFailed[l] = 1;
          if (remapErrors[l].empty()) remapErrors[l] = "layer " + std::to_string(l) +
              ", path " + std::to_string(r) + ": MissingChartOrSection";
          outPaths3D[l][r] = PrintPath();
        }
      }
#pragma omp atomic
      totalPaths += int(paths2D[l].size());
    } catch (...) { parallelErrors5.capture(l); }
  }
  parallelErrors5.rethrow();

  const auto wallEnd = std::chrono::high_resolution_clock::now();
  const double wallTotalMs =
      std::chrono::duration<double, std::milli>(wallEnd - wallStart).count();

  std::ostringstream out;
  out << "Projected " << totalPaths << " 3D PrintPaths in " << std::fixed
      << std::setprecision(2) << wallTotalMs << " ms";
  outStatus = out.str();


  std::vector<int> failed;
  for (int l = 0; l < totalLayers; l++)
    if (remapFailed[l])
      failed.push_back(l);
  if (outBadLayers)
    *outBadLayers = failed;
  if (!failed.empty()) {
    outStatus = "failed: 2D-to-3D mapping on " + remapErrors[failed.front()];
    zTsUtils::printDebug("[projectPrintPaths] Some paths could not be mapped to their matching surface charts; their invalid projected graphs were discarded.");
  }
  return failed.empty();
}

ZSPACE_INLINE bool zPrintSyn::linkPrintPaths(const PrintPathArray &paths3D,
                                  PrintPath &outCombined,
                                  std::string &outStatus) {
  outCombined = PrintPath();

  bool hasHeights = false, hasNormals = false, hasTrackIds = false, hasFeatures = false;
  bool hasMultiplePaths = false;
  for (size_t layer = 0; layer < paths3D.size(); ++layer) {
    const auto &ringPaths = paths3D[layer];
    if (ringPaths.size() > 1) {
      hasMultiplePaths = true;
      zTsUtils::printDebug("[linkPrintPaths] layer=" + std::to_string(layer) +
          " MultiplePrintPaths: pathCount=" + std::to_string(ringPaths.size()) +
          "; paths remain disconnected in the combined graph.");
    }
    for (const PrintPath &path : ringPaths) {
      if (!path.heights.empty())
        hasHeights = true;
      if (!path.normals.empty())
        hasNormals = true;
      if (!path.trackIds.empty())
        hasTrackIds = true;
      if (!path.featureFlags.empty()) {
        if (path.featureFlags.size() != size_t(zFnGraph(const_cast<zObjGraph &>(path.graph)).numVertices())) {
          outStatus = "failed: feature count does not match path vertices on layer " + std::to_string(layer);
          return false;
        }
        hasFeatures = true;
      }
    }
  }

  zPointArray combinedPositions;
  zFloatArray combinedWidths;
  zFloatArray combinedHeights;
  zVectorArray combinedNormals;
  zFloatArray combinedTrackIds;
  zIntArray combinedEdges;
  zIntArray combinedFeatures;

  for (const auto &ringPaths : paths3D) {
    for (const PrintPath &layerPath : ringPaths) {

      zFnGraph fnGraph(const_cast<zObjGraph &>(layerPath.graph));
      zPointArray positions;
      fnGraph.getVertexPositions(positions);
      if (positions.empty())
        continue;
      const int base = static_cast<int>(combinedPositions.size());
      if (hasMultiplePaths) {
        zIntArray edges;
        fnGraph.getEdgeData(edges);
        for (int id : edges) combinedEdges.push_back(base + id);
      }

      combinedPositions.insert(combinedPositions.end(), positions.begin(),
                               positions.end());
      combinedWidths.insert(combinedWidths.end(), layerPath.widths.begin(),
                            layerPath.widths.end());
      if (hasFeatures) {
        if (!layerPath.featureFlags.empty())
          combinedFeatures.insert(combinedFeatures.end(), layerPath.featureFlags.begin(), layerPath.featureFlags.end());
        else combinedFeatures.insert(combinedFeatures.end(), positions.size(), 0);
      }
      if (hasHeights) {
        if (!layerPath.heights.empty())
          combinedHeights.insert(combinedHeights.end(), layerPath.heights.begin(),
                                 layerPath.heights.end());
        else
          combinedHeights.insert(combinedHeights.end(), positions.size(), 0.0f);
      }
      if (hasNormals) {
        if (!layerPath.normals.empty())
          combinedNormals.insert(combinedNormals.end(), layerPath.normals.begin(),
                                 layerPath.normals.end());
        else
          combinedNormals.insert(combinedNormals.end(), positions.size(),
                                 zVector(0, 0, 1));
      }
      if (hasTrackIds) {
        if (!layerPath.trackIds.empty())
          combinedTrackIds.insert(combinedTrackIds.end(),
                                  layerPath.trackIds.begin(),
                                  layerPath.trackIds.end());
        else
          combinedTrackIds.insert(combinedTrackIds.end(), positions.size(), 0.0f);
      }
    }
  }

  const int total = static_cast<int>(combinedPositions.size());
  if (total < 2) {
    outStatus = "failed: fewer than 2 print-path vertices to link.";
    zTsUtils::printDebug("[linkPrintPaths] Not enough usable path vertices remain to form a connection.");
    return false;
  }

  // Preserve disjoint curves when there is no explicit travel path between
  // components. Keep the existing inter-layer link for single-path layers.
  if (!hasMultiplePaths) {
    combinedEdges.reserve((total - 1) * 2);
    for (int i = 0; i + 1 < total; i++) {
      combinedEdges.push_back(i);
      combinedEdges.push_back(i + 1);
    }
  }

  zFnGraph fnOut(outCombined.graph);
  fnOut.create(combinedPositions, combinedEdges);
  outCombined.widths = combinedWidths;
  outCombined.featureFlags = combinedFeatures;
  if (hasHeights)
    outCombined.heights = combinedHeights;
  if (hasNormals)
    outCombined.normals = combinedNormals;
  if (hasTrackIds)
    outCombined.trackIds = combinedTrackIds;

  outStatus = "Combined " + std::to_string(total) + " print-path vertices" +
      (hasMultiplePaths ? "; disconnected curves retained without travel links." : " into one print path.");

  return true;
}

ZSPACE_INLINE bool zPrintSyn::computePipeline(const zObjMeshArray &charts2D,
    const zObjMeshArray &sections3D, const zFloatArray *vertexLayerHeights,
    SeamMapping seamMapping, const PrintSynSetting &settings,
    PrintPathArray &outPaths2D, PrintPathArray &outPaths3D,
    PrintPath &outCombined, std::string &outStatus, bool useParallel,
    const std::vector<zFloatArray> *perLayerVertexHeights,
    const std::vector<zVectorArray> *perLayerVertexNormals) {
  return computePipeline(charts2D,sections3D,vertexLayerHeights,seamMapping,settings,
      outPaths2D,outPaths3D,outCombined,outStatus,useParallel,perLayerVertexHeights,perLayerVertexNormals,nullptr);
}

ZSPACE_INLINE bool zPrintSyn::computePipeline(const zObjMeshArray &charts2D,
    const zObjMeshArray &sections3D, const zFloatArray *vertexLayerHeights,
    SeamMapping seamMapping, const PrintSynSetting &settings,
    PrintPathArray &outPaths2D, PrintPathArray &outPaths3D,
    PrintPath &outCombined, std::string &outStatus, bool useParallel,
    const std::vector<zFloatArray> *perLayerVertexHeights,
    const std::vector<zVectorArray> *perLayerVertexNormals,PrintPathArray *outSurfaceContours2D) {
  if(outSurfaceContours2D)outSurfaceContours2D->clear();
  outPaths3D.clear();
  outCombined = PrintPath();
  const size_t n = outPaths2D.size();
  if (n == 0 || charts2D.size() < n || sections3D.size() < n) {
    outStatus = "failed: paths and charts must have one entry per layer.";
    return false;
  }
  std::vector<int> badLayers;
  const auto stopAt = [&](bool ok, const std::string &step) {
    if (ok && badLayers.empty()) return false;
    if(outSurfaceContours2D)outSurfaceContours2D->clear();
    std::string layers;
    for (size_t i = 0; i < badLayers.size(); ++i)
      layers += (i ? "," : "") + std::to_string(badLayers[i]);
    const std::string reason = ok ? "problems on layer(s) " + layers
        : outStatus + (layers.empty() ? "" : " (layers " + layers + ")");
    if (step != "projectPrintPaths" && step != "resamplePrintPaths3D" && step != "linkPrintPaths") {
      std::string projectStatus;
      if (!projectPrintPaths(outPaths2D, charts2D, sections3D, vertexLayerHeights,
                             outPaths3D, projectStatus, nullptr, useParallel, perLayerVertexHeights, perLayerVertexNormals)) outPaths3D.clear();
    }
    outCombined = PrintPath();
    outStatus = "stopped at " + step + ": " + reason + " - the paths are kept as this step left them.";
    return true;
  };
  if (stopAt(resamplePrintPaths(outPaths2D, settings, outStatus, &badLayers, useParallel),
              "resamplePrintPaths")) return false;
  const bool aligned = seamMapping == SeamMapping::SharedTopology
      ? alignStartPoints(outPaths2D, charts2D, settings.SeamOffset, settings.SeamSize, outStatus, &badLayers, useParallel)
      : alignStartPoints(outPaths2D, charts2D, sections3D, settings.SeamOffset, settings.SeamSize, outStatus, &badLayers, useParallel);
  if (stopAt(aligned, "alignStartPoints")) return false;
  if (stopAt(projectPrintPaths(outPaths2D, charts2D, sections3D, vertexLayerHeights,
                              outPaths3D, outStatus, &badLayers, useParallel, perLayerVertexHeights, perLayerVertexNormals), "projectPrintPaths")) return false;
  if(outSurfaceContours2D)*outSurfaceContours2D=outPaths2D;
  badLayers.clear();
  if (stopAt(resamplePrintPaths3D(outPaths2D,outPaths3D,charts2D,sections3D,settings,
      outStatus,useParallel,vertexLayerHeights,perLayerVertexHeights,perLayerVertexNormals),"resamplePrintPaths3D")) return false;
  if (stopAt(linkPrintPaths(outPaths3D, outCombined, outStatus), "linkPrintPaths")) return false;
  outStatus = "computePostProcess: " + std::to_string(n) + " layers processed. " + outStatus;
  return true;
}

} // namespace zSpace
