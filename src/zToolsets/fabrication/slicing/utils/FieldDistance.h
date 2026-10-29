#pragma once

#include <zspace/zToolsets/fabrication/slicing/utils/Utils.h>

namespace zSpace { namespace zTsDetail {
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

// Shared scalar primitive for width/path fields and categorical ownership.
// Preserve the original float arithmetic, including the short-edge threshold.
inline float segmentDistanceSquared(float px, float py, const FastSegment &seg) {
  float t = ((px - seg.x0) * seg.dx + (py - seg.y0) * seg.dy) * seg.lenSqInv;
  t = (std::max)(0.0f, (std::min)(1.0f, t));
  const float projX = seg.x0 + t * seg.dx;
  const float projY = seg.y0 + t * seg.dy;
  return (px - projX) * (px - projX) + (py - projY) * (py - projY);
}


} } // namespace zSpace::zTsDetail
