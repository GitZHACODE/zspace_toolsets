#pragma once
// Frozen pre-optimization scalar classifier: independent parity oracle.
#include <zspace/slicing.h>
#include <algorithm>
#include <cmath>
#include <limits>
namespace face_tag_reference {
using namespace zSpace;
using namespace zSpace::zTsUtils;
using namespace zSpace::zTsField;
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


namespace {
template<class Work> bool withFaceTagQuery(const FaceTagSource &source, Work &&work) {
  zPointArray boundary;
  zIntArray edges;
  zFnGraph(const_cast<zObjGraph &>(source.boundary)).getVertexPositions(boundary);
  zFnGraph(const_cast<zObjGraph &>(source.boundary)).getEdgeData(edges);
  if (source.edgeTags.size() != edges.size()/2 || source.edgeWidthIds.size() != edges.size()/2)
    return false;
  std::vector<FastSegment> boundaryEdges, slotEdges;
  computeFastSegments(source.boundary, boundaryEdges);
  computeFastSegments(source.slot, slotEdges);
  std::vector<FieldSegment> groups;
  loopWidthSegments(source.boundary, source.edgeWidthIds, source.widths, groups);
  std::vector<std::vector<FastSegment>> groupEdges, ribs, ribSlots;
  for (const auto &g : groups) { groupEdges.emplace_back(); computeFastSegments(g.graph, groupEdges.back()); }
  for (const auto &g : source.bracing) {
    ribs.emplace_back(); computeFastSegments(g.graph, ribs.back());
    zObjGraph slot; ribSlots.emplace_back();
    if (g.trim >= 0 && g.trim <= 1 && createSlot(g.graph, slot, g.trim, g.width*1.5))
      computeFastSegments(slot, ribSlots.back());
  }
  zIntArray groupIds = source.edgeWidthIds;
  std::sort(groupIds.begin(), groupIds.end());
  groupIds.erase(std::unique(groupIds.begin(), groupIds.end()), groupIds.end());
  const auto distance = [&](const zPoint &p, const std::vector<FastSegment> &segments) {
    zScalarArray d; computeDF(&p, 1, segments, d); return double(d[0]);
  };
  const auto label = [&](const zPoint &p) {
    zScalarArray sign; computeSDF(&p, 1, boundaryEdges, sign);
    double outer = -1e30, inner = -1e30; int group = -1;
    for (size_t g = 0; g < groups.size(); ++g) {
      const double d = distance(p, groupEdges[g]) * (sign[0] <= 0 ? -1 : 1);
      const double half = .5 * printedWidth(groups[g].width, true, source.doubleBoundary, source.overlap);
      if (d + half > outer) { outer = d + half; group = int(g); }
      inner = (std::max)(inner, d + groups[g].width - half);
    }
    if (group < 0) return 0;
    if (source.doubleBoundary) {
      double bracing = 1e30;
      for (size_t r = 0; r < ribs.size(); ++r) {
        const auto &g = source.bracing[r];
        const double half = .5 * printedWidth(g.width, false, true, source.overlap);
        double value = distance(p, ribs[r]) - (.5*g.width-half);
        if (!ribSlots[r].empty()) value = (std::max)(value, -(distance(p, ribSlots[r])-half));
        bracing = (std::min)(bracing, value);
      }
      inner = (std::max)(inner, -bracing);
      if (-inner > outer + 1e-8) return 0;
      if (!slotEdges.empty() && source.slotWidth > 0 &&
          -(distance(p, slotEdges) - .5*printedWidth(source.slotWidth, true, true, source.overlap)) > outer + 1e-8)
        return 0;
    }
    // Resolve the owner inside the winning width field, keeping widths and
    // side identity independent. Deterministic ties use the smaller tag.
    double best = 1e30; int tag = 0;
    for (size_t e = 0; e < boundaryEdges.size(); ++e) {
      if (source.edgeWidthIds[e] != groupIds[group]) continue;
      const double d = distance(p, std::vector<FastSegment>{boundaryEdges[e]});
      const int t = source.edgeTags[e];
      if (d < best-1e-8 || (std::abs(d-best) <= 1e-8 && (!tag || t < tag))) { best=d; tag=t; }
    }
    return isSideFaceTag(tag) ? tag : 0;
  };
  return work(label);
}
}

inline bool queryFaceTags(const FaceTagSource &source, const zPointArray &points, zIntArray &tags) {
  tags.clear();
  return withFaceTagQuery(source,[&](const auto &label) {
    tags.reserve(points.size());for(const auto &point:points)tags.push_back(label(point));return true;
  });
}

inline bool labelFaceTagSegments(PrintPath &path, const FaceTagSource &source, bool splitTransitions) {
  return withFaceTagQuery(source,[&](const auto &label) {
  std::vector<GraphPolyline> curves;
  if (!graphPolylines(path.graph, curves) || curves.size()!=1) return false;
  auto ids = curves[0].Vertices;
  if (ids.size()>1 && ids.front()==ids.back()) ids.pop_back();
  const bool closed=curves[0].Closed;
  zPointArray points; zFnGraph(path.graph).getVertexPositions(points);
  if (points.size()!=path.widths.size() || ids.size()<2) return false;
  double coordinateScale=1;
  for(const auto &p:points)coordinateScale=(std::max)(coordinateScale,(std::max)(std::abs(double(p.x)),std::abs(double(p.y))));
  const double anchorTolerance=4*std::numeric_limits<float>::epsilon()*coordinateScale;
  PrintPath result; zPointArray output; zIntArray connects;
  const auto append = [&](int a, int b, double t) {
    // Intentional new contract: original endpoints are copied exactly. The
    // independent historical classifier and interior transition math stay intact.
    const auto p=t==0?points[a]:t==1?points[b]:points[a]+(points[b]-points[a])*t;
    const auto value=[&](float a,float b){return t==0?a:t==1?b:float(a+(b-a)*t);};
    if (!output.empty() && vecDistSq(output.back(),p)<1e-18) return;
    output.push_back(p);
    result.featureFlags.push_back(path.featureFlags.empty()?0:t<=1e-8?path.featureFlags[a]:t>=1-1e-8?path.featureFlags[b]:0);
    result.widths.push_back(value(path.widths[a],path.widths[b]));
    if (!path.heights.empty()) result.heights.push_back(value(path.heights[a],path.heights[b]));
    if (!path.normals.empty()) { auto n=t==0?path.normals[a]:t==1?path.normals[b]:path.normals[a]+(path.normals[b]-path.normals[a])*t; if(n.length()>0) n.normalize(); result.normals.push_back(n); }
    if (!path.trackIds.empty()) result.trackIds.push_back((path.trackIds[a]+(path.trackIds[b]-path.trackIds[a])*t)>=.5 ? 1.f:0.f);
  };
  append(ids[0],ids[0],0);
  for (size_t e=0;e<ids.size()-(closed?0:1);++e) {
    const int a=ids[e], b=ids[(e+1)%ids.size()];
    const auto at=[&](double t){ return label(points[a]+(points[b]-points[a])*t); };
    // Raw contour edges are grid-scale; subdivision also catches two
    // transitions on one edge near a slot/corner. IDs are never interpolated.
    double previous=1e-7; int last=at(previous);
    if(splitTransitions) for(int q=1;q<=8;++q) {
      const double t=q==8 ? 1-1e-7 : double(q)/8;
      const int current=at(t);
      if(current!=last) {
        double lo=previous,hi=t;
        for(int it=0;it<30;++it) { const double mid=.5*(lo+hi); if(at(mid)==last) lo=mid; else hi=mid; }
        double cut=.5*(lo+hi);
        const double length=std::sqrt(vecDistSq(points[a],points[b]));
        // Reuse a geometric feature when a categorical crossing differs only
        // by float precision. Distinct recorded geometric features stay intact.
        if(path.featureFlags.size()==points.size()) {
          if(path.featureFlags[a] && cut*length<=anchorTolerance)cut=0;
          else if(path.featureFlags[b] && (1-cut)*length<=anchorTolerance)cut=1;
        }
        append(a,b,cut);
      }
      previous=t; last=current;
    }
    if (!closed || e+1<ids.size()) append(a,b,1);
  }
  result.faceTags.assign(output.size(),0);
  for(size_t i=0;i<output.size()-(closed?0:1);++i) {
    const int a=int(i),b=int((i+1)%output.size());
    connects.push_back(a);connects.push_back(b);
    const int tag=label((output[a]+output[b])*.5);
    result.edgeFaceTags.push_back(tag);
    for(int v:{a,b}) if(tag && (!result.faceTags[v] || tag<result.faceTags[v])) result.faceTags[v]=tag;
  }
  zFnGraph(result.graph).create(output,connects); path=std::move(result); return true;
  });
}


} // namespace face_tag_reference
