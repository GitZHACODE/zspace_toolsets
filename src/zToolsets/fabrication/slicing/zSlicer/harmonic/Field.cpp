#ifndef EIGEN_DONT_PARALLELIZE
#define EIGEN_DONT_PARALLELIZE
#endif

#include <zspace/zToolsets/fabrication/slicing/utils/Parallel.h>
#include <zspace/zToolsets/fabrication/slicing/zSlicer/HarmonicSlicer.h>
#include <zspace/zToolsets/fabrication/slicing/utils/Utils.h>

#ifndef ZSPACE_INLINE
#define ZSPACE_INLINE ZSPACE_TOOLSETS_INLINE
#endif

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <unordered_map>
#include <utility>
#include <vector>

namespace zSpace {

using namespace zTsUtils;

// The SDF field kernels (segment distances, width/path fields, trim slots)
// are shared with TopologySlicer/zPrintSyn - see Utils.h. Here they run once
// per layer on its o_IsoMeshes2D chart (computeContours()), exactly as TopologySlicer
// runs them on its unrolled layers.

namespace {

// Splits every edge of a graph into pieces no longer than step (the original
// vertices are kept), so a curve mapped point-by-point through a chart
// follows the chart's bending instead of only its end points.
void subdivideGraphEdges(const zPointArray &points, const zIntArray &edges, double step, zPointArray &outPoints,
                         zIntArray &outEdges) {
  outPoints = points;
  outEdges.clear();
  for (size_t e = 0; e + 1 < edges.size(); e += 2) {
    const int a = edges[e], b = edges[e + 1];
    const zPoint &pa = points[a], &pb = points[b];
    const double length = std::sqrt(vecDistSq(pa, pb));
    const int pieces = (step > 1e-9) ? (std::max)(1, static_cast<int>(std::ceil(length / step))) : 1;
    int previous = a;
    for (int i = 1; i < pieces; ++i) {
      const double t = static_cast<double>(i) / pieces;
      outPoints.emplace_back(pa.x + t * (pb.x - pa.x), pa.y + t * (pb.y - pa.y), pa.z + t * (pb.z - pa.z));
      const int current = static_cast<int>(outPoints.size()) - 1;
      outEdges.push_back(previous);
      outEdges.push_back(current);
      previous = current;
    }
    outEdges.push_back(previous);
    outEdges.push_back(b);
  }
}

} // namespace

//--------------------------
//---- PUBLIC METHODS: STEP 7 (PER-LAYER PATH/WIDTH FIELD)
//--------------------------

// Builds every layer's path/width/track SDF fields in its o_IsoMeshes2D chart, on
// one grid shared by all layers. (u,v) is used only to put the bracing and
// the boundary slot anchor onto each layer: their floor (u,v) is located on
// the layer's o_IsoMeshUV and realized with the same barycentric coordinates
// on o_IsoMeshes2D.
ZSPACE_INLINE bool HarmonicSlicer::computeField(bool useParallel) {
  const bool success = computeFieldImpl(useParallel);
  if (success) zTsUtils::printSuccess("[HarmonicSlicer::computeField] Success");
  else zTsUtils::printError("[HarmonicSlicer::computeField] Failed");
  return success;
}

ZSPACE_INLINE bool HarmonicSlicer::computeFieldImpl(bool useParallel) {
  if (!zTsField::validateFieldSettings(i_Settings)) {
    _Status = "failed: invalid physical field parameters (minimum 1e-4 mm).";
    return false;
  }
  clearFields();

  if (!i_SlotPoint) {
    _Status = "failed: i_SlotPoint is required. Call setSlotAnchor() before computing.";
    printDebug("[computeFieldImpl] " + _Status);
    return false;
  }

  const int totalLayers = static_cast<int>(o_IsoMeshes2D.size());
  if (totalLayers == 0 || o_IsoMeshUV.size() != o_IsoMeshes2D.size() ||
      o_IsoMeshes.size() != o_IsoMeshes2D.size() || o_IsoMeshBoundary2D.size() != o_IsoMeshes2D.size()) {
    _Status = "failed: no unrolled layers. Call computeContours() first.";
    printDebug("[computeFieldImpl] " + _Status);
    return false;
  }
  if (i_Settings.BoundaryWidth.empty()) {
    _Status = "failed: BoundaryWidth is empty - one entry per width group is needed.";
    printDebug("[computeFieldImpl] " + _Status);
    return false;
  }
  const double cellSize = static_cast<double>(i_Settings.FieldCellSize);
  const double scaleUp = (i_Settings.FieldScaleUp > 1e-6f) ? static_cast<double>(i_Settings.FieldScaleUp) : 1.0;

  // The bottom cap and all its actual side-boundary components.
  zObjMesh &bottomMesh = o_IsoMeshes[0];

  // Bottom (u,v) of a world-coordinate graph: every vertex located on the
  // bottom layer (nearest point), its (u,v) interpolated from o_IsoMeshUV[0].
  const auto floorUVGraph = [&](const zObjGraph &worldGraph, zObjGraph &outUVGraph) -> bool {
    ParameterGraph params;
    zPointArray uvPts;
    if (!parameterizeGraph(worldGraph, bottomMesh, params) ||
        !interpolateParameterGraph(params, o_IsoMeshUV[0], uvPts))
      return false;
    for (zPoint &p : uvPts) p.z = 0.0;
    zIntArray edges = params.edgeConnects;
    zFnGraph(outUVGraph).create(uvPts, edges);
    return true;
  };

  // Bracing (WORLD coordinates, see i_InputBracing), its edges first
  // subdivided (2 cells) so each layer's (u,v) placement can bend them.
  zObjGraph bracingUV;
  bool hasBracing = false;
  if (i_InputBracing) {
    zFnGraph fnSource(*i_InputBracing);
    zPointArray worldPts;
    zIntArray edgeConnects;
    fnSource.getVertexPositions(worldPts);
    fnSource.getEdgeData(edgeConnects);
    if (!worldPts.empty() && !edgeConnects.empty()) {
      zPointArray pts;
      zIntArray edges;
      subdivideGraphEdges(worldPts, edgeConnects, 2.0 * cellSize, pts, edges);
      zObjGraph bracingWorld;
      zFnGraph(bracingWorld).create(pts, edges);
      hasBracing = floorUVGraph(bracingWorld, bracingUV);
      if (!hasBracing)
        printDebug("[computeFieldImpl] layer=0 could not locate the bracing on the floor - no bracing.");
    }
  }

  // Locate the required i_SlotPoint on the world bottom loop and map to UV.
  zObjGraph slotUV;
  bool hasSlotAnchor = false;
  std::string slotAnchorFailure = "NoUsableBottomLoop";
  if (i_Settings.DoubleBoundary) {
    zTsField::LoopAnchor floorAnchor;
    zPointArray bottomLoop;
    zIntArray bottomEdges;
    if (!o_IsoMeshBoundary.empty() &&
        zTsField::locateGraphLoopAnchor(o_IsoMeshBoundary[0], *i_SlotPoint,
                                       bottomLoop, bottomEdges, floorAnchor)) {
      zPointArray pts = {floorAnchor.point, floorAnchor.point}; // a point as a zero-length edge
      zIntArray edges = {0, 1};
      zObjGraph anchorWorld;
      zFnGraph(anchorWorld).create(pts, edges);
      hasSlotAnchor = floorUVGraph(anchorWorld, slotUV);
      if (!hasSlotAnchor) slotAnchorFailure = "BottomAnchorUVMappingFailed";
      if (!hasSlotAnchor)
        printDebug("[computeFieldImpl] layer=0 could not locate the slot anchor on the floor - no boundary slot.");
    }
  }

  // Per layer, as in TopologySlicer: the bottom (u,v) located on the
  // layer's o_IsoMeshUV (triangle + barycentric coordinates; triangles with
  // an invalid (u,v) are skipped) and those same coordinates realized on
  // o_IsoMeshes2D.
  o_BracingGraphs2D.assign(totalLayers, zObjGraph());
  zPointArray slotPoints2D(totalLayers);
  std::vector<char> slotValid(totalLayers, 0);
  zTsDetail::ParallelErrors parallelErrors1(totalLayers);
#pragma omp parallel for schedule(dynamic) if(useParallel)
  for (int k = 0; k < totalLayers; ++k) {
    try {
      zObjMesh &unrolled = o_IsoMeshes2D[k];
      zFnMesh fnUnroll(unrolled);
      if (fnUnroll.numVertices() == 0 || static_cast<int>(o_IsoMeshUV[k].size()) != fnUnroll.numVertices()) continue;
      zPointArray layerUV = o_IsoMeshUV[k];
      for (zPoint &p : layerUV) p.z = 0.0;
      if (hasBracing) {
        ParameterGraph params;
        // Transported UV is not guaranteed bijective or identical across layers.
        // Placement may use the nearest valid UV triangle; final path projection
        // still uses strict containment in the actual 2D chart.
        if (parameterizeGraph(bracingUV, unrolled, params, &layerUV, true)) {
          realizeParameterGraph(params, unrolled, o_BracingGraphs2D[k]);
        }
      }
      if (hasSlotAnchor) {
        ParameterGraph params;
        zPointArray unrollPositions, located;
        fnUnroll.getVertexPositions(unrollPositions);
        if (parameterizeGraph(slotUV, unrolled, params, &layerUV, true) &&
            interpolateParameterGraph(params, unrollPositions, located) && !located.empty()) {
          slotPoints2D[k] = located[0];
          slotValid[k] = 1;
        }
      }
    } catch (...) { parallelErrors1.capture(k); }
  }
  parallelErrors1.rethrow();

  // One grid over every layer's o_IsoMeshes2D chart and bracing.
  zDomain<zPoint> gridBounds;
  int resX = 0, resY = 0;
  {
    zPoint minBB((std::numeric_limits<double>::max)(), (std::numeric_limits<double>::max)(), 0.0);
    zPoint maxBB(-(std::numeric_limits<double>::max)(), -(std::numeric_limits<double>::max)(), 0.0);
    const auto accumulate = [&](const zPoint &p) {
      minBB.x = (std::min)(minBB.x, p.x);
      minBB.y = (std::min)(minBB.y, p.y);
      maxBB.x = (std::max)(maxBB.x, p.x);
      maxBB.y = (std::max)(maxBB.y, p.y);
    };
    for (int k = 0; k < totalLayers; ++k) {
      zFnMesh fn(o_IsoMeshes2D[k]);
      if (fn.numVertices() > 0) {
        zPoint meshMin, meshMax;
        fn.getBounds(meshMin, meshMax);
        accumulate(meshMin);
        accumulate(meshMax);
      }
      zPointArray bracingPts;
      zFnGraph(o_BracingGraphs2D[k]).getVertexPositions(bracingPts);
      for (const zPoint &p : bracingPts) accumulate(p);
    }
    if (!zTsField::computeGridBounds(minBB, maxBB, scaleUp, cellSize, gridBounds, resX, resY)) {
      _Status = "failed: degenerate bounding box of the unrolled layers.";
      printDebug("[computeFieldImpl] " + _Status);
      return false;
    }
  }
  {
    zFnMeshScalarField fnField(o_ScalarField);
    fnField.clear();
    fnField.create(gridBounds.min, gridBounds.max, resX, resY, 1, true, false);
  }
  const int gridVertexCount = zFnMesh(o_ScalarField).numVertices();

  o_PathFields.assign(totalLayers, zScalarArray());
  o_WidthFields.assign(totalLayers, zScalarArray());
  o_TrackFields.assign(totalLayers, zScalarArray());
  o_FaceTagSources.assign(totalLayers, zTsField::FaceTagSource());
  std::vector<int> segmentCounts(totalLayers, 0);
  std::vector<char> layerPrinted(totalLayers, 0), layerSlotted(totalLayers, 0);

  // The grid is only read, so every thread shares it.
  zTsDetail::ParallelErrors parallelErrors2(totalLayers);
#pragma omp parallel for schedule(dynamic) if(useParallel)
  for (int k = 0; k < totalLayers; ++k) {
    try {
      zObjGraph &boundary2D = o_IsoMeshBoundary2D[k];
      if (zFnGraph(boundary2D).numEdges() == 0 || zFnMesh(o_IsoMeshes2D[k]).numPolygons() == 0) {
        if (i_Settings.DoubleBoundary)
          printDebug("[computeFieldImpl] layer=" + std::to_string(k) + " BoundarySlot skipped: EmptyBoundaryOrChart.");
        // Nothing to print on this layer: a field with no zero-isocontour.
        o_PathFields[k].assign(gridVertexCount, 1.0);
        o_WidthFields[k].assign(gridVertexCount, 0.0);
        o_TrackFields[k].assign(gridVertexCount, 0.0);
        continue;
      }
      const bool isEven = (k % 2 == 0);
      const zIntArray &edgeIds = o_IsoMeshBoundaryWidthId[k];

      // Distances use every actual side edge, including open components.
      std::vector<zTsField::FieldSegment> segments;
      segmentCounts[k] = zTsField::loopWidthSegments(boundary2D, edgeIds, i_Settings.BoundaryWidth, segments);
      if (segments.empty()) {
        if (i_Settings.DoubleBoundary)
          printDebug("[computeFieldImpl] layer=" + std::to_string(k) + " BoundarySlot skipped: NoBoundarySegments.");
        o_PathFields[k].assign(gridVertexCount, 1.0);
        o_WidthFields[k].assign(gridVertexCount, 0.0);
        o_TrackFields[k].assign(gridVertexCount, 0.0);
        continue;
      }

      // Boundary slot: the anchor's 2D image re-located on this layer's loop,
      // +/- SlotOffset * w by parity, w = width of the edge there. The
      // 2D loop turns like the world loop seen from +Z (computeContours()), so it
      // gives the sense itself.
      zObjGraph boundarySlot;
      float slotWidth = 0.0f;
      std::string slotResult = !hasSlotAnchor ? slotAnchorFailure : "LayerAnchorMappingFailed";
      if (slotValid[k]) {
        zPointArray loop2D;
        zIntArray loopEdges;
        zTsField::LoopAnchor anchor;
        slotResult = "NoUsableLayerLoop: boundary is open, branched, empty or invalid";
        if (zTsField::locateGraphLoopAnchor(boundary2D, slotPoints2D[k], loop2D, loopEdges, anchor)) {
          const int edge = loopEdges[anchor.edge];
          slotWidth = zTsField::widthOfGroup(
              i_Settings.BoundaryWidth, (edge < static_cast<int>(edgeIds.size())) ? edgeIds[edge] : 0);
          layerSlotted[k] = zTsField::createParitySlot(loop2D, anchor, isEven, i_Settings.SlotOffset * slotWidth,
                                                       static_cast<zScalar>(slotWidth), boundarySlot);
          slotResult = layerSlotted[k] ? "Created" : "CreationFailed: invalid width, anchor or tangent";
        }
      }

      if (i_Settings.DoubleBoundary && !layerSlotted[k])
        printDebug("[computeFieldImpl] layer=" + std::to_string(k) + " BoundarySlot: " + slotResult);

      if (zFnGraph(o_BracingGraphs2D[k]).numEdges() > 0)
        if (!zTsField::appendBracingSegments(o_BracingGraphs2D[k], i_Settings,
              i_Settings.BracingSlot, segments, k)) continue;

      if (std::any_of(i_FaceTags.begin(), i_FaceTags.end(), [](int t){ return t > 3; })) {
        auto &source = o_FaceTagSources[k]; source.boundary = boundary2D;
        source.slot = boundarySlot; source.slotWidth = slotWidth;
        source.edgeTags = o_IsoMeshBoundaryFaceTags[k]; source.edgeWidthIds = edgeIds;
        source.widths = i_Settings.BoundaryWidth; source.doubleBoundary = i_Settings.DoubleBoundary;
        source.overlap = i_Settings.OverlapWidth;
        for (const auto &segment : segments) if (segment.type != 0) source.bracing.push_back(segment);
      }
      std::vector<zScalarArray> segDistances;
      zScalarArray boundarySigns;
      zTsField::computeSegmentDistances(o_ScalarField, segments, segDistances);
      // Filtered side edges may be open. Only the complete chart defines inside;
      // non-side mesh boundaries contribute no printed boundary-distance segment.
      zTsField::computeMeshDomainSigns(o_ScalarField, o_IsoMeshes2D[k], boundarySigns);
      zTsField::computeWidthField(o_ScalarField, segDistances, segments, i_Settings.DoubleBoundary,
                                  i_Settings.OverlapWidth, &o_WidthFields[k]);
      zTsField::computePathField(o_ScalarField, boundary2D, segments, segDistances, i_Settings.DoubleBoundary,
                                 i_Settings.OverlapWidth, nullptr, &o_PathFields[k], &o_TrackFields[k],
                                 &boundarySlot, static_cast<zScalar>(slotWidth), &boundarySigns);
      layerPrinted[k] = 1;
    } catch (...) { parallelErrors2.capture(k); }
  }
  parallelErrors2.rethrow();

  int printed = 0, slotted = 0, maxSegments = 0, bracingLayers = 0;
  for (int k = 0; k < totalLayers; ++k) {
    printed += layerPrinted[k];
    slotted += layerSlotted[k];
    maxSegments = (std::max)(maxSegments, segmentCounts[k]);
    if (zFnGraph(o_BracingGraphs2D[k]).numEdges() > 0) ++bracingLayers;
  }

  const bool fielded = printed > 0;
  std::ostringstream out;
  out << "per-layer fields on " << resX << "x" << resY << " cells (cell size=" << cellSize << "): " << printed << "/"
      << totalLayers << " layers, up to " << maxSegments << " boundary segments per layer";
  if (hasBracing) out << ", bracing on " << bracingLayers << " layers";
  if (i_Settings.DoubleBoundary) out << ", boundary slot on " << slotted << " layers";
  if (printed < totalLayers) out << " | WARNING: " << (totalLayers - printed) << " layer(s) without usable side edges or a chart";
  _Status = out.str();
  if (!fielded || printed < totalLayers) printDebug("[computeFieldImpl] " + _Status);

  return fielded;
}

ZSPACE_INLINE bool HarmonicSlicer::computePrintPath(bool useParallel) {
  const bool success = computePrintPathImpl(useParallel);
  if (success) updatePrintStatistics(false);
  else resetPrintStatistics();
  if (success) zTsUtils::printSuccess("[HarmonicSlicer::computePrintPath] Success");
  else zTsUtils::printError("[HarmonicSlicer::computePrintPath] Failed");
  return success;
}

ZSPACE_INLINE bool HarmonicSlicer::computePrintPathImpl(bool useParallel) {
  o_PrintPaths2D.clear();
  o_PrintPaths3D.clear();
  o_ResultantPrintPath = PrintPath();
  resetPrintStatistics();
  const int n = int(o_PathFields.size());
  if (n == 0 || o_WidthFields.size() != size_t(n) || o_TrackFields.size() != size_t(n) ||
      o_IsoMeshes2D.size() != size_t(n) || o_IsoMeshes.size() != size_t(n) ||
      o_VertexLayerHeights.size() != size_t(n)) {
    _Status = "failed: fields and charts must have one entry per layer. Call computeField() first.";
    zTsUtils::printDebug("[computePrintPathImpl] " + _Status);
    return false;
  }
  for (int k = 0; k < n; ++k) {
    if (o_VertexLayerHeights[k].size() != size_t(zFnMesh(o_IsoMeshes[k]).numVertices())) {
      _Status = "failed: harmonic layer heights must match their section vertices. Recompute contours.";
      return false;
    }
  }
  o_PrintPaths2D.resize(n);
  o_PrintPaths3D.resize(n);
  std::vector<char> failed(n, 0);
  std::vector<std::vector<std::string>> pathDebug(n);
  // Flush any lazy face-list synchronization before workers share the grid.
  (void)zFnMesh(o_ScalarField).numVertices();
  zTsDetail::ParallelErrors parallelErrors3(n);
#pragma omp parallel for schedule(dynamic) if(useParallel)
  for (int k = 0; k < n; ++k) {
    try {
      if (!zTsField::extractPrintPaths(o_ScalarField, o_PathFields[k], o_WidthFields[k],
                                      o_TrackFields[k], o_PrintPaths2D[k], k, &pathDebug[k])) {
        failed[k] = 1;
        continue;
      }
      if (k < int(o_FaceTagSources.size()) && !o_FaceTagSources[k].edgeTags.empty())
        for (auto &path : o_PrintPaths2D[k])
          if (!zTsField::labelFaceTagSegments(path, o_FaceTagSources[k], false)) throw std::runtime_error("Invalid FaceTag provenance");
      if (!zTsField::mapPrintPaths(o_PrintPaths2D[k], o_IsoMeshes2D[k], o_IsoMeshes[k],
                                  &o_VertexLayerHeights[k], o_PrintPaths3D[k], k)) failed[k] = 1;
    } catch (...) { parallelErrors3.capture(k); }
  }
  parallelErrors3.rethrow();
  for (int k = 0; k < n; ++k) {
    for (const std::string &message : pathDebug[k]) {
      const size_t end = message.find(']');
      zTsUtils::printDebug(message.substr(0, end + 1) + " layer=" + std::to_string(k) + message.substr(end + 1));
    }
  }
  std::string failedLayers;
  for (int k = 0; k < n; ++k)
    if (failed[k]) failedLayers += (failedLayers.empty() ? "" : ",") + std::to_string(k);
  _Status = "Processed print paths for " + std::to_string(n) + " layers.";
  if (!failedLayers.empty()) {
    _Status += " Extraction or mapping failed on layer(s): " + failedLayers;
    return false;
  }
  return true;
}

} // namespace zSpace
