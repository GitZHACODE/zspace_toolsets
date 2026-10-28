#ifndef EIGEN_DONT_PARALLELIZE
#define EIGEN_DONT_PARALLELIZE
#endif

#include <zspace/zToolsets/fabrication/slicing/utils/Parallel.h>
#include <zspace/zToolsets/fabrication/slicing/zSlicer/TopologySlicer.h>
#include <zspace/zToolsets/fabrication/slicing/utils/Utils.h>

#ifndef ZSPACE_INLINE
#define ZSPACE_INLINE ZSPACE_TOOLSETS_INLINE
#endif

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <limits>
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
//---- PRIVATE HELPERS
//--------------------------
// The SDF field kernels (segment distances, width/path fields, trim slots)
// are shared with HarmonicSlicer - see Utils.h.

namespace {

// Computes a cell-aligned bounding rectangle and resolution for the unrolled
// meshes.
bool _getboundingrect(const zObjMeshArray &unrollMeshes,
                      zDomain<zPoint> &outBounds, double scaleFactor = 1.00,
                      double targetCellSize = 0.05, int *outResX = nullptr,
                      int *outResY = nullptr) {
  if (unrollMeshes.empty())
    return false;

  zPoint minBB((std::numeric_limits<double>::max)(),
               (std::numeric_limits<double>::max)(), 0.0);
  zPoint maxBB(-(std::numeric_limits<double>::max)(),
               -(std::numeric_limits<double>::max)(), 0.0);

  for (const auto &mesh : unrollMeshes) {
    zFnMesh fn(const_cast<zObjMesh &>(mesh));
    if (fn.numVertices() == 0)
      continue;
    zPoint meshMin, meshMax;
    fn.getBounds(meshMin, meshMax);
    minBB.x = (std::min)(minBB.x, meshMin.x);
    minBB.y = (std::min)(minBB.y, meshMin.y);
    maxBB.x = (std::max)(maxBB.x, meshMax.x);
    maxBB.y = (std::max)(maxBB.y, meshMax.y);
  }

  int resX = 0, resY = 0;
  if (!zTsField::computeGridBounds(minBB, maxBB, scaleFactor, targetCellSize,
                                   outBounds, resX, resY))
    return false;
  if (outResX)
    *outResX = resX;
  if (outResY)
    *outResY = resY;
  return true;
}

} // anonymous namespace

//--------------------------
//---- PUBLIC METHODS: STEP 3 (2D SDF SCALAR FIELD)
//--------------------------

// Builds the SDF path/width/track fields per layer. Stops here - isocontour
// extraction into PrintPaths is zPrintSyn::extractTrackPaths now.
ZSPACE_INLINE bool TopologySlicer::computeField(bool useParallel) {
  const bool success = computeFieldImpl(useParallel);
  if (success) zTsUtils::printSuccess("[TopologySlicer::computeField] Success");
  else zTsUtils::printError("[TopologySlicer::computeField] Failed");
  return success;
}

ZSPACE_INLINE bool TopologySlicer::computeFieldImpl(bool useParallel) {
  if (!zTsField::validateFieldSettings(i_Settings)) {
    _Status = "failed: invalid physical field parameters (minimum 1e-4 mm).";
    return false;
  }
  const zObjMeshArray &sectionMeshes2D = o_IsoMeshes2D;
  zObjMeshScalarField &outScalarField = o_ScalarField;
  std::vector<zScalarArray> &outPathFields = o_PathFields;
  std::vector<zScalarArray> &outWidthFields = o_WidthFields;
  std::vector<zScalarArray> &outTrackFields = o_TrackFields;
  outPathFields.clear();
  outWidthFields.clear();
  outTrackFields.clear();
  o_FaceTagSources.clear();

  if (!i_SlotPoint) {
    _Status = "failed: i_SlotPoint is required. Call setSlotAnchor() before computing.";
    printDebug("[computeFieldImpl] " + _Status);
    return false;
  }

  const int totalLayers = int(sectionMeshes2D.size());
  if (totalLayers == 0 || o_IsoMeshBoundary2D.empty()) {
    _Status =
        "failed: No unrolled 2D meshes available. Call computeUnroll() first.";
    printDebug("[computeFieldImpl] " + _Status);
    return false;
  }

  outPathFields.assign(totalLayers, zScalarArray());
  outWidthFields.assign(totalLayers, zScalarArray());
  outTrackFields.assign(totalLayers, zScalarArray());
  o_FaceTagSources.resize(totalLayers);

  const float printOverlapWidth = i_Settings.OverlapWidth;
  const float targetCellSize =
      i_Settings.FieldCellSize;

  zDomain<zPoint> layerFieldBB;
  int resX = 0, resY = 0;
  if (!_getboundingrect(sectionMeshes2D, layerFieldBB, i_Settings.FieldScaleUp,
                        targetCellSize, &resX, &resY)) {
    _Status =
        "failed: Could not compute bounding rectangle for unrolled meshes.";
    printDebug("[computeFieldImpl] " + _Status);
    return false;
  }

  const std::vector<float> &edgeWidths = i_Settings.BoundaryWidth;

  // Boundary trim slot (DoubleBoundary): the anchor is located once on the
  // floor section's boundary, transferred by barycentric coordinates
  // through the 3D floor section onto every layer's unrolled chart (the
  // same transfer zPrintSyn::alignStartPoints() uses for the seam), and
  // re-located there on that layer's own boundary; offset +/-
  // SlotOffset * w by layer parity.
  zTsField::LoopAnchor slotAnchor;
  zPointArray slotPoints2D;
  std::vector<char> slotPointValid;
  float slotWidth = 0.0f;
  if (i_Settings.DoubleBoundary) {
    zPointArray floorLoop;
    zFnGraph(o_IsoMeshBoundary[0]).getVertexPositions(floorLoop);
    if (zTsField::locateSlotAnchor(floorLoop, *i_SlotPoint, true,
                                   "TopologySlicer::computeField",
                                   slotAnchor)) {
      // w = width of the floor edge holding the anchor.
      zIntArray floorEdgeIds;
      _layerEdgeWidthIds(0, totalLayers, floorEdgeIds);
      slotWidth = zTsField::widthOfGroup(
          edgeWidths, (slotAnchor.edge < int(floorEdgeIds.size()))
                          ? floorEdgeIds[slotAnchor.edge]
                          : 0);

      std::vector<const zObjMesh *> charts;
      for (int l = 0; l < totalLayers; l++)
        charts.push_back(l < int(sectionMeshes2D.size()) ? &sectionMeshes2D[l]
                                                         : nullptr);
      if (!zTsUtils::mapPointToMeshes(slotAnchor.point, o_IsoMeshes[0],
                                      charts, slotPoints2D, slotPointValid, nullptr, useParallel))
        printDebug("[computeFieldImpl] could not transfer the "
                   "slot anchor onto the layer charts - no boundary slot.");
    }
  }

  // Initialize reusable scalar field topology
  {
    zFnMeshScalarField fnField(outScalarField);
    fnField.clear();
    fnField.create(layerFieldBB.min, layerFieldBB.max, resX, resY, 1, true,
                   false);
  }

  const auto wallStart = std::chrono::high_resolution_clock::now();

  // Multi-threaded layer-by-layer computation (OpenMP); the grid is only
  // read, so every thread shares outScalarField.

  zTsDetail::ParallelErrors parallelErrors1(totalLayers);
#pragma omp parallel for schedule(dynamic) if(useParallel)
  for (int l = 0; l < totalLayers; l++) {
    try {
      zObjGraph unrollBracingGraph;
      bool hasValidBracing = false;
      if (l < int(o_BracingGraphs2D.size())) {
        unrollBracingGraph = o_BracingGraphs2D[l];
        zFnGraph fnBracingCheck(unrollBracingGraph);
        if (fnBracingCheck.numVertices() > 0 && fnBracingCheck.numEdges() > 0) {
          hasValidBracing = true;
        }
      }

      // Boundary segments: runs of the same wall-face width group along
      // this layer's loop (non-adjacent faces of one group stay separate).
      std::vector<zTsField::FieldSegment> segmentGraphs;
      zIntArray layerEdgeIds;
      _layerEdgeWidthIds(l, totalLayers, layerEdgeIds);
      zTsField::loopWidthSegments(o_IsoMeshBoundary2D[l], layerEdgeIds,
                                  edgeWidths, segmentGraphs);

      zObjGraph boundarySlot;
      if (l < int(slotPointValid.size()) && slotPointValid[l]) {
        zPointArray layerLoop2D, layerLoop3D;
        zFnGraph(o_IsoMeshBoundary2D[l]).getVertexPositions(layerLoop2D);
        zFnGraph(o_IsoMeshBoundary[l]).getVertexPositions(layerLoop3D);
        zTsField::LoopAnchor layerAnchor;
        if (zTsField::locateLoopAnchor(layerLoop2D, slotPoints2D[l],
                                       layerAnchor, &layerLoop3D))
          zTsField::createParitySlot(layerLoop2D, layerAnchor, l % 2 == 0,
                                     i_Settings.SlotOffset * slotWidth,
                                     slotWidth, boundarySlot);
      }

      if (hasValidBracing) {
        if (!zTsField::appendBracingSegments(unrollBracingGraph, i_Settings,
              i_Settings.BracingSlot, segmentGraphs, l)) continue;
      }

      if (std::any_of(i_FaceTags.begin(), i_FaceTags.end(), [](int t){ return t > 3; })) {
        auto &source = o_FaceTagSources[l];
        source.boundary = o_IsoMeshBoundary2D[l]; source.slot = boundarySlot;
        _layerEdgeFaceTags(l, totalLayers, source.edgeTags);
        source.edgeWidthIds = layerEdgeIds; source.widths = edgeWidths;
        source.doubleBoundary = i_Settings.DoubleBoundary;
        source.overlap = printOverlapWidth; source.slotWidth = slotWidth;
        for (const auto &segment : segmentGraphs) if (segment.type != 0) source.bracing.push_back(segment);
      }

      std::vector<zScalarArray> segDistances;
      zScalarArray boundarySigns;
      zTsField::computeSegmentDistances(outScalarField, segmentGraphs,
                                        segDistances, &boundarySigns);

      zTsField::computeWidthField(outScalarField, segDistances, segmentGraphs,
                                  i_Settings.DoubleBoundary, printOverlapWidth,
                                  &outWidthFields[l]);

      zTsField::computePathField(outScalarField, o_IsoMeshBoundary2D[l],
                                 segmentGraphs, segDistances,
                                 i_Settings.DoubleBoundary, printOverlapWidth,
                                 nullptr, &outPathFields[l],
                                 &outTrackFields[l], &boundarySlot, slotWidth, &boundarySigns);
    } catch (...) { parallelErrors1.capture(l); }
  }
  parallelErrors1.rethrow();


  const auto wallEnd = std::chrono::high_resolution_clock::now();
  const double wallTotalMs =
      std::chrono::duration<double, std::milli>(wallEnd - wallStart).count();

  std::ostringstream out;
  out << "Compute " << totalLayers << " layers (CellSize: " << std::fixed
      << std::setprecision(2) << targetCellSize << ", Res: " << resX << "x"
      << resY << ") in " << std::setprecision(2) << wallTotalMs << " ms";
  _Status = out.str();

  return true;
}

//--------------------------
//---- PROTECTED METHODS: WIDTH GROUPS
//--------------------------

// Layer `layer` sits at the same normalized arc length on every V-loop
// (computeSlice()); edge j of its loop spans the V-loops of vertices j and
// j + 1 and takes the group of the wall face between them at that height,
// from _VLoopSegmentWidthIds. Loop vertex j is V-loop j unless the layer's
// loop came out reversed (_alignBoundaryVertices() keeps it
// counter-clockwise), which is detected from vertex 1.
ZSPACE_INLINE void TopologySlicer::_layerEdgeWidthIds(int layer, int totalLayers,
                                                     zIntArray &outEdgeIds) const {
  outEdgeIds.clear();
  if (layer < 0 || layer >= int(o_IsoMeshBoundary.size()))
    return;
  zPointArray loop;
  zFnGraph(const_cast<zObjGraph &>(o_IsoMeshBoundary[layer])).getVertexPositions(loop);
  const int n = int(loop.size());
  if (n < 3 || n != int(_VLoops.size()) || n != int(_VLoopSegmentWidthIds.size()))
    return;
  const double ratio = (totalLayers > 1) ? double(layer) / (totalLayers - 1) : 0.0;

  // Group of each V-loop's segment at this height.
  zIntArray columnGroup(n, 0);
  for (int c = 0; c < n; c++) {
    const zScalarArray &arc = _VScalar[c];
    const zIntArray &groups = _VLoopSegmentWidthIds[c];
    if (arc.size() < 2 || groups.empty())
      continue;
    const double s = ratio * arc.back();
    int k = int(std::upper_bound(arc.begin(), arc.end(), s) - arc.begin()) - 1;
    k = (std::max)(0, (std::min)(k, int(groups.size()) - 1));
    columnGroup[c] = groups[k];
  }

  // Loop vertex -> V-loop: identity, or (n - j) % n if the loop was reversed.
  bool reversed = false;
  {
    zPoint onNext, onPrev;
    samplePolyline(_VLoops[1], _VScalar[1], ratio, onNext);
    samplePolyline(_VLoops[n - 1], _VScalar[n - 1], ratio, onPrev);
    reversed = vecDistSq(loop[1], onPrev) < vecDistSq(loop[1], onNext);
  }
  const auto columnOf = [&](int j) { return reversed ? (n - j) % n : j; };

  outEdgeIds.assign(n, 0);
  for (int j = 0; j < n; j++) {
    const int a = columnOf(j), b = columnOf((j + 1) % n);
    // The face strip between V-loops a and a + 1 is recorded on V-loop a.
    outEdgeIds[j] = columnGroup[(b == (a + 1) % n) ? a : b];
  }
}

ZSPACE_INLINE void TopologySlicer::_layerEdgeFaceTags(int layer, int totalLayers,
                                                     zIntArray &outEdgeIds) const {
  outEdgeIds.clear();
  if (layer < 0 || layer >= int(o_IsoMeshBoundary.size()))
    return;
  zPointArray loop;
  zFnGraph(const_cast<zObjGraph &>(o_IsoMeshBoundary[layer])).getVertexPositions(loop);
  const int n = int(loop.size());
  if (n < 3 || n != int(_VLoops.size()) || n != int(_VLoopSegmentFaceTags.size()))
    return;
  const double ratio = (totalLayers > 1) ? double(layer) / (totalLayers - 1) : 0.0;

  // Group of each V-loop's segment at this height.
  zIntArray columnGroup(n, 0);
  for (int c = 0; c < n; c++) {
    const zScalarArray &arc = _VScalar[c];
    const zIntArray &groups = _VLoopSegmentFaceTags[c];
    if (arc.size() < 2 || groups.empty())
      continue;
    const double s = ratio * arc.back();
    int k = int(std::upper_bound(arc.begin(), arc.end(), s) - arc.begin()) - 1;
    k = (std::max)(0, (std::min)(k, int(groups.size()) - 1));
    columnGroup[c] = groups[k];
  }

  // Loop vertex -> V-loop: identity, or (n - j) % n if the loop was reversed.
  bool reversed = false;
  {
    zPoint onNext, onPrev;
    samplePolyline(_VLoops[1], _VScalar[1], ratio, onNext);
    samplePolyline(_VLoops[n - 1], _VScalar[n - 1], ratio, onPrev);
    reversed = vecDistSq(loop[1], onPrev) < vecDistSq(loop[1], onNext);
  }
  const auto columnOf = [&](int j) { return reversed ? (n - j) % n : j; };

  outEdgeIds.assign(n, 0);
  for (int j = 0; j < n; j++) {
    const int a = columnOf(j), b = columnOf((j + 1) % n);
    // The face strip between V-loops a and a + 1 is recorded on V-loop a.
    outEdgeIds[j] = columnGroup[(b == (a + 1) % n) ? a : b];
  }
}

ZSPACE_INLINE bool TopologySlicer::computePrintPath(bool useParallel) {
  const bool success = computePrintPathImpl(useParallel);
  if (success) updatePrintStatistics(false);
  else resetPrintStatistics();
  if (success) zTsUtils::printSuccess("[TopologySlicer::computePrintPath] Success");
  else zTsUtils::printError("[TopologySlicer::computePrintPath] Failed");
  return success;
}

ZSPACE_INLINE bool TopologySlicer::computePrintPathImpl(bool useParallel) {
  o_PrintPaths2D.clear();
  o_PrintPaths3D.clear();
  o_ResultantPrintPath = PrintPath();
  resetPrintStatistics();
  const int n = int(o_PathFields.size());
  if (n == 0 || o_WidthFields.size() != size_t(n) || o_TrackFields.size() != size_t(n) ||
      o_IsoMeshes2D.size() != size_t(n) || o_IsoMeshes.size() != size_t(n)) {
    _Status = "failed: fields and charts must have one entry per layer. Call computeField() first.";
    zTsUtils::printDebug("[computePrintPathImpl] " + _Status);
    return false;
  }
  o_PrintPaths2D.resize(n);
  o_PrintPaths3D.resize(n);
  std::vector<char> failed(n, 0);
  std::vector<std::vector<std::string>> pathDebug(n);
  // Flush any lazy face-list synchronization before workers share the grid.
  (void)zFnMesh(o_ScalarField).numVertices();
  zTsDetail::ParallelErrors parallelErrors2(n);
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
                                  &_VLength, o_PrintPaths3D[k], k)) failed[k] = 1;
    } catch (...) { parallelErrors2.capture(k); }
  }
  parallelErrors2.rethrow();
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
