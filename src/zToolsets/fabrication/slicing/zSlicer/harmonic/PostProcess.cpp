#ifndef EIGEN_DONT_PARALLELIZE
#define EIGEN_DONT_PARALLELIZE
#endif

#include <zspace/zToolsets/fabrication/slicing/zSlicer/HarmonicSlicer.h>

#ifndef ZSPACE_INLINE
#define ZSPACE_INLINE ZSPACE_TOOLSETS_INLINE
#endif

#include <sstream>
#include <string>
#include <vector>

namespace zSpace {

using namespace zTsUtils;

//--------------------------
//---- PUBLIC METHODS: zSlicer
//--------------------------

ZSPACE_INLINE bool HarmonicSlicer::computeSlicer(bool useParallel) {
  const bool success = computeSlicerImpl(useParallel);
  if (success) zTsUtils::printSuccess("[HarmonicSlicer::computeSlicer] Success");
  else zTsUtils::printError("[HarmonicSlicer::computeSlicer] Failed");
  return success;
}

ZSPACE_INLINE bool HarmonicSlicer::computeSlicerImpl(bool useParallel) {
  if (!computeBasicStructure(useParallel)) return false;
  if (!computeHarmonicField()) return false;
  if (!computeVolumetricParameterization(useParallel)) return false;
  if (!computeContours(useParallel)) return false;
  return computeField(useParallel);
}

// Every layer has its own path/width/track field on the shared o_ScalarField
// grid, in its own o_IsoMeshes2D chart. The layers share no topology, so the seam is
// carried between them in world coordinates (zPrintSyn's charts/sections
// overload of alignStartPoints()), and every layer's 2D paths go to world
// through its own o_IsoMeshes2D / o_IsoMeshes pair, interpolating that layer's
// vertex heights after resampling and seam alignment.
//
// Stops at the first zPrintSyn step that fails or reports a problem on any
// layer: o_PrintPaths2D keeps the paths as that step left them, projected to
// o_PrintPaths3D for inspection (no o_ResultantPrintPath), _Status names the
// step and its layers, and false is returned.
ZSPACE_INLINE bool HarmonicSlicer::computePostProcess(const PrintSynSetting &settings, bool useParallel) {
  const bool success = computePostProcessImpl(settings, useParallel);
  if (success) updatePrintStatistics(true);
  else resetPrintStatistics();
  if (success) zTsUtils::printSuccess("[HarmonicSlicer::computePostProcess] Success");
  else zTsUtils::printError("[HarmonicSlicer::computePostProcess] Failed");
  return success;
}

ZSPACE_INLINE bool HarmonicSlicer::computePostProcessImpl(const PrintSynSetting &settings, bool useParallel) {
  if (!computePrintPath(useParallel)) return false;

  const bool ok = zPrintSyn::computePipeline(o_IsoMeshes2D, o_IsoMeshes, nullptr,
      zPrintSyn::SeamMapping::WorldNearest, settings, o_PrintPaths2D, o_PrintPaths3D,
       o_ResultantPrintPath, _Status, useParallel, &o_VertexLayerHeights, &o_VertexLayerNormals);
  if (!ok) zTsUtils::printDebug("[computePostProcessImpl] Post-processing stopped because a helper failed or reported invalid paths; partial results are retained.");
  return ok;
}

} // namespace zSpace
