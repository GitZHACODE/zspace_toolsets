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

#include <zspace/zToolsets/fabrication/slicing/zSlicer/TopologySlicer.h>

#ifndef ZSPACE_INLINE
#define ZSPACE_INLINE ZSPACE_TOOLSETS_INLINE
#endif

#include <string>
#include <vector>

namespace zSpace {

//--------------------------
//---- PUBLIC METHODS: zSlicer
//--------------------------

ZSPACE_INLINE bool TopologySlicer::computeSlicer(bool useParallel) {
  zTsUtils::ScopedSlicingTiming methodTiming("TopologySlicer.computeSlicer");
  const bool success = computeSlicerImpl(useParallel);
  if (success) zTsUtils::printSuccess("[TopologySlicer::computeSlicer] Success");
  else zTsUtils::printError("[TopologySlicer::computeSlicer] Failed");
  return success;
}

ZSPACE_INLINE bool TopologySlicer::computeSlicerImpl(bool useParallel) {
  o_PrintPaths2D.clear();
  o_PrintPaths3D.clear();
  o_ResultantPrintPath = PrintPath();
  resetPrintStatistics();
  if (!computeBasicStructure(useParallel))
    return false;
  if (!computeSlice(useParallel))
    return false;
  if (!computeUnroll(useParallel))
    return false;
  return computeField(useParallel);
}

// The unrolled section meshes share the U-mesh's topology, so they are the
// align charts and, paired 1:1 with the 3D section meshes, the projection
// charts at once; every layer has its own path/width/track field on the
// shared o_ScalarField grid.
//
// Stops at the first zPrintSyn step that fails or reports a problem on any
// layer: o_PrintPaths2D keeps the paths as that step left them, projected to
// o_PrintPaths3D for inspection (no o_ResultantPrintPath), _Status names the
// step and its layers, and false is returned.
ZSPACE_INLINE bool TopologySlicer::computePostProcess(const PrintSynSetting &settings, bool useParallel) {
  zTsUtils::ScopedSlicingTiming methodTiming("TopologySlicer.computePostProcess");
  const bool success = computePostProcessImpl(settings, useParallel);
  if (success) updatePrintStatistics(true);
  else resetPrintStatistics();
  if (success) zTsUtils::printSuccess("[TopologySlicer::computePostProcess] Success");
  else zTsUtils::printError("[TopologySlicer::computePostProcess] Failed");
  return success;
}

ZSPACE_INLINE bool TopologySlicer::computePostProcessImpl(const PrintSynSetting &settings, bool useParallel) {
  if (!computePrintPath(useParallel)) return false;

  const int offset = getPrintSectionOffset();
  if (offset < 0 || size_t(offset) >= o_IsoMeshes.size() ||
      o_IsoMeshes2D.size() != o_IsoMeshes.size() ||
      o_FaceTagSources.size() != o_IsoMeshes.size()) {
    _Status = "failed: inconsistent printable TCP section window.";
    return false;
  }
  // Section zero is the substrate. Every path uses the matching section above
  // it, including the exact top section, with provenance sliced in lockstep.
  const zObjMeshArray charts(o_IsoMeshes2D.begin() + offset, o_IsoMeshes2D.end());
  const zObjMeshArray sections(o_IsoMeshes.begin() + offset, o_IsoMeshes.end());
  const std::vector<zTsField::FaceTagSource> sources(o_FaceTagSources.begin() + offset,
                                                  o_FaceTagSources.end());
  const bool ok = zPrintSyn::computePipeline(charts, sections, &_VLength,
      zPrintSyn::SeamMapping::SharedTopology, settings, o_PrintPaths2D, o_PrintPaths3D,
      o_ResultantPrintPath, _Status, useParallel, nullptr, nullptr, &sources);
  if (!ok) zTsUtils::printDebug("[computePostProcessImpl] Post-processing stopped because a helper failed or reported invalid paths; partial results are retained.");
  return ok;
}

} // namespace zSpace
