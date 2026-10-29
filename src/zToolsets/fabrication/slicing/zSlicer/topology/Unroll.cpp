// This file is part of zspace, a simple C++ collection of geometry
// data-structures & algorithms, data analysis & visualization framework.
//
// Copyright (C) 2019 ZSPACE
//
// This Source Code Form is subject to the terms of the MIT License
// If a copy of the MIT License was not distributed with this file, You can
// obtain one at https://opensource.org/licenses/MIT.
//
// Author : Vishu Bhooshan <vishu.bhooshan@zaha-hadid.com>
//

#ifndef EIGEN_DONT_PARALLELIZE
#define EIGEN_DONT_PARALLELIZE
#endif

#include <zspace/zToolsets/fabrication/slicing/utils/Parallel.h>
#include <zspace/zToolsets/fabrication/slicing/zSlicer/TopologySlicer.h>

#ifndef ZSPACE_INLINE
#define ZSPACE_INLINE ZSPACE_TOOLSETS_INLINE
#endif

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <utility>
#include <vector>
#if defined(_OPENMP)
#include <omp.h>
#endif

namespace zSpace {

using namespace zTsUtils;

//--------------------------
//---- PUBLIC METHODS: STEP 2 (2D UNROLLING)
//--------------------------

// Face-by-face unfolding (zTsUnroll): one dual-graph BSF order from layer 0,
// shared by every layer (all have the U-mesh's topology).
ZSPACE_INLINE bool TopologySlicer::computeUnroll(bool useParallel) {
  zTsUtils::ScopedSlicingTiming methodTiming("TopologySlicer.computeUnroll");
  const bool success = computeUnrollImpl(useParallel);
  if (success) zTsUtils::printSuccess("[TopologySlicer::computeUnroll] Success");
  else zTsUtils::printError("[TopologySlicer::computeUnroll] Failed");
  return success;
}

ZSPACE_INLINE bool TopologySlicer::computeUnrollImpl(bool useParallel) {
  const zObjMeshArray &sectionMeshes = o_IsoMeshes;
  zObjMeshArray &outSectionMeshes2D = o_IsoMeshes2D;
  outSectionMeshes2D.clear();
  o_IsoMeshBoundary2D.clear();
  o_BracingGraphs2D.clear();

  const int totalLayers = int(sectionMeshes.size());
  if (totalLayers == 0) {
    _Status = "failed: No section meshes available for unroll.";
    printDebug("[computeUnrollImpl] " + _Status);
    return false;
  }

  zIntPairArray bsfHalfEdgePairs;
  if (zFnMesh(const_cast<zObjMesh &>(sectionMeshes[0])).numPolygons() > 1) {
    if (!zTsUnroll::bsfUnfoldOrder(const_cast<zObjMesh &>(sectionMeshes[0]),
                           bsfHalfEdgePairs)) {
      _Status = "failed: Could not compute dual graph for unroll.";
      printDebug("[computeUnrollImpl] " + _Status);
      return false;
    }
  }

  const auto wallStart = std::chrono::high_resolution_clock::now();

  outSectionMeshes2D.assign(totalLayers, zObjMesh());
  o_IsoMeshBoundary2D.assign(totalLayers, zObjGraph());
  if (!_BracingParams.empty()) {
    o_BracingGraphs2D.assign(totalLayers, zObjGraph());
  }

  zTsDetail::ParallelErrors parallelErrors1(totalLayers);
#pragma omp parallel for schedule(dynamic) if(useParallel)
  for (int l = 0; l < totalLayers; l++) {
    try {
      if (zTsUnroll::unfoldMesh(const_cast<zObjMesh &>(sectionMeshes[l]), bsfHalfEdgePairs,
                      outSectionMeshes2D[l])) {
        zFnMesh fnUnrollMesh(outSectionMeshes2D[l]);
        zPointArray unrollMeshPos;
        fnUnrollMesh.getVertexPositions(unrollMeshPos);

        zFnGraph fn3DGraph(o_IsoMeshBoundary[l]);
        const int N = fn3DGraph.numVertices();
        zIntArray edgeConnects;
        fn3DGraph.getEdgeData(edgeConnects);

        if (N > 0 && N <= (int)unrollMeshPos.size()) {
          zPointArray unrollBoundPositions(unrollMeshPos.begin(),
                                         unrollMeshPos.begin() + N);
          zFnGraph fnUnrollGraph(o_IsoMeshBoundary2D[l]);
          fnUnrollGraph.create(unrollBoundPositions, edgeConnects);

          if (!_BracingParams.empty()) {
            realizeParameterGraph(
                _BracingParams, outSectionMeshes2D[l], o_BracingGraphs2D[l]);
          }
        }
      }
    } catch (...) { parallelErrors1.capture(l); }
  }
  parallelErrors1.rethrow();

  const auto wallEnd = std::chrono::high_resolution_clock::now();
  const double wallTotalMs =
      std::chrono::duration<double, std::milli>(wallEnd - wallStart).count();

  std::ostringstream out;
  out << "Computed unroll: " << outSectionMeshes2D.size() << " layers in "
      << std::fixed << std::setprecision(2) << wallTotalMs << " ms";
  _Status = out.str();

  return true;
}

} // namespace zSpace
