#ifndef ZSPACE_TS_GEOMETRY_PRINTSYN_H
#define ZSPACE_TS_GEOMETRY_PRINTSYN_H

#pragma once

#include <zspace/zToolsets/zToolsetsExport.h>
#include <zspace/zToolsets/fabrication/slicing/utils/Utils.h>
#include <zspace/interface.h>
#include <zspace/io.h>
#include <zspace/zInterface/functionsets/zFnGraph.h>
#include <zspace/zInterface/functionsets/zFnMesh.h>
#include <zspace/zInterface/functionsets/zFnMeshField.h>
#include <zspace/zInterface/objects/zObjGraph.h>
#include <zspace/zInterface/objects/zObjMesh.h>
#include <zspace/zInterface/objects/zObjMeshField.h>

#include <functional>
#include <string>
#include <vector>

namespace zSpace {

/** \addtogroup zToolsets
 *	\brief Toolsets for geometry, fabrication and computation.
 *  @{
 */

/*! \struct PrintSynSetting
 *	\brief Configuration settings for 2D toolpath cleanup and resampling.
 */
struct PrintSynSetting {
  // Sampling uses 2D chart arc length; projection preserves those targets.
  // Tagged paths share per-FaceTag count budgets across layers. The production
  // pipeline does not run a second resampling pass after mapping into 3D.
  float ResampleLength = 0.1f;
  // Retained for settings/layout compatibility. The restored sharp-feature
  // sampler does not use this setting or add curvature-support targets.
  double SurfaceChordTolerance = 0.1;
  // Angle between the incident edges, in degrees. Strictly less than this
  // threshold selects an inner reconstruction/sampling support; equality is
  // excluded (straight = 180). Uniform outer runs and junction feature selection
  // do not add targets at these geometric supports.
  float FeatureAngleTol = 90.0f;
  float SeamOffset = 0.0f;
  float SeamSize = 0.0f;
  // Independent nominal chart-arc spacing for inner/untagged runs, in mm.
  // Appended to preserve all preceding field offsets; rebuild native consumers.
  float InnerResampleLength = 20.0f;
};

/*! \class zPrintSyn
 *	\brief Stateless toolpath synthesis - a method class of static functions
 * only (never instantiated): per-layer SDF fields -> 2D print paths -> 3D
 * print paths -> one bottom-to-top print path.
 *
 *  Every input is passed in and only read for the duration of the call;
 *  every result goes to an out parameter. Each slicer operator's
 *  computePostProcess() runs computePipeline() on its own paths/charts/
 *  sections, selecting the seam mapping strategy. The stage order is:
 *    resamplePrintPaths -> alignStartPoints ->
 *    projectPrintPaths -> linkPrintPaths
 *  FaceTag paths use fixed-count sampling after seam alignment instead of
 *  generic 2D resampling. Neither route resamples the projected 3D targets.
 *  Each method prints its own "[zPrintSyn::...]" success/failure line and
 *  writes the same message to outStatus. outBadLayers (optional) gets the
 *  layers the step reports a problem on (sorted; empty if none) - a step can
 *  succeed with bad layers, see each method.
 *	\since 1.0.0
 */
class ZSPACE_TOOLSETS zPrintSyn {
public:
  zPrintSyn() = delete;

  // Current geometric corners (open endpoints excluded), using FeatureAngleTol.
  static zIntArray cornerVertexIds(const PrintPath &path, float angleThreshold);
  // Distinct incident outer side tags (3..6), using actual graph adjacency.
  // Track metadata, when present, must mark the vertex outer (>=0.5). Angle
  // threshold is retained for source/ABI compatibility and is ignored.
  static zIntArray featureVertexIds(const PrintPath &path, float fallbackAngleThreshold);
  // Replace featureFlags with aligned OuterJunctionFeature/0 metadata. Invalid
  // edge/track arrays return false without modifying the existing metadata.
  static bool refreshJunctionFeatures(PrintPath &path);
  // Reconstruct XY corner groups from two non-corner vertices on each side.
  // Parallel/degenerate groups or groups without enough support are preserved.
  // Production uses this only for inner/untagged runs, never outer runs.
  static bool reconstructCorners(PrintPath &path, float angleThreshold);


  enum class SeamMapping { SharedTopology, WorldNearest };

  // Runs the four post-processing stages, stopping on failure or bad layers. Partial paths
  // are retained for inspection; no combined path is exposed on failure.
  // perLayerVertexHeights overrides vertexLayerHeights for independent layer meshes.
  static bool computePipeline(const zObjMeshArray &charts2D,
                               const zObjMeshArray &sections3D,
                               const zFloatArray *vertexLayerHeights,
                               SeamMapping seamMapping, const PrintSynSetting &settings,
                               PrintPathArray &outPaths2D, PrintPathArray &outPaths3D,
                               PrintPath &outCombined, std::string &outStatus, bool useParallel = true,
                               const std::vector<zFloatArray> *perLayerVertexHeights = nullptr,
                               const std::vector<zVectorArray> *perLayerVertexNormals = nullptr,
                               const std::vector<zTsField::FaceTagSource> *faceTagSources = nullptr);

  // Same production result; optional owned COPY of the successfully mapped,
  // seam-trimmed processed 2D contours corresponding to the final 3D targets. This is
  // internal validation geometry, never an additional fabrication target list.
  // The pointer is explicit to keep existing calls and their ABI unambiguous.
  // Its output is cleared on entry and on any pipeline failure.
  static bool computePipeline(const zObjMeshArray &charts2D,
                               const zObjMeshArray &sections3D,
                               const zFloatArray *vertexLayerHeights,
                               SeamMapping seamMapping, const PrintSynSetting &settings,
                               PrintPathArray &outPaths2D, PrintPathArray &outPaths3D,
                               PrintPath &outCombined, std::string &outStatus, bool useParallel,
                               const std::vector<zFloatArray> *perLayerVertexHeights,
                               const std::vector<zVectorArray> *perLayerVertexNormals,
                               PrintPathArray *outSurfaceContours2D,
                               const std::vector<zTsField::FaceTagSource> *faceTagSources = nullptr);

  // Geometric corners and categorical transitions closer than 0.01 mm along
  // the chart or printable section coalesce before global budget allocation.
  // Exact corner coordinates win over transition coordinates. Dense ordinary
  // contour vertices define the reference arc. Each outer piece is sampled
  // uniformly from endpoint to endpoint with shared per-tag count budgets.
  // Tag0 runs reconstruct privately with fixed shared endpoints, then use
  // independent InnerResampleLength steps; their counts may vary across layers.
  static bool resampleFaceTaggedPaths(PrintPathArray &paths, const std::vector<zTsField::FaceTagSource> &sources,
      const PrintSynSetting &settings, std::string &status,
      const zObjMeshArray *charts2D = nullptr, const zObjMeshArray *sections3D = nullptr);
  // Independent layer preparation may run in parallel; the global per-tag
  // budget reduction retains layer/source order. The legacy overload is serial.
  static bool resampleFaceTaggedPaths(PrintPathArray &paths, const std::vector<zTsField::FaceTagSource> &sources,
      const PrintSynSetting &settings, std::string &status,
      const zObjMeshArray *charts2D, const zObjMeshArray *sections3D, bool useParallel);

  /*! \brief Step 1: resamples (corner-preserving), glitch-filters and
   * orients every 2D print path in place. Untagged paths reconstruct corners
   * and use InnerResampleLength; labeled outer paths do not reconstruct.
   * Only outer junctions are features. Merges targets closer than 0.01 mm.
   * Bad layers: a first path that was
   * closed before and is open after. */
  static bool resamplePrintPaths(PrintPathArray &ioPaths2D,
                                 const PrintSynSetting &settings,
                                 std::string &outStatus,
                                 std::vector<int> *outBadLayers = nullptr, bool useParallel = true);

  /*! \brief Step 2: aligns every layer's ring start point to the seed
   * layer's start (carried barycentrically through the align charts) and
   * cuts the seam gap. alignCharts2D: one chart per layer, all sharing ONE
   * topology. Bad layers: the seed layer if its start cannot be placed, and
   * every layer that keeps its unaligned start. */
  static bool alignStartPoints(PrintPathArray &ioPaths2D,
                               const zObjMeshArray &alignCharts2D,
                               double seamOffset, double seamSize,
                               std::string &outStatus,
                               std::vector<int> *outBadLayers = nullptr, bool useParallel = true);

  /*! \brief Step 3, single shared align chart (every layer uses
   * alignChart2D). */
  static bool alignStartPoints(PrintPathArray &ioPaths2D,
                               const zObjMesh &alignChart2D, double seamOffset,
                               double seamSize, std::string &outStatus,
                               std::vector<int> *outBadLayers = nullptr, bool useParallel = true);

  /*! \brief Step 3 for layers that share no topology: the seed layer's start
   * is lifted to world coordinates through its charts2D/sections3D pair
   * (same vertex order and faces - projectPrintPaths()' pair), located on
   * every layer's 3D section (nearest point) and brought back to that
   * layer's chart the same way. */
  static bool alignStartPoints(PrintPathArray &ioPaths2D,
                               const zObjMeshArray &charts2D,
                               const zObjMeshArray &sections3D,
                               double seamOffset, double seamSize,
                               std::string &outStatus,
                               std::vector<int> *outBadLayers = nullptr, bool useParallel = true);

  /*! \brief Step 3: projects each layer's 2D paths through its charts2D[l] /
   * sections3D[l] pair (same vertex order and faces) onto the 3D section.
   * vertexLayerHeights (optional, may be null) is passed straight to
   * zTsUtils::remapGraph(). perLayerVertexHeights, when supplied, overrides it
   * with one table per section, with one finite nonnegative height per vertex.
   * perLayerVertexNormals optionally supplies recovered normals in that same
   * section-vertex order; empty tables preserve face normals (e.g. fixed caps).
   * Preserves vertex order/count and featureFlags. Returns false on any
   * unmappable path; partial valid paths remain available for inspection.
   * Bad layers: a path that could not be remapped. */
  static bool projectPrintPaths(const PrintPathArray &paths2D,
                                const zObjMeshArray &charts2D,
                                const zObjMeshArray &sections3D,
                                const zFloatArray *vertexLayerHeights,
                                PrintPathArray &outPaths3D,
                                std::string &outStatus,
                                std::vector<int> *outBadLayers = nullptr, bool useParallel = true,
                                const std::vector<zFloatArray> *perLayerVertexHeights = nullptr,
                                const std::vector<zVectorArray> *perLayerVertexNormals = nullptr);

  /*! \brief Step 4: combines all 3D paths. Single-path layers are linked
   * bottom-to-top; disjoint curves remain separate without explicit travel. */
  static bool linkPrintPaths(const PrintPathArray &paths3D,
                             PrintPath &outCombined, std::string &outStatus);

  // Explicit standalone 3D resampling utility; not called by computePipeline.
  // After reconstructed corners are mapped, keeps
  // exact mandatory targets; ordinary points have half-spacing clearance from
  // adjacent features/endpoints and are placed on the current section. Updates
  // the matching 2D targets/attributes using the same surface barycentrics.
  static bool resamplePrintPaths3D(PrintPathArray &paths2D, PrintPathArray &paths3D,
      const zObjMeshArray &charts2D, const zObjMeshArray &sections3D,
      const PrintSynSetting &settings, std::string &status, bool useParallel=true,
      const zFloatArray *vertexHeights=nullptr,
      const std::vector<zFloatArray> *perLayerHeights=nullptr,
      const std::vector<zVectorArray> *perLayerNormals=nullptr);

private:
  // Carries the seed layer's start point to every layer: outMapped[l] /
  // outValid[l] per layer; false if the seed cannot be located at all.
  using SeedMapper = std::function<bool(int seedLayer, const zPoint &seed,
                                        zPointArray &outMapped,
                                        std::vector<char> &outValid)>;

  // Shared body of every alignStartPoints() overload. The per-path helpers
  // live in zPrintSyn.cpp.
  static bool alignStartPointsWith(PrintPathArray &ioPaths,
                                   const SeedMapper &mapSeed,
                                   double seamOffset, double seamSize,
                                   std::string &outStatus,
                                   std::vector<int> *outBadLayers, bool useParallel = true);

  // The two align-chart overloads (alignCharts[l] is layer l's chart, all
  // sharing one topology).
  static bool alignStartPointsOnCharts(PrintPathArray &ioPaths,
                                       const std::vector<const zObjMesh *> &alignCharts,
                                       double seamOffset, double seamSize,
                                       std::string &outStatus,
                                       std::vector<int> *outBadLayers, bool useParallel = true);
};

/** @}*/

} // namespace zSpace


#endif // ZSPACE_TS_GEOMETRY_PRINTSYN_H
