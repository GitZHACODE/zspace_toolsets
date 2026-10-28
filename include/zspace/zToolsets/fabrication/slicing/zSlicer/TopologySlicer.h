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

#ifndef ZSPACE_TS_3DP_TOPOLOGY_SLICER_H
#define ZSPACE_TS_3DP_TOPOLOGY_SLICER_H

#pragma once

#ifndef EIGEN_DONT_PARALLELIZE
#define EIGEN_DONT_PARALLELIZE
#endif

#include <zspace/zToolsets/zToolsetsExport.h>
#include <zspace/zToolsets/fabrication/slicing/utils/Utils.h>
#include <zspace/zToolsets/fabrication/slicing/zSlicer/zSlicer.h>

#include <zspace/interface.h>
#include <zspace/io.h>
#include <zspace/zInterface/functionsets/zFnMeshField.h>
#include <zspace/zInterface/objects/zObjGraph.h>
#include <zspace/zInterface/objects/zObjMesh.h>
#include <zspace/zInterface/objects/zObjMeshField.h>


#include <string>
#include <vector>

namespace zSpace {

/** \addtogroup zToolsets
 *	\brief Collection of toolsets for geometry, analysis, and fabrication.
 *  @{
 */

/** \addtogroup zTsGeometry
 *	\brief Toolsets for geometry generation, slicing, and manipulation.
 *  @{
 */

/*! \class TopologySlicer
 *	\brief Topology-driven slicing of a closed quad mesh whose print floor
 * and ceiling are given by vertex tags (setInputMesh(mesh, tags)), plus SDF
 * path/width field generation per layer.
 *
 *  Structure (computeBasicStructure()), with no user-picked start vertex:
 *  - U-mesh: the faces whose vertices are all Bottom-tagged (the floor cap;
 *    one connected patch with a single boundary loop).
 *  - U-loop: the U-mesh boundary, counter-clockwise, starting at the cap
 *    corner (vertex in exactly one cap face) nearest the world origin - model
 *    the intended start corner closest to the origin to control it.
 *  - V-loops: from every U-loop vertex, its one edge leaving the floor, then
 *    straight up the quad edge loop until the first Top-tagged vertex; each
 *    column has its own segment count.
 *  Layers interpolate every column at the same normalized arc length.
 *  U-mesh interior vertices are not supported.
 */
class ZSPACE_TOOLSETS TopologySlicer : public zSlicer {
protected:
  //--------------------------
  //---- PROTECTED ATTRIBUTES
  //--------------------------

  //---- INPUTS
  zObjMesh *i_BlockMesh;
  zIntArray i_FaceTags;
  zIntArray i_FaceWidthIds;
  zObjGraph *i_BracingGraph;
  const zPoint *i_SlotPoint = nullptr;
  SlicerSettings i_Settings;

  //---- INTERNAL
  zIntArray _VertexTags;
  zObjMesh _UMesh;
  zObjGraphArray _VLoops;
  std::vector<zScalarArray> _VScalar;
  zFloatArray _VLength;
  std::vector<zIntArray> _VLoopSegmentWidthIds, _VLoopSegmentFaceTags;
  ParameterGraph _BracingParams;
  std::string _Status;

  //---- OUTPUTS
  zObjMeshArray o_IsoMeshes;
  zObjGraphArray o_IsoMeshBoundary;

  zObjMeshArray o_IsoMeshes2D;
  zObjGraphArray o_IsoMeshBoundary2D;
  zObjGraphArray o_BracingGraphs2D;

  zObjMeshScalarField o_ScalarField;
  std::vector<zScalarArray> o_PathFields;
  std::vector<zScalarArray> o_WidthFields;
  std::vector<zScalarArray> o_TrackFields;
  std::vector<zTsField::FaceTagSource> o_FaceTagSources;

  //---- OUTPUTS : computePostProcess() (cleared by computeSlicer())
  PrintPathArray o_PrintPaths2D;
  PrintPathArray o_PrintPaths3D;
  zDomain<float> o_LayerHeightDomain = zDomain<float>(0.0f, 0.0f);
  PrintPath o_ResultantPrintPath;

public:
  //--------------------------
  //---- CONSTRUCTOR & DESTRUCTOR
  //--------------------------

  TopologySlicer() {
    i_BlockMesh = nullptr;
    i_BracingGraph = nullptr;
    _Status = "Ready";
  }

  explicit TopologySlicer(const SlicerSettings &inSettings) : TopologySlicer() {
    i_Settings = inSettings;
  }

  ~TopologySlicer() {}

  //--------------------------
  //---- SET METHODS
  //--------------------------

  void setSettings(const SlicerSettings &inSettings) override {
    i_Settings = inSettings;
  }
  const SlicerSettings &getSettings() const override { return i_Settings; }
  SlicerSettings &getSettings() override { return i_Settings; }

  void setInputMesh(zObjMesh &inMesh, const zIntArray &inFaceTags,
                    const zIntArray &inFaceWidthIds = zIntArray()) override {
    i_BlockMesh = &inMesh;
    i_FaceTags = inFaceTags;
    i_FaceWidthIds = inFaceWidthIds;
    _VertexTags.clear();
  }
  const std::vector<zTsField::FaceTagSource> &getFaceTagSources() const { return o_FaceTagSources; }
  const zIntArray &getFaceTags() const { return i_FaceTags; }
  const zIntArray &getFaceWidthIds() const { return i_FaceWidthIds; }
  // Derived by computeBasicStructure(); empty before it.
  const zIntArray &getVertexTags() const { return _VertexTags; }
  void setBracingGraph(zObjGraph &inGraph) override {
    i_BracingGraph = &inGraph;
  }
  void setSlotAnchor(const zPoint &inAnchor) override {
    i_SlotPoint = &inAnchor;

  }
  void clearSlotAnchor() override {
    i_SlotPoint = nullptr;

  }

  void setStatus(const std::string &inStatus) { _Status = inStatus; }
  const std::string &getStatus() const override { return _Status; }

  //--------------------------
  //---- GET METHODS
  //--------------------------

  bool isLoaded() const override { return i_BlockMesh != nullptr; }
  int getLayerCount() const { return int(o_IsoMeshBoundary.size()); }

  zObjMesh &getInputMesh() { return *i_BlockMesh; }
  const zObjMesh &getInputMesh() const { return *i_BlockMesh; }

  zObjMesh &getUMesh() { return _UMesh; }
  const zObjMesh &getUMesh() const { return _UMesh; }

  zObjGraphArray &getVLoops() { return _VLoops; }
  const zObjGraphArray &getVLoops() const { return _VLoops; }

  std::vector<zScalarArray> &getVScalar() { return _VScalar; }
  const std::vector<zScalarArray> &getVScalar() const { return _VScalar; }

  zObjGraphArray &getIsoMeshBoundary() { return o_IsoMeshBoundary; }
  const zObjGraphArray &getIsoMeshBoundary() const { return o_IsoMeshBoundary; }

  zObjGraphArray &getIsoMeshBoundary2D() { return o_IsoMeshBoundary2D; }
  const zObjGraphArray &getIsoMeshBoundary2D() const { return o_IsoMeshBoundary2D; }
  zObjGraphArray &getBracingGraphs2D() { return o_BracingGraphs2D; }
  const zObjGraphArray &getBracingGraphs2D() const { return o_BracingGraphs2D; }


  zObjMeshArray &getIsoMeshes() { return o_IsoMeshes; }
  const zObjMeshArray &getIsoMeshes() const { return o_IsoMeshes; }
  zObjMeshArray &getIsoMeshes2D() { return o_IsoMeshes2D; }
  const zObjMeshArray &getIsoMeshes2D() const { return o_IsoMeshes2D; }

  zObjMeshScalarField &getScalarField() { return o_ScalarField; }
  const zObjMeshScalarField &getScalarField() const { return o_ScalarField; }
  std::vector<zScalarArray> &getPathFields() { return o_PathFields; }
  const std::vector<zScalarArray> &getPathFields() const {
    return o_PathFields;
  }
  std::vector<zScalarArray> &getWidthFields() { return o_WidthFields; }
  const std::vector<zScalarArray> &getWidthFields() const {
    return o_WidthFields;
  }
  std::vector<zScalarArray> &getTrackFields() { return o_TrackFields; }
  const std::vector<zScalarArray> &getTrackFields() const {
    return o_TrackFields;
  }
  const zFloatArray &getVertexLayerHeights() const {
    return _VLength;
  }
  const zDomain<float> &getLayerHeightDomain() const {
    return o_LayerHeightDomain;
  }

  const ParameterGraph &getBracingParams() const { return _BracingParams; }

  /*! \brief Width group of every edge of layer `layer`'s boundary loop
   * (o_IsoMeshBoundary[layer], edge j = vertex j -> j + 1): the group of the
   * wall face between the two V-loops the edge spans, at that layer's
   * height. Requires computeSlice(). */
  void _layerEdgeWidthIds(int layer, int totalLayers, zIntArray &outEdgeIds) const;
  void _layerEdgeFaceTags(int layer, int totalLayers, zIntArray &outEdgeIds) const;

  //--------------------------
  //---- MAIN COMPUTE METHODS
  //--------------------------

  /*! \brief From the input mesh and its vertex tags: U-mesh (Bottom faces),
   * U-loop, V-loops with per-vertex values,
   * bracing parameterization - see the class comment. Fails (with a status
   * naming the offending vertex/face) on non-quad faces, missing tags, a
   * disconnected or multi-loop floor, or a column that does not reach Top
   * through regular (valence-4) vertices. */
  bool computeBasicStructure(bool useParallel = true) override;

  /*! \brief Computes layer section meshes and boundary/bracing graphs. */
  bool computeSlice(bool useParallel = true);

  /*! \brief 2D unrolling of section meshes and their boundary/bracing graphs.
   */
  bool computeUnroll(bool useParallel = true);

  /*! \brief Computes the 2D SDF path/width/track fields. */
  bool computeField(bool useParallel = true);

  /*! \brief Full pipeline: computeBasicStructure() -> computeSlice() ->
   * computeUnroll() -> computeField(). */
  bool computeSlicer(bool useParallel = true) override;

  // Extract each layer in parallel and map directly to its 3D section.
  bool computePrintPath(bool useParallel = true) override;

  /*! \brief Extracts paths and runs the remaining PrintSyn stages. */
  bool computePostProcess(const PrintSynSetting &settings, bool useParallel = true) override;

  const PrintPathArray &getPrintPaths2D() const override { return o_PrintPaths2D; }
  const PrintPathArray &getPrintPaths3D() const override { return o_PrintPaths3D; }
  const PrintPath &getResultantPrintPath() const override { return o_ResultantPrintPath; }
private:
  // Implementation helpers own detailed diagnostics.
  bool computeFieldImpl(bool useParallel);
  bool computePrintPathImpl(bool useParallel);
  bool computeSlicerImpl(bool useParallel);
  bool computePostProcessImpl(const PrintSynSetting &settings, bool useParallel);
  bool computeBasicStructureImpl(bool useParallel);
  bool computeSliceImpl(bool useParallel);
  bool computeUnrollImpl(bool useParallel);

};

/** @}*/

/** @}*/

} // namespace zSpace


#endif // ZSPACE_TS_3DP_TOPOLOGY_SLICER_H
