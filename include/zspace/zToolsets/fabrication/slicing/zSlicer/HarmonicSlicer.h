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

#ifndef ZSPACE_TS_GEOMETRY_HARMONICSLICER_H
#define ZSPACE_TS_GEOMETRY_HARMONICSLICER_H

#pragma once

#ifndef EIGEN_DONT_PARALLELIZE
#define EIGEN_DONT_PARALLELIZE
#endif

#include <zspace/zToolsets/zToolsetsExport.h>
#include <zspace/zToolsets/fabrication/slicing/utils/Utils.h>
#include <zspace/zToolsets/fabrication/slicing/zSlicer/zSlicer.h>

#include <zspace/interface.h>
#include <zspace/interface.h>
#include <zspace/zInterface/functionsets/zFnGraph.h>
#include <zspace/zInterface/functionsets/zFnMesh.h>
#include <zspace/zInterface/functionsets/zFnMeshField.h>
#include <zspace/zInterface/objects/zObjGraph.h>
#include <zspace/zInterface/objects/zObjMesh.h>
#include <zspace/zInterface/objects/zObjMeshField.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

#include <array>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace zSpace {

struct TriangleSoup {
  zPointArray Positions;
  std::vector<std::array<int, 3>> Triangles;
};

struct SurfaceOperators {
  Eigen::SparseMatrix<double> Laplacian;
  Eigen::VectorXd Mass;
  double MeanEdgeLength = 0.0;
};



// IDs are zero-based array indices; solver arrays remain contiguous.
class zTetrahedron {
public:
  struct TetVertex {
    zPoint Position;
    zIntArray EdgeIds, FaceIds, VolumeIds;
  };
  struct TetEdge {
    std::array<int, 2> VertexIds = {{-1, -1}};
    zIntArray FaceIds, VolumeIds;
  };
  struct TetFace {
    bool Valid = false;
    double Area = 0.0;
    std::array<int, 3> VertexIds = {{-1, -1, -1}};
    std::array<int, 3> EdgeIds = {{-1, -1, -1}};
    std::array<int, 2> VolumeIds = {{-1, -1}};
  };
  struct TetVolume {
    bool Valid = false;
    double Volume = 0.0;
    std::array<int, 4> VertexIds = {{-1, -1, -1, -1}};
    std::array<int, 6> EdgeIds = {{-1, -1, -1, -1, -1, -1}};
    std::array<int, 4> FaceIds = {{-1, -1, -1, -1}};
  };

  std::vector<TetVertex> Vertices;
  std::vector<TetEdge> Edges;
  std::vector<TetFace> Faces;
  std::vector<TetVolume> Volumes;
  
  zIntArray SurfaceFaceIds;
  const TetVertex &vertex(size_t id) const { return Vertices.at(id); }
  const TetVolume &volume(size_t id) const { return Volumes.at(id); }
  const std::array<int, 4> &volumeVertices(size_t id) const { return Volumes.at(id).VertexIds; }
};

class HarmonicField : public zTetrahedron {
public:
  // Surface cache for tetrahedralization and harmonic-field operations.
  TriangleSoup Surface;
  enum class TetFaceTag { Interior, Top, Bottom, Side };
  // Each array is indexed by its corresponding geometric element.
  zIntArray vertextag;
  zPointArray UVW;
  std::vector<bool> UVWValid;
  std::vector<TetFaceTag> FaceTags;
  zIntArray FaceWidthIds, FaceRegionTags;
  std::vector<Eigen::Matrix3d> FaceBasisGradients;
  std::vector<Eigen::Matrix<double, 3, 4>> VolumeBasisGradients;
  std::vector<Eigen::Matrix4d> VolumeStiffness;
  std::vector<Eigen::Vector3d> FaceFieldGradients, VolumeFieldGradients;
  // Temporary facet-marker data consumed when building topology.
  zIntArray SurfaceWidthIds, SurfaceRegionTags;
  std::vector<TetFaceTag> SurfaceTags;
  std::pair<zFloatArray, zFloatArray> HeatValue;
  bool FieldValid = false;
};

// Keys: sorted face vertex ids / zTsUtils::makeEdgeKey().
using FaceAdjacency = std::unordered_map<std::array<int, 3>, std::array<int, 2>, zTsUtils::TriangleKeyHash>;

using SurfaceEdgeAdjacency = std::unordered_map<unsigned long long, std::array<int, 2>>;

using EdgeWidthIds = std::unordered_map<unsigned long long, int>;

class ZSPACE_TOOLSETS HarmonicSlicer : public zSlicer {
protected:
  zObjMesh *i_BlockMesh;
  zIntArray i_FaceTags;
  zIntArray i_FaceWidthIds;
  zIntArray _VertexTags;
  SlicerSettings i_Settings;
  zObjGraph *i_InputBracing = nullptr;
  const zPoint *i_SlotPoint = nullptr;


  HarmonicField _Tet;

  zPointArray _BottomUV;
  zPoint _BottomUVScale = zPoint(1.0, 1.0, 0.0);
  zPoint _BottomUVOrigin = zPoint(0.0, 0.0, 0.0);
  zObjGraphArray o_TraceGraphs;

  std::string _Status;

  zObjMeshArray o_IsoMeshes;
  zObjGraphArray o_IsoMeshBoundary;
  std::vector<zPointArray> o_IsoMeshUV;
  std::vector<zFloatArray> o_VertexLayerHeights;
  std::vector<EdgeWidthIds> _HarmonicIsoEdgeWidthIds;

  zObjMeshArray o_IsoMeshes2D;
  zObjGraphArray o_IsoMeshBoundary2D;
  std::vector<zIntArray> o_IsoMeshBoundaryWidthId, o_IsoMeshBoundaryFaceTags;

  zObjMeshScalarField o_ScalarField;
  std::vector<zScalarArray> o_PathFields;
  std::vector<zScalarArray> o_WidthFields;
  std::vector<zScalarArray> o_TrackFields;
  std::vector<zTsField::FaceTagSource> o_FaceTagSources;
  zObjGraphArray o_BracingGraphs2D;

  PrintPathArray o_PrintPaths2D;
  PrintPathArray o_PrintPaths3D;
  PrintPath o_ResultantPrintPath;

public:
  HarmonicSlicer() {
    i_BlockMesh = nullptr;
    _Status = "Ready";
  }
  explicit HarmonicSlicer(const SlicerSettings &inSettings) : HarmonicSlicer() {
    i_Settings = inSettings;
  }
  ~HarmonicSlicer() {}

  void setInputMesh(zObjMesh &inMesh, const zIntArray &inFaceTags,
                    const zIntArray &inFaceWidthIds = zIntArray()) override {
    i_BlockMesh = &inMesh;
    i_FaceTags = inFaceTags;
    i_FaceWidthIds = inFaceWidthIds;
    clearStructure();
  }
  const std::vector<zTsField::FaceTagSource> &getFaceTagSources() const { return o_FaceTagSources; }
  const zIntArray &getFaceTags() const { return i_FaceTags; }
  const zIntArray &getFaceWidthIds() const { return i_FaceWidthIds; }
  const zIntArray &getVertexTags() const { return _VertexTags; }
  void setSettings(const SlicerSettings &inSettings) override {
    i_Settings = inSettings;
    clearStructure();
  }
  void setBracingGraph(zObjGraph &inGraph) override { i_InputBracing = &inGraph; clearFields(); }
  void clearBracingGraph() { i_InputBracing = nullptr; clearFields(); }
  void setSlotAnchor(const zPoint &inAnchor) override { i_SlotPoint = &inAnchor; clearFields(); }
  void clearSlotAnchor() override { i_SlotPoint = nullptr; clearFields(); }

  SlicerSettings &getSettings() override { return i_Settings; }
  const SlicerSettings &getSettings() const override { return i_Settings; }
  const std::string &getStatus() const override { return _Status; }

  bool isLoaded() const override { return i_BlockMesh != nullptr; }

  const std::pair<zFloatArray, zFloatArray> &getHeatValue() const { return _Tet.HeatValue; }

  const HarmonicField &getTetVolume() const { return _Tet; }
  zFloatArray getHarmonicField() const {
    if (!_Tet.FieldValid) return zFloatArray();
    zFloatArray field(_Tet.UVW.size());
    for (size_t i = 0; i < _Tet.UVW.size(); ++i) field[i] = _Tet.UVW[i].z;
    return field;
  }
  zObjMeshArray &getIsoMeshes() { return o_IsoMeshes; }
  const zObjMeshArray &getIsoMeshes() const { return o_IsoMeshes; }
  const std::vector<zPointArray> &getIsoMeshUV() const { return o_IsoMeshUV; }
  // One height per iso-mesh vertex, in model units: volume-weighted local
  // delta-w / |grad(w)|, interpolated at intersections (also on both caps).
  const std::vector<zFloatArray> &getVertexLayerHeights() const { return o_VertexLayerHeights; }
  const zObjGraphArray &getIsoMeshBoundary() const { return o_IsoMeshBoundary; }

  const zObjMeshArray &getIsoMeshes2D() const { return o_IsoMeshes2D; }
  const zObjGraphArray &getIsoMeshBoundary2D() const { return o_IsoMeshBoundary2D; }
  const std::vector<zIntArray> &getIsoMeshBoundaryWidthId() const { return o_IsoMeshBoundaryWidthId; }

  const zObjMeshScalarField &getScalarField() const { return o_ScalarField; }
  const std::vector<zScalarArray> &getPathFields() const { return o_PathFields; }
  const std::vector<zScalarArray> &getWidthFields() const { return o_WidthFields; }
  const std::vector<zScalarArray> &getTrackFields() const { return o_TrackFields; }
  const zObjGraphArray &getBracingGraphs2D() const { return o_BracingGraphs2D; }

  const zPointArray &getVolumetricParameterization() const { return _Tet.UVW; }
  const std::vector<bool> &getVolumetricParameterizationValidity() const { return _Tet.UVWValid; }
  zObjGraphArray &getTraceGraphs() { return o_TraceGraphs; }
  const zObjGraphArray &getTraceGraphs() const { return o_TraceGraphs; }

  bool computeBasicStructure(bool useParallel = true) override;
  bool computeHeatDistances();
  bool computeHarmonicField();
  bool computeVolumetricParameterization(bool useParallel = true);
  // Generates matched 3D/2D iso meshes and boundaries; parallelizes by layer.
  bool computeContours(bool useParallel = true);
  bool computeField(bool useParallel = true);
  bool computeSlicer(bool useParallel = true) override;
  // Extract each layer in parallel and map directly to its 3D section.
  bool computePrintPath(bool useParallel = true) override;

  bool computePostProcess(const PrintSynSetting &settings, bool useParallel = true) override;

  const PrintPathArray &getPrintPaths2D() const override { return o_PrintPaths2D; }
  const PrintPathArray &getPrintPaths3D() const override { return o_PrintPaths3D; }
  const PrintPath &getResultantPrintPath() const override { return o_ResultantPrintPath; }

private:
  bool columnLengths(double &outShortest, double &outLongest) const;

  void clearHeat();
  void clearPrintPaths();
  void clearFields();
  void clearUnroll();
  void clearContours();
  // Clears (u,v) while retaining the harmonic w field for a tracing rerun.
  void clearParameterization();
  void clearStructure();


  void buildSurfaceSoup();
  bool buildTetVolume();
  bool buildTetTopology();
  void updateFieldGradients();
  void buildTetElements(int &degenerateTets, bool useParallel);

  static bool buildSurfaceOperators(const TriangleSoup &soup,
                                    SurfaceOperators &ops,
                                    int &degenerateTriangles);
  static bool computeHeatGeodesics(const TriangleSoup &soup,
                                   const SurfaceOperators &ops,
                                   Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> &heatSolver,
                                   const zIntArray &sourceIds,
                                   zFloatArray &outDistance,
                                   std::string &error);

  struct TetDirichletSystem {
    std::vector<int> FreeIndex;
    zFloatArray FixedField;
    Eigen::VectorXd BoundaryRhs;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> Solver;
  };
  bool prepareTetDirichlet(const std::vector<double> &fixedValue,
                           TetDirichletSystem &system, std::string &error) const;
  bool solveTetDirichlet(TetDirichletSystem &system, zFloatArray &outField,
                         std::string &error,
                         const std::vector<Eigen::Vector3d> *targetGradients = nullptr) const;
  std::vector<double> capDirichletValues() const;
  bool equaliseField(TetDirichletSystem &system, zFloatArray &ioField,
                      std::string &outReport) const;
  Eigen::Vector3d tetGradient(size_t tet, const zFloatArray &field) const;

  bool buildBottomParameterization();
  struct VolumeTraceData {
    std::vector<int> Seeds;
    FaceAdjacency Adjacency;
    std::vector<Eigen::Vector3d> Gradient;
  };
  struct WallTraceData {
    std::vector<int> Seeds;
    SurfaceEdgeAdjacency Adjacency;
    std::vector<Eigen::Vector3d> Gradient;
    std::vector<std::vector<int>> IncidentTriangles;
  };
  bool traceToBottom(int vertexId, const VolumeTraceData &data,
                     zPoint &outUV, zPointArray *outPath = nullptr,
                     std::string *failureReason = nullptr) const;
  static bool computeTriElement(const std::vector<zTetrahedron::TetVertex> &vertices,
                                const std::array<int, 3> &tri,
                                zTetrahedron::TetFace &element, Eigen::Matrix3d &gradients);
  bool traceWallToBottom(int vertexId, const WallTraceData &data,
                         zPoint &outUV, zPointArray *outPath = nullptr,
                         std::string *failureReason = nullptr) const;

  struct ContourData;
  bool computeVertexLayerHeights(int numLayers, zFloatArray &heights) const;
  int crossingVertex(int a, int b, float iso, const zFloatArray &vertexHeights, ContourData &slice);
  static bool
  appendTriangle(zIntArray &connects, int a, int b, int c,
                 const zPointArray &positions, const zVector &lowToHigh,
                 std::unordered_set<std::array<int, 3>, zTsUtils::TriangleKeyHash>
                     &uniqueTriangles);
  struct WallWidthLookup;
  void extractIsoMesh(float iso, zObjMesh &mesh, zPointArray &uv,
                      const zFloatArray &vertexHeights, zFloatArray &heights,
                      const WallWidthLookup &lookup, EdgeWidthIds &edgeIds, std::vector<std::string> &debug);
  void extractCapMesh(VertexTag tag, zObjMesh &mesh, zPointArray &uv,
                      const zFloatArray &vertexHeights, zFloatArray &heights,
                      const WallWidthLookup &lookup, EdgeWidthIds &edgeIds, std::vector<std::string> &debug) const;
  struct WallWidthLookup {
    std::unordered_map<std::array<int, 3>, int, zTsUtils::TriangleKeyHash> Triangles;
    std::unordered_map<unsigned long long, int> Edges;
  };
  void extractWidthId(const zIntArray &connects, const zIntArray &counts,
                              const std::vector<std::pair<int, int>> &vertexTetEdges,
                              const WallWidthLookup &lookup, EdgeWidthIds &edgeIds, std::vector<std::string> &debug) const;
  bool unrollIsoMesh(zObjMesh &isoMesh, const zPointArray &uv, zObjMesh &mesh2D,
                     zPointArray &positions3D, zPointArray &positions2D) const;
  static void createIsoBoundaries(const zPointArray &positions3D,
                                   const zPointArray &positions2D, const zIntArray &connects,
                                   const zIntArray &counts,
                                   const EdgeWidthIds &widthIds, zObjGraph &boundary3D,
                                   zObjGraph &boundary2D, zIntArray &boundaryWidthIds, std::vector<std::string> &debug);
private:
  // Implementation helpers own detailed diagnostics; public compute methods report results.
  bool computeFieldImpl(bool useParallel);
  bool computePrintPathImpl(bool useParallel);
  bool computeHarmonicFieldImpl();
  bool computeHeatDistancesImpl();
  bool computeSlicerImpl(bool useParallel);
  bool computePostProcessImpl(const PrintSynSetting &settings, bool useParallel);
  bool computeBasicStructureImpl(bool useParallel);
  bool computeVolumetricParameterizationImpl(bool useParallel);

};

} // namespace zSpace


#endif // ZSPACE_TS_GEOMETRY_HARMONICSLICER_H
