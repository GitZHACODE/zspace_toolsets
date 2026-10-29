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

#ifndef ZSPACE_TS_3DP_UTILS_H
#define ZSPACE_TS_3DP_UTILS_H

#pragma once

// The one shared header of the ts3dp toolset: common types, geometry/topology
// and coloring utilities (zTsUtils), 2D SDF field kernels (zTsField) and USD
// file IO (zTsIO). Implementations stay split over utils/*.cpp.

#ifndef EIGEN_DONT_PARALLELIZE
#define EIGEN_DONT_PARALLELIZE
#endif

#include <zspace/zToolsets/zToolsetsExport.h>
#include <zspace/zToolsets/fabrication/slicing/utils/SlicingTiming.h>
#include <memory>

#include <zspace/interface.h>
#include <zspace/io.h>
#include <zspace/zInterface/functionsets/zFnGraph.h>
#include <zspace/zInterface/functionsets/zFnMesh.h>
#include <zspace/zInterface/functionsets/zFnMeshField.h>
#include <zspace/zInterface/objects/zObjGraph.h>
#include <zspace/zInterface/objects/zObjMesh.h>
#include <zspace/zInterface/objects/zObjMeshField.h>

#include <Eigen/Core>

#include <array>
#include <string>
#include <vector>

namespace zSpace {

/** \addtogroup zToolsets
 *	\brief Collection of toolsets for geometry, analysis, and fabrication.
 *  @{
 */

/** \addtogroup zTsSlicing
 *	\brief 3D printing toolset: slicer operators, print-path post-processing.
 *  @{
 */

//--------------------------
//---- COMMON CONTAINER TYPEDEFS
//--------------------------

typedef std::vector<zObjMesh> zObjMeshArray;
typedef std::vector<zObjGraph> zObjGraphArray;
typedef std::vector<zObjGraph *> zObjGraphPointerArray;
// zObjMeshScalarFieldArray intentionally not typedef'd here - collides with
// zObjMeshField.h's own definition. Include zObjMeshField.h directly instead.

//--------------------------
//---- COMMON DATA STRUCTURES
//--------------------------

/*! \enum VertexTag
 *	\brief Tags a boundary loop / capped-mesh vertex as bottom, top, side
 *  wall, or unconstrained interior (tet-volume marker convention).
 */
enum VertexTag : int { Interior = 0, Bottom = 1, Top = 2, Side = 3 };
inline bool isSideFaceTag(int tag) { return tag >= 3 && tag <= 6; }
// Input convention: one tag per input-mesh FACE (Bottom = print floor, Top =
// print ceiling, Side = wall), from USD FaceTag attributes
// (zTsIO::readSlicerInputUSD()). The slicers derive the per-vertex tags they
// work with from them - zTsUtils::vertexTagsFromFaceTags().

// Independent fabrication bounds, never inferred from measured output stats.
struct PrintHeightOptimizationSettings {
  bool Enabled = false; // Harmonic only; bounds use the mesh's model units.
  double Tolerance = 0.0; // Fabrication allowance in model units; measured values are never clamped.
  double Minimum = 6.0, Maximum = 12.0;
  int Iterations = 16; // Total outer correction budget across seed, feedback and count candidates.
  bool RedistributeIsoValues = true;
  bool AllowLayerCountChange = true;
  bool OptimizeSurfaces = true; // False retains scheduling only; never enables the old path-target Newton solver.
  double MaxDisplacement = 0.0; // 0: two maximum target heights, from the relaxed seed.
};

/*! \struct SlicerSettings
 *	\brief One settings set for both slicer operators. Each operator reads
 *  only what it needs: shared plain names, T* for Topology, H* for Harmonic.
 *  PrintSynSetting lives in zPrintSyn.h.
 */
struct SlicerSettings {
  //---- Shared: SDF print-field inputs (zTsField)
  // Outer-wall bead width per width group: a wall face tagged side_<n>
  // ("side" = group 0) prints BoundaryWidth[n]; groups past the end use the
  // last entry. Both slicers split each boundary loop into runs of the same
  // group (zTsField::loopWidthSegments()) - see zTsSlicing::getWidthGroupCount().
  std::vector<float> BoundaryWidth = {0.5f};
  float BracingWidth = 0.5f;
  // Optional per-component widths, ordered by smallest input vertex index.
  zFloatArray BracingCurveWidths;
  float OverlapWidth = 0.0f;
  bool DoubleBoundary = false;
  // Boundary trim slot, cut only with DoubleBoundary (it joins the inner and
  // outer loop into one chain; a single wall is already one closed loop).
  // Placed from the slot anchor point (zSlicer::setSlotAnchor(),
  // world coordinates): the anchor's nearest point on the floor boundary,
  // then SlotOffset * w along the loop - even layers counter-
  // clockwise, odd layers clockwise (seen from +Z), w the boundary width
  // there. Both slicers place it the same way.
  float SlotOffset = 1.0f;
  // Base bracing position (0..1). Even/odd layers offset by +/- SlotOffset
  // times the rib width in arc length; open curves clamp, closed curves wrap.
  float BracingSlot = 0.5f;
  float FieldScaleUp = 1.05f;
  // Grid cell size, used exactly as given - keep it comfortably below the
  // narrowest bead width or the zero-isocontour may not resolve.
  float FieldCellSize = 0.05f;
  // Layer spacing - the only layering input: ceil(L / LayerSpacing) + 1
  // layers (floor and ceiling included, zTsUtils::layersFromSpacing()), L
  // the shortest Bottom->Top distance. TopologySlicer: the shortest path
  // along the input mesh's edges (zTsUtils::shortestBottomToTopLength());
  // HarmonicSlicer: the shortest heat-method geodesic on the wall.
  float LayerSpacing = 1.0f;

  //---- Topology only (TopologySlicer): none - the boundary segments come
  // from the face width groups.

  //---- Harmonic only (HarmonicSlicer)
  // Tetrahedral target edge length relative to the input surface's mean
  // triangle edge. Must be finite and in (0,1]. 1 restores the original
  // density; 0.5 refines the volume without modifying the input mesh.
  double _HTetEdgeLengthScale = 1.0;
  // Constrained relaxation of intermediate harmonic layers before unrolling.
  // Caps remain fixed; boundary vertices slide on their input side patches.
  // 0 iterations retains the original geometry, height and face-normal path.
  int _HIsoRelaxIterations = 8;
  double _HIsoRelaxStrength = 0.25; // finite (0,1]
  // Model units; 0 uses a limit derived from neighbouring-layer separation.
  double _HIsoRelaxMaxDeviation = 0.0;
  PrintHeightOptimizationSettings HPrintHeight;
  // Maximum Euler tracing step as a fraction of the current simplex's
  // minimum altitude (volume and wall traces). Must be in (0, 1].
  double _HTraceStepFraction = 0.2;
  int _HTraceMaxSteps = 100000; // includes interior steps and boundary transitions
  // Gradient equalisation of the slicing field: after the harmonic solve,
  // this many least-squares re-solves pull |grad| toward its volume-weighted
  // mean along the current flow direction (same Dirichlet caps), so the
  // evenly spaced iso values give more even layer heights. 0 = the plain
  // harmonic field. The layer count does not depend on it.
  int _HEqualiseIterations = 2;
  // ARAP (as-rigid-as-possible) iterations unrolling every slice at world
  // scale (HarmonicSlicer::computeContours()). 0 = keep the start layout.
  int _HArapIterations = 10;
  // ARAP start layout of every slice, which the result is also rigidly
  // aligned back to: 0 = planar projection of the slice itself
  // (zTsUnroll::planarProjection(), its top view - default), 1 = its (u,v)
  // in the floor's physical frame.
  int _HUnrollStart = 0;
};

/*! \struct ParameterGraph
 *	\brief Graph parameterized by 2D/3D mesh face barycentric coordinates.
 */
struct ParameterGraph {
  struct Vertex {
    int v0, v1, v2; // Triangle 3 vertex indices
    float u, v, w;  // Barycentric coordinate weights
  };

  std::vector<Vertex> vertices;
  zIntArray edgeConnects;

  bool empty() const { return vertices.empty() || edgeConnects.empty(); }
  void clear() {
    vertices.clear();
    edgeConnects.clear();
  }
};

/*! \struct PrintPath
 *	\brief Shared toolpath exchange type between the slicer operators and
 * zPrintSyn.
 */
struct PrintPath {
  // Legacy geometric bits remain available to standalone reconstruction users.
  // Production features are outer FaceTag junctions; no storage-layout change.
  enum Feature { AngleFeature = 1, IntersectionFeature = 2, OuterJunctionFeature = 4 };
  zObjGraph graph;
  zFloatArray widths;
  zFloatArray heights;
  zVectorArray normals;
  // Per-vertex inner(0)/outer(1) tag, parallel to widths. Populated by
  // zPrintSyn::extractTrackPaths for a double-layer boundary ring. A single
  // perimeter is outer (1); an empty array means no track metadata is supplied.
  zFloatArray trackIds;
  zIntArray faceTags, edgeFaceTags;
  struct FaceTagSegment { int tag=0; zIntArray vertices; };
  std::vector<FaceTagSegment> faceTagSegments;
  // Optional per-vertex bit mask, parallel to graph vertices when populated.
  // Production stages refresh OuterJunctionFeature from incident side ownership;
  // mapping carries aligned metadata. Legacy angle bits do not select features.
  zIntArray featureFlags;
};

// Print paths per layer: [layer][path], as produced by the zPrintSyn
// functions and stored by each slicer operator's computePostProcess().
typedef std::vector<std::vector<PrintPath>> PrintPathArray;

//--------------------------
//---- COMMON GEOMETRY / TOPOLOGY UTILITIES        (utils/Utils.cpp)
//---- FIELD-MESH COLORING UTILITIES               (utils/Color.cpp)
//--------------------------

namespace zTsUtils {

// Loft one 16-sided ellipse per print vertex: full width w, full height h.
// The TCP is its top; its center is -h/2 and bottom is -h along the height axis.
// Sections are normal to the neighbor chord, with the print normal projected
// into the section's height axis. Chords keep a consistent forward direction.
// Each graph component/layer is independent; open paths receive end caps.
// Zero dimensions are retained (collapsed sections); no artificial thickness.
// Failure clears outMesh and reports the layer/path through outError.
ZSPACE_TOOLSETS bool buildPrintMesh(const PrintPathArray &layers, zObjMesh &outMesh,
                                   std::string *outError = nullptr);

//---- Thread-Safe Console Logging

ZSPACE_TOOLSETS void printSuccess(const std::string &msg);
ZSPACE_TOOLSETS void printDebug(const std::string &msg);
ZSPACE_TOOLSETS void printError(const std::string &msg);

//---- Basic Vector Utilities

ZSPACE_TOOLSETS zVector vecSub(const zPoint &a, const zPoint &b);
ZSPACE_TOOLSETS double vecDistSq(const zPoint &a, const zPoint &b);
ZSPACE_TOOLSETS double vecDot(const zVector &a, const zVector &b);
ZSPACE_TOOLSETS zVector vecCross(const zVector &a, const zVector &b);
ZSPACE_TOOLSETS unsigned long long makeEdgeKey(int u, int v);
ZSPACE_TOOLSETS Eigen::Vector3d toEigenVec(const zPoint &p);

// Twice the area of triangle abc.
ZSPACE_TOOLSETS double triangleArea2(const zPoint &a, const zPoint &b, const zPoint &c);

// Mean edge length over the non-degenerate triangles (0 if there are none).
ZSPACE_TOOLSETS double meanEdgeLength(const zPointArray &positions, const std::vector<std::array<int, 3>> &triangles);

// Hash of a vertex-id triple (e.g. a sorted triangle or tet face).
struct TriangleKeyHash {
  size_t operator()(const std::array<int, 3> &key) const noexcept {
    size_t seed = 0;
    for (int v : key) seed = seed * 1000003u + static_cast<size_t>(static_cast<unsigned int>(v));
    return seed;
  }
};

// Angle (radians) at cur between the edges to prev and next. Returns false
// (outAngleRad untouched) if either edge is degenerate (< 1e-3 mm long).
ZSPACE_TOOLSETS bool turnAngleAt(const zPoint &prev, const zPoint &cur, const zPoint &next,
                 double &outAngleRad);

// Signed area (shoelace formula) of a closed 2D polyline in the XY plane;
// positive = counter-clockwise, negative = clockwise.
ZSPACE_TOOLSETS double signedArea2D(const zPointArray &positions);

/*! \brief Ordered boundary loops of a polygon mesh given as counts/connects
 * (vertex ids into positions), each following its faces' own winding.
 * outLoop gets the longest one (vertex ids); returns the number of loops. */
ZSPACE_TOOLSETS int outerBoundaryLoop(const zIntArray &connects, const zIntArray &counts,
                      const zPointArray &positions, zIntArray &outLoop);

//---- HalfEdge Mesh Topology Walk Utilities

ZSPACE_TOOLSETS zItMeshHalfEdge nextEdgeRing(zItMeshHalfEdge he);
ZSPACE_TOOLSETS zItMeshHalfEdge nextEdgeLoop(zItMeshHalfEdge he);
ZSPACE_TOOLSETS zItMeshHalfEdge zUpEdge(zItMeshVertex v);
ZSPACE_TOOLSETS void getFaceVerticesFromHalfedge(zItMeshHalfEdge &heStart, bool forward,
                                 zIntArray &faceVerts);

//---- Combine Multiple Graph Objects

ZSPACE_TOOLSETS void combineGraphObjects(const zObjGraphArray &sources, zObjGraph &outGraph);

//---- Single-Chain Graph Closed/Open Detection

/*! \brief Infers closed/open from a wraparound edge (vertex 0 <-> N-1). */
ZSPACE_TOOLSETS bool isClosedGraph(const zObjGraph &graph);

//---- Input Tags

/*! \brief Per-vertex VertexTag from one VertexTag per mesh face (in the
 * mesh's polygon order): a vertex of any Bottom face is Bottom, of any Top
 * face Top, every other vertex Side. False (outError says why) if the
 * counts differ, a face tag is not Bottom/Top/Side, there is no Bottom or no
 * Top face, or a vertex lies on both a Bottom and a Top face (no wall
 * between floor and ceiling there). */
ZSPACE_TOOLSETS bool vertexTagsFromFaceTags(zObjMesh &mesh, const zIntArray &faceTags,
                            zIntArray &outVertexTags,
                            std::string *outError = nullptr);

/*! \brief Number of width groups the side faces use: the largest width id
 * of a Side face + 1 (faceWidthIds empty = every face group 0), 0 without
 * side faces - the length SlicerSettings::BoundaryWidth needs. */
ZSPACE_TOOLSETS int widthGroupCount(const zIntArray &faceTags, const zIntArray &faceWidthIds);

//---- Layering Shared by the Slicer Operators

/*! \brief Length of the shortest path along the given edges (pairs of
 * vertex ids) from any Bottom-tagged vertex to any Top-tagged vertex
 * (multi-source Dijkstra). -1 if the tags don't match the vertex count or no
 * such path exists. */
ZSPACE_TOOLSETS double shortestBottomToTopLength(const zPointArray &positions,
                                 const zIntArray &edgePairs,
                                 const zIntArray &vertexTags);

/*! \brief Same, over the mesh's own (polygon) edges. */
ZSPACE_TOOLSETS double shortestBottomToTopLength(zObjMesh &mesh, const zIntArray &vertexTags);

/*! \brief Number of layers both slicer operators cut for a shortest
 * Bottom->Top length and a LayerSpacing: ceil(length / spacing) + 1 (floor
 * and ceiling included, >= 2). 0 if length or spacing is not positive. */
ZSPACE_TOOLSETS int layersFromSpacing(double length, double spacing);

//---- Universal Polyline Sampling by Arc-Length with Optional Tangent Vector

/*! \brief Orders one connected, unbranched graph into polyline points and
 * cumulative 3D arc lengths (first length = 0). A closed loop repeats its
 * first point at the end and includes the closing edge's length.
 * Starts at the lowest-index endpoint, or vertex 0 for a loop; loop
 * direction follows the first incident edge in the input edge order.
 * False, with both outputs cleared, for empty/invalid data, isolated
 * vertices, multiple components, branching, duplicate/self edges,
 * zero-length edges or non-finite positions/lengths. Does not test geometric
 * self-intersections. Does not modify the source graph. */
ZSPACE_TOOLSETS bool graphtopolyline(const zObjGraph &graph, zPointArray &outPoints,
                     zScalarArray &outArcLengths);

struct GraphPolyline {
  zIntArray Vertices, Edges;
  bool Closed = false;
};
// Orders each connected chain/loop by real edges; rejects branching topology.
ZSPACE_TOOLSETS bool graphPolylines(const zObjGraph &graph, std::vector<GraphPolyline> &outCurves);

// minSegmentLength controls interpolation on very short arc segments;
// the default preserves existing callers, anchored XY slots use 1e-9 mm.
ZSPACE_TOOLSETS bool samplePolyline(const zPoint *positions, const zScalar *scalars,
                    size_t count, double t, zPoint &outPoint,
                    zVector *outTangent = nullptr,
                    double minSegmentLength = 1e-3);
ZSPACE_TOOLSETS bool samplePolyline(const zObjGraph &polylineGraph, const zScalarArray &scalars,
                    double t, zPoint &outPoint, zVector *outTangent = nullptr);

//---- Universal Mesh Barycentric Parameterizer (2D unrolled and 3D space meshes)

/*! \brief Locates every graph vertex on inMesh's triangles (nearest point,
 * 2D or 3D) as barycentric coordinates. meshPositions (optional): locate
 * against these per-vertex positions instead of inMesh's own (same vertex
 * order) - e.g. a per-vertex (u,v) attribute, to find where a (u,v) point
 * lies on the mesh. allowNearest=false rejects points outside the domain,
 * non-finite points and ambiguous overlapping chart triangles. */
ZSPACE_TOOLSETS bool parameterizeGraph(const zObjGraph &inGraph, zObjMesh &inMesh,
                       ParameterGraph &outParams,
                       const zPointArray *meshPositions = nullptr,
                       bool allowNearest = true, std::string *outError = nullptr,
                       const zPointArray *diagnosticTargetPositions = nullptr);

/*! \brief Interpolates a per-mesh-vertex value array at every parameterized
 * graph vertex (outValues 1:1 with inParams.vertices). False if a referenced
 * vertex is out of range. */
ZSPACE_TOOLSETS bool interpolateParameterGraph(const ParameterGraph &inParams,
                               const zPointArray &vertexValues,
                               zPointArray &outValues);

//---- Universal Barycentric Mesh Graph Reconstructor

ZSPACE_TOOLSETS bool realizeParameterGraph(const ParameterGraph &inParams, zObjMesh &targetMesh,
                           zObjGraph &outGraph,
                           zVectorArray *outNormals = nullptr,
                           const zFloatArray *vertexScalars = nullptr,
                           zFloatArray *outHeights = nullptr,
                           const zVectorArray *vertexNormals = nullptr);

//---- Barycentric Mesh-to-Mesh Graph Projection

ZSPACE_TOOLSETS bool remapGraph(zObjGraph &graph, zObjMesh &inMesh2D,
                zObjMesh &projectionMesh3D,
                zVectorArray *mappedNormals = nullptr,
                const zFloatArray *vertexHeights = nullptr,
                zFloatArray *mappedHeights = nullptr, std::string *outError = nullptr,
                const zVectorArray *vertexNormals = nullptr);

//---- Barycentric Point Transfer (shared by both slicer operators and zPrintSyn)

/*! \brief Transfers one point by barycentric coordinates: located once on
 * fromMesh (nearest point on its triangles, 2D or 3D), realized on every
 * toMeshes[i] (each with fromMesh's topology; null entries are skipped).
 * outValid[i] = 0 where nothing was realized. outDistance (optional): from
 * point to its located position on fromMesh, 0 on/inside it. False if point
 * cannot be located on fromMesh. */
ZSPACE_TOOLSETS bool mapPointToMeshes(const zPoint &point, zObjMesh &fromMesh,
                      const std::vector<const zObjMesh *> &toMeshes,
                      zPointArray &outPoints, std::vector<char> &outValid,
                      double *outDistance = nullptr, bool useParallel = true);

/*! \brief mapPointToMeshes() onto a single mesh. */
ZSPACE_TOOLSETS bool mapPoint(const zPoint &point, zObjMesh &fromMesh, zObjMesh &toMesh,
              zPoint &outPoint, double *outDistance = nullptr);

//---- Toolpath Pipe Mesh Construction

/*! \brief Builds a 3D extruded pipe mesh from a print path. */
ZSPACE_TOOLSETS void buildPathMesh(const PrintPath &path, float defaultLayerHeight,
                   zObjMesh &outMesh);

//---- Color Format Conversion

ZSPACE_TOOLSETS zColor hslToRgb(float h, float s, float l);

//---- Distinct Shuffled Color Palette Generator

ZSPACE_TOOLSETS void generateColorPalette(int numColors, std::vector<zColor> &outColors,
                          unsigned int seed = 1337, float saturation = 0.85f,
                          float lightness = 0.50f);

//---- Field-Mesh Coloring: Positive Values Solid, Negative Values a
//---- Gradient Between Two Input Colors (near-zero -> most-negative)

ZSPACE_TOOLSETS void colorPathField(
    const zScalarArray &values, std::vector<zColor> &outColors,
    const zColor &positiveColor = zColor(1.0f, 1.0f, 1.0f, 1.0f),
    const zColor &negativeColorNearZero = zColor(0.0f, 40.0f / 255.0f,
                                                 240.0f / 255.0f, 1.0f),
    const zColor &negativeColorAtMin = zColor(180.0f / 255.0f, 200.0f / 255.0f,
                                              1.0f, 1.0f));
ZSPACE_TOOLSETS void colorPathField(
    zObjMesh &mesh, const zScalarArray &values,
    const zColor &positiveColor = zColor(1.0f, 1.0f, 1.0f, 1.0f),
    const zColor &negativeColorNearZero = zColor(0.0f, 40.0f / 255.0f,
                                                 240.0f / 255.0f, 1.0f),
    const zColor &negativeColorAtMin = zColor(180.0f / 255.0f, 200.0f / 255.0f,
                                              1.0f, 1.0f));

//---- Field-Mesh Coloring: Discrete Lookup by Index Into a Generated Color
//---- List (one distinct color per unique value)

// paletteWidths (optional): build the palette from these widths instead of
// from values, so several fields / meshes colored with the same list get the
// same color for the same width even when each holds only some of them.
ZSPACE_TOOLSETS void colorWidthField(const zScalarArray &values,
                     std::vector<zColor> &outColors,
                     const zScalarArray *paletteWidths = nullptr);
ZSPACE_TOOLSETS void colorWidthField(zObjMesh &mesh, const zScalarArray &values,
                     const zScalarArray *paletteWidths = nullptr);

//---- Generic Mesh/Graph Transform & Segmentation

/*! \brief Translates all vertex positions of a mesh. */
ZSPACE_TOOLSETS void translateMesh(zObjMesh &mesh, const zVector &offset);

/*! \brief Translates all vertex positions of a graph. */
ZSPACE_TOOLSETS void translateGraph(zObjGraph &graph, const zVector &offset);

/*! \brief Segments a graph at specified split vertex indices. */
ZSPACE_TOOLSETS bool segmentGraph(const zObjGraph &sourceBoundaryGraph,
                  const std::vector<int> &splitIndices,
                  std::vector<zObjGraph> &outBoundarySegments);

} // namespace zTsUtils

//--------------------------
//---- SECTION UNROLLING             (utils/Unroll.cpp)
//--------------------------
// Both ways of unrolling a 3D section mesh; each slicer operator calls the
// one that suits its sections: TopologySlicer the face-by-face unfolding
// (one order shared by every layer, all layers having the U-mesh's
// topology), HarmonicSlicer the ARAP unrolling (every layer its own mesh).

namespace zTsUnroll {

/*! \brief Face-by-face unfolding order of a polygon mesh: a breadth-first
 * spanning tree of its dual graph from the highest-valence face, as one
 * (parent, child) shared-half-edge pair per tree edge. False for a mesh of
 * one face or a disconnected dual graph. */
ZSPACE_TOOLSETS bool bsfUnfoldOrder(zObjMesh &mesh, zIntPairArray &outHalfEdgePairs);

/*! \brief Unfolds sourceMesh face by face along bsfHalfEdgePairs (from
 * bsfUnfoldOrder() on this mesh or on one of the same topology): each face is
 * laid isometrically against its already placed parent edge, every vertex
 * is placed once. outUnrollMesh: a copy of sourceMesh with the unrolled (z=0)
 * positions. Exact on a developable mesh only. */
ZSPACE_TOOLSETS bool unfoldMesh(zObjMesh &sourceMesh, const zIntPairArray &bsfHalfEdgePairs,
                zObjMesh &outUnrollMesh);

/*! \struct ArapStats
 *	\brief Length distortion of an unrolling: over its triangles, the
 * largest and the area-weighted mean of max(|s1 - 1|, |s2 - 1|) (s1, s2 the
 * singular values of the triangle's 3D -> 2D map; 0 = lengths kept), and the
 * number of flipped triangles.
 */
struct ArapStats {
  double stretchMax = 0.0;
  double stretchMean = 0.0;
  int flipped = 0;
};

/*! \brief As-rigid-as-possible unrolling (Liu et al. 2008, local/global
 * iterations, cotangent weights, one rotation per triangle) of a 3D mesh at
 * world scale: the triangle shapes come from mesh3D (polygons are fanned),
 * start2D (one point per mesh vertex, z ignored) is only the starting layout
 * - e.g. planarProjection(), or a (u,v) chart - whose handedness the result
 * keeps. A tiny pull to start2D only fixes each component's free
 * translation. stats is measured on the result. False (out2D = start2D) for
 * bad input or a failed solve. */
ZSPACE_TOOLSETS bool arapUnroll(zObjMesh &mesh3D, const zPointArray &start2D, int iterations,
                 zPointArray &out2D, ArapStats &stats);

/*! \brief Default ARAP start: mesh3D's vertices projected onto its
 * best-fit plane, the plane turned onto XY by the smallest rotation taking
 * its normal (oriented to +Z) to +Z, about the centroid. A near-horizontal
 * mesh comes out as its plain top view (z dropped); a tilted one keeps its
 * in-plane direction. */
ZSPACE_TOOLSETS void planarProjection(zObjMesh &mesh3D, zPointArray &out2D);

/*! \brief Rotates and translates points (no scaling, no mirroring) onto
 * target in the least-squares sense, 1:1 by index (XY plane, z = 0). */
ZSPACE_TOOLSETS void alignRigid2D(zPointArray &points, const zPointArray &target);

} // namespace zTsUnroll

//--------------------------
//---- SHARED 2D SDF PRINT-FIELD KERNELS          (utils/Field.cpp)
//--------------------------
// Used by every slicer operator's computeField() (TopologySlicer per unrolled
// layer, HarmonicSlicer per unrolled layer) and by zPrintSyn's isocontour
// extraction. All geometry is read in the XY plane.

namespace zTsField {

/*! \struct FieldSegment
 *	\brief A boundary/Slicer_bracing graph with its print width, offset and trim
 * ratio - one input curve of the path/width fields. trim is read for
 * bracing only; the boundary slot is anchored (createAnchoredSlot()).
 */
struct FieldSegment {
  zObjGraph graph;
  int type = 0;         // 0: Outer Boundary, 1: Internal Bracing
  zScalar width = 0.0;  // Print bead width
  zScalar offset = 0.0; // Distance offset
  zScalar trim = 0.0;   // Trim slot ratio

  FieldSegment() = default;
  FieldSegment(const zObjGraph &inGraph, int inType, zScalar inWidth,
               zScalar inOffset, zScalar inTrim = 0.0)
      : graph(inGraph), type(inType), width(inWidth), offset(inOffset),
        trim(inTrim) {}
};

// Compact continuous field provenance, in the layer's flattened chart.
// Geometry fields remain grouped by width; these labels never alter the CSG.
struct FaceTagSource {
  zObjGraph boundary, slot;
  zIntArray edgeTags, edgeWidthIds;
  std::vector<float> widths;
  bool doubleBoundary = false;
  float overlap = 0, slotWidth = 0;
  std::vector<FieldSegment> bracing;
};

// Immutable snapshot of the source's prepared distance fields and original
// edge identities. Safe for concurrent reads; construct again after changing
// source geometry/settings. Copies share the same prepared storage.
class ZSPACE_TOOLSETS FaceTagQuery {
public:
  FaceTagQuery() = default;
  explicit FaceTagQuery(const FaceTagSource &source);
  bool valid() const;
  int label(const zPoint &point) const;
  bool query(const zPointArray &points, zIntArray &tags) const;
private:
  struct Impl;
  std::shared_ptr<const Impl> _impl;
};

// Classify using the same boundary/slot distance fields as path generation.
// Inserts exact category-transition anchors and records directed edge labels.
ZSPACE_TOOLSETS bool labelFaceTagSegments(PrintPath &path, const FaceTagSource &source, bool splitTransitions = true);
ZSPACE_TOOLSETS bool labelFaceTagSegments(PrintPath &path, const FaceTagQuery &query, bool splitTransitions = true);
ZSPACE_TOOLSETS bool queryFaceTags(const FaceTagSource &source, const zPointArray &points, zIntArray &tags);

/*! \brief BoundaryWidth entry of width group id: the last entry for ids
 * past the end; empty tables and widths below 1e-4 mm are invalid. */
// Minimum accepted physical parameter in internal millimetres (not a geometry tolerance).
inline constexpr float MinPhysicalParameter = 1e-4f;
ZSPACE_TOOLSETS bool validateFieldSettings(const SlicerSettings &settings);
ZSPACE_TOOLSETS float widthOfGroup(const std::vector<float> &widths, int id);
ZSPACE_TOOLSETS bool appendBracingSegments(const zObjGraph &graph, const SlicerSettings &settings,
                           zScalar trim, std::vector<FieldSegment> &segments, int layer);

/*! \brief The width value computeWidthField() stores for a segment of the
 * given width: the width itself for a single-perimeter boundary, else (a
 * double-boundary bead or a bracing rib) 0.5 * width + 0.5 * overlapWidth,
 * must be finite and at least 1e-4 mm. */
ZSPACE_TOOLSETS zScalar printedWidth(zScalar width, bool isBoundary, bool doubleBoundary,
                     zScalar overlapWidth);

/*! \brief Groups actual boundary edges by width ID without adding connections.
 * Supports multiple loops, open chains and branches. edgeWidthIds follows
 * getEdgeData() order. Open boundaries require a separate domain sign field.
 * Appends one FieldSegment per width group and returns the count added. */
ZSPACE_TOOLSETS int loopWidthSegments(const zObjGraph &loop, const zIntArray &edgeWidthIds,
                      const std::vector<float> &widths,
                      std::vector<FieldSegment> &outSegments);

/*! \brief Builds one XY-perpendicular trim slot per connected polyline/loop
 * at clamped 3D arc ratio t, collecting all slots in outGraph. trimLength
 * is the positive half-length. Branching topology and non-finite parameters
 * are rejected; degenerate curves are skipped. False if no slot is created. */
ZSPACE_TOOLSETS bool createSlot(const zObjGraph &sourceGraph, zObjGraph &outGraph,
                             zScalar t, zScalar trimLength);

/*! \struct LoopAnchor
 *	\brief A point on a closed loop (edge i runs from vertex i to i + 1, the
 * last back to 0): its edge, the parameter on it and the loop's sense seen
 * from +Z. Found by locateLoopAnchor(), read by createAnchoredSlot().
 */
struct LoopAnchor {
  int edge = -1;          // -1 = not located
  double t = 0.0;         // 0..1 along edge
  zPoint point;           // the located point on the loop
  bool ccw = true;        // increasing index runs counter-clockwise seen from +Z
  double distance = 0.0;  // from the anchor point to the loop
  bool valid() const { return edge >= 0; }
};

/*! \brief Nearest point on the closed loop to anchor. The sense (ccw) is
 * read from senseLoop3D - the same loop in world coordinates, same vertex
 * order - when given (e.g. loop is its unrolled/chart image), else from loop
 * itself. False for fewer than 3 points. */
ZSPACE_TOOLSETS bool locateLoopAnchor(const zPointArray &loop, const zPoint &anchor,
                      LoopAnchor &out,
                      const zPointArray *senseLoop3D = nullptr);
// Selects the nearest actual closed component; edge IDs refer to source graph.
ZSPACE_TOOLSETS bool locateGraphLoopAnchor(const zObjGraph &graph, const zPoint &point,
                           zPointArray &outLoop, zIntArray &outEdgeIds,
                           LoopAnchor &outAnchor);

/*! \brief Boundary trim slot from an anchor: starts at anchor's point on
 * loop2D (the field-plane image of the loop it was located on, same vertex
 * order), moves offset along loop2D (positive = counter-clockwise seen from
 * +Z in 3D, negative = clockwise) and cuts perpendicular to the loop there,
 * halfLength to each side. outCenter (optional) gets the slot centre. False
 * for an invalid anchor or a degenerate loop. */
ZSPACE_TOOLSETS bool createAnchoredSlot(const zPointArray &loop2D, const LoopAnchor &anchor,
                        double offset, zScalar halfLength, zObjGraph &outSlot,
                        zPoint *outCenter = nullptr);

/*! \brief Slot-anchor step of both slicer operators' computeField()
 * (DoubleBoundary): locateLoopAnchor() on the floor boundary floorLoop3D
 * (world coordinates), logged as "[source] ...": red if hasAnchor is false
 * (missing input returns false; there is no default anchor) or nothing could be located, green with
 * the anchor's edge, parameter and distance otherwise. Each operator then
 * transfers out.point onto its layer charts with
 * zTsUtils::mapPointToMeshes() - the same barycentric transfer
 * zPrintSyn::alignStartPoints() uses for the seam - and re-locates it on
 * each chart's loop. */
ZSPACE_TOOLSETS bool locateSlotAnchor(const zPointArray &floorLoop3D, const zPoint &anchor,
                      bool hasAnchor, const std::string &source,
                      LoopAnchor &out);

/*! \brief One layer's boundary slot in both slicer operators: even layers
 * offsetLength counter-clockwise of anchor, odd layers clockwise (seen from
 * +Z), half-length w - see createAnchoredSlot(). offsetLength is
 * SlicerSettings::SlotOffset * w. False (empty outSlot) for an
 * invalid anchor. */
ZSPACE_TOOLSETS bool createParitySlot(const zPointArray &loop2D, const LoopAnchor &anchor,
                      bool isEven, double offsetLength, zScalar w,
                      zObjGraph &outSlot);

/*! \brief Expands [minBB, maxBB] about its centre by scaleFactor and snaps it
 * to a whole number of cellSize cells. False for an empty/degenerate box. */
ZSPACE_TOOLSETS bool computeGridBounds(const zPoint &minBB, const zPoint &maxBB,
                       double scaleFactor, double cellSize,
                       zDomain<zPoint> &outBounds, int &outResX,
                       int &outResY);

// The field functions below read only the grid's vertex positions and
// faces (the field is never written - results go to the output arrays), so
// one grid can be shared read-only by every thread.

/*! \brief Per-field-vertex unsigned distance to each segment graph.
 * Optional outBoundarySigns reuses the prepared type-0 edges for winding
 * (-1 inside, +1 outside). Use it only when those directed edges partition
 * the complete closed boundary exactly once (e.g. loopWidthSegments()).
 * Bracing edges are excluded. Recompute after geometry/grid changes. */
ZSPACE_TOOLSETS void computeSegmentDistances(const zObjMeshScalarField &field,
                             const std::vector<FieldSegment> &segments,
                             std::vector<zScalarArray> &outSegDistances,
                             zScalarArray *outBoundarySigns = nullptr);

// Inside/outside of the complete chart; independent of filtered side edges.
ZSPACE_TOOLSETS void computeMeshDomainSigns(const zObjMeshScalarField &field, const zObjMesh &mesh,
                             zScalarArray &outSigns);

/*! \brief Assigns each field vertex the print width (and index) of its
 * dominating (nearest, width-normalized) segment. */
ZSPACE_TOOLSETS void computeWidthField(const zObjMeshScalarField &field,
                       const std::vector<zScalarArray> &segDistances,
                       const std::vector<FieldSegment> &segments,
                       bool doubleBoundary, zScalar overlapWidth,
                       zScalarArray *outWidthField = nullptr,
                       zScalarArray *outIdField = nullptr);

/*! \brief Combines outer wall and internal bracing fields (with trim slots)
 * into the final SDF path field; the zero-isocontour is the print path.
 * Bracing ribs are cut perpendicular at their trim ratio (half-length
 * 1.5 * width). The boundary is cut only with doubleBoundary, by
 * boundarySlot (from createAnchoredSlot(), half-length 1.0 * width) with a
 * capsule sized from boundarySlotWidth, joining the inner and outer loop
 * into one chain; null or empty = no boundary slot. A single-perimeter wall
 * stays one closed loop. boundaryPolygon is the closed loop used for the
 * inside/outside test. outTrackIdField tags inner(0)/outer(1) for
 * doubleBoundary. Requested outputs are written directly; invalid/empty
 * input produces a positive path field (no contour), zero track IDs and
 * clears outTrimSlotsGraph. Optional boundarySigns supplies precomputed winding signs
 * for this grid and boundary (computeSegmentDistances()); null or a size
 * mismatch uses the boundaryPolygon query instead. */
ZSPACE_TOOLSETS void computePathField(const zObjMeshScalarField &field,
                      const zObjGraph &boundaryPolygon,
                      const std::vector<FieldSegment> &segments,
                      const std::vector<zScalarArray> &segDistances,
                      bool doubleBoundary, zScalar overlapWidth,
                      zObjGraph *outTrimSlotsGraph = nullptr,
                      zScalarArray *outPathField = nullptr,
                      zScalarArray *outTrackIdField = nullptr,
                      const zObjGraph *boundarySlot = nullptr,
                      zScalar boundarySlotWidth = 0.0,
                      const zScalarArray *boundarySigns = nullptr);

/*! \brief Raw marching-triangles zero-isocontour of pathVals over the field
 * grid (n-gon faces are triangle-fanned), interpolating widthVals and
 * thresholding trackIdVals (optional, may be empty) at every crossing.
 * Reads the arrays directly, never the field's own stored copy. */
// Single-layer extraction and direct chart-to-world mapping; no PrintSyn processing.
ZSPACE_TOOLSETS bool extractPrintPaths(const zObjMeshScalarField &fieldGrid,
    const zScalarArray &pathField, const zScalarArray &widthField,
    const zScalarArray &trackField, std::vector<PrintPath> &outPaths, int layer = -1,
    std::vector<std::string> *debug = nullptr);
ZSPACE_TOOLSETS bool mapPrintPaths(const std::vector<PrintPath> &paths2D,
    zObjMesh &chart2D, zObjMesh &mesh3D, const zFloatArray *vertexHeights,
    std::vector<PrintPath> &outPaths3D, int layer = -1,
    const zVectorArray *vertexNormals = nullptr);

ZSPACE_TOOLSETS void extractRawIsoContour(const zObjMeshScalarField &field,
                          const zScalarArray &pathVals,
                          const zScalarArray &widthVals,
                          const zScalarArray &trackIdVals,
                          zPointArray &outPositions, zIntArray &outEdgeConnects,
                          zFloatArray &outVertexWidths,
                          zFloatArray &outVertexTrackIds);

} // namespace zTsField

//--------------------------
//---- USD IO                           (utils/IO.cpp)
//--------------------------
// The only file IO in ts3dp. The slicer operators, zPrintSyn and zTsSlicing never
// call these - callers (sketches, tools) load/save with them and exchange
// in-memory objects only.
//
// Face tags travel as "usemtl" names: every "f" line takes the material of
// the last "usemtl" before it - "bottom", "top", "side" or "side_<n>" (wall,
// width group n -> SlicerSettings::BoundaryWidth[n]); stored in the USD input.
// Vertex colors on "v" lines are ignored.

namespace zTsIO {

// Copy an input USD and replace /TDP_gcode (ordered arrays), /TDP_printmesh
// (rectangle loft Mesh) and /TDP_printpath (linear BasisCurves).
// Custom arrays printNormals/printWidths/printHeight correspond 1:1 to points.
// Path coordinates/lengths are millimetres; export converts to the source stage units.
// Supply processed layerPaths to keep visualization separate across layers and
// components (no travel connectors). Omission renders the supplied path alone.
ZSPACE_TOOLSETS bool exportPrintPathUSD(const PrintPath &path, const std::string &inputFilePath,
                        const std::string &filePath, int layerCount,
                        std::string *outError = nullptr, int blockid = 0,
                        const PrintPathArray *layerPaths = nullptr);

// SettingsJSON preserves the existing {"slicer": {...}, "printSyn": {...}} schema.
// Returned coordinates and dimensional settings are millimetres after object transforms.
// MetersPerUnit records the source file unit; ratios, angles, tags and IDs are unchanged.
struct SlicerUSDInput {
  zObjMesh Block;
  zIntArray FaceTags;
  zIntArray FaceWidthIds; // Per-face widthid into settings BoundaryWidth.
  zIntArray BracingWidthIds; // Per-curve widthid into settings BracingWidth.
  zObjGraphArray BracingCurves; // One graph per USD curve, in authored order; coordinates in mm.
  zFloatArray BracingPrintWidths;
  zPoint StartPoint;
  std::string SettingsJSON;
  // Root customLayerData double[] printheightdomain, converted to mm. Empty
  // when absent. Input target; export overwrites metadata with measured output.
  std::vector<double> PrintHeightDomain;
  double MetersPerUnit = 0.01;
  std::string UpAxis;
};
// Requires ZSPACE_TS3DP_OPENUSD and matching OpenUSD headers/libraries.
// Paths must identify the geometry prim itself; missing/wrong schemas are errors.
// Failure clears the output and returns a diagnostic through outError.
ZSPACE_TOOLSETS bool readSlicerInputUSD(const std::string &filePath, SlicerUSDInput &out,
                       std::string *outError = nullptr,
                       const std::string &blockPath = "/AG_block",
                       const std::string &bracingPath = "/Slicer_bracing",
                       const std::string &startPath = "/Slicer_slot");

} // namespace zTsIO

/** @}*/

/** @}*/

} // namespace zSpace


#endif // ZSPACE_TS_3DP_UTILS_H
