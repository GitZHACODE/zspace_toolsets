#ifndef EIGEN_DONT_PARALLELIZE
#define EIGEN_DONT_PARALLELIZE
#endif

#include <zspace/zToolsets/fabrication/slicing/utils/Parallel.h>
#include <zspace/zToolsets/fabrication/slicing/zSlicer/HarmonicSlicer.h>

#ifndef ZSPACE_INLINE
#define ZSPACE_INLINE ZSPACE_TOOLSETS_INLINE
#endif

#ifndef TETLIBRARY
#define TETLIBRARY
#endif
#include <tetgen.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <limits>
#include <sstream>
#include <unordered_map>

namespace zSpace {

using namespace zTsUtils;

// Input face tags/width IDs become TetGen facet markers. Rebuild the refined
// surface with shared volume vertex IDs, then build full edge/face incidence.
ZSPACE_INLINE bool HarmonicSlicer::buildTetVolume() {
  tetgenio input;
  tetgenio output;
  input.firstnumber = 0;
  input.numberofpoints = static_cast<int>(_Tet.Surface.Positions.size());
  input.pointlist = new REAL[input.numberofpoints * 3];
  input.pointmarkerlist = new int[input.numberofpoints];
  for (int i = 0; i < input.numberofpoints; ++i) {
    const zPoint &point = _Tet.Surface.Positions[i];
    input.pointlist[i * 3 + 0] = point.x;
    input.pointlist[i * 3 + 1] = point.y;
    input.pointlist[i * 3 + 2] = point.z;
    input.pointmarkerlist[i] = _VertexTags[i];
  }

  input.numberoffacets = static_cast<int>(_Tet.Surface.Triangles.size());
  input.facetlist = new tetgenio::facet[input.numberoffacets];
  input.facetmarkerlist = new int[input.numberoffacets];
  std::vector<std::pair<int,int>> sideGroups;
  const auto faceGroup = [&](int t) {
    const auto tag = _Tet.SurfaceTags[t];
    if (tag == HarmonicField::TetFaceTag::Bottom) return 0;
    if (tag == HarmonicField::TetFaceTag::Top) return 1;
    const auto key = std::make_pair(_Tet.SurfaceWidthIds[t], _Tet.SurfaceRegionTags[t]);
    auto it = std::find(sideGroups.begin(), sideGroups.end(), key);
    if (it == sideGroups.end()) { sideGroups.push_back(key); return 1 + int(sideGroups.size()); }
    return 2 + int(it - sideGroups.begin());
  };
  for (size_t t = 0; t < _Tet.SurfaceTags.size(); ++t) {
    if (_Tet.SurfaceTags[t] == HarmonicField::TetFaceTag::Interior ||
        (_Tet.SurfaceTags[t] == HarmonicField::TetFaceTag::Side && _Tet.SurfaceWidthIds[t] < 0)) {
      _Status = "Input surface has an invalid face tag or side width ID.";
      printDebug("[buildTetVolume] " + _Status);
      return false;
    }
  }
  for (int i = 0; i < input.numberoffacets; ++i) {
    tetgenio::facet &facet = input.facetlist[i];
    facet.numberofpolygons = 1;
    facet.polygonlist = new tetgenio::polygon[1];
    facet.numberofholes = 0;
    facet.holelist = nullptr;
    tetgenio::polygon &polygon = facet.polygonlist[0];
    polygon.numberofvertices = 3;
    polygon.vertexlist = new int[3];
    for (int k = 0; k < 3; ++k) polygon.vertexlist[k] = _Tet.Surface.Triangles[i][k];
    input.facetmarkerlist[i] = faceGroup(i) + 1; // 0 means unmarked in TetGen
  }

  // Bound tetrahedron size/quality so the interior isn't just a handful of
  // large Delaunay tets spanning the whole volume: "-q1.4" caps the
  // radius-edge ratio (element shape), "-a<volume>" caps each tet's volume,
  // both driving TetGen to insert interior Steiner points. Target volume is
  // that of a regular tetrahedron whose edge equals the input mesh's own
  // mean edge length, so tet density roughly matches the input surface's
  // own resolution. Without "-Y" the quality/volume bounds may also split
  // boundary facets/edges - those surface points are tagged from the facet
  // markers and the surface rebuilt below.
  const double targetEdgeLength = meanEdgeLength(_Tet.Surface.Positions, _Tet.Surface.Triangles);
  std::string switchesStr = "pQ";
  if (targetEdgeLength > 1e-6) {
    const double targetVolume = (targetEdgeLength * targetEdgeLength * targetEdgeLength) / (6.0 * std::sqrt(2.0));
    char volumeText[64];
    std::snprintf(volumeText, sizeof(volumeText), "%.10f", targetVolume);
    switchesStr = std::string("pq1.4a") + volumeText + "Q";
  }
  std::vector<char> switches(switchesStr.begin(), switchesStr.end());
  switches.push_back('\0');

  try {
    tetrahedralize(switches.data(), &input, &output);
  } catch (...) {
    _Status = "computeBasicStructure failed: TetGen threw while tetrahedralizing the capped mesh.";
    printDebug("[computeBasicStructureImpl] " + _Status);
    return false;
  }

  if (output.numberofpoints <= 0 || output.numberoftetrahedra <= 0 || output.numberofcorners < 4) {
    _Status = "computeBasicStructure failed: TetGen returned an empty or invalid tetrahedralization.";
    printDebug("[computeBasicStructureImpl] " + _Status);
    return false;
  }
  if (output.numberoftrifaces <= 0 || !output.trifacelist || !output.trifacemarkerlist) {
    _Status = "computeBasicStructure failed: TetGen returned no marked boundary faces.";
    printDebug("[computeBasicStructureImpl] " + _Status);
    return false;
  }

  const int outputCount = output.numberofpoints;
  const int inputCount = input.numberofpoints;
  if (outputCount < inputCount) {
    _Status = "TetGen returned fewer vertices than the input; input vertex correspondence cannot be preserved.";
    printDebug("[computeBasicStructureImpl] " + _Status);
    return false;
  }
  zPointArray outPositions(outputCount);
  for (int i = 0; i < outputCount; ++i)
    outPositions[i] = zPoint(output.pointlist[i * 3 + 0], output.pointlist[i * 3 + 1], output.pointlist[i * 3 + 2]);
  std::vector<std::array<int, 4>> outTets(output.numberoftetrahedra);
  for (int i = 0; i < output.numberoftetrahedra; ++i) {
    const int *tet = output.tetrahedronlist + (i * output.numberofcorners);
    outTets[i] = {tet[0], tet[1], tet[2], tet[3]};
  }

  // Boundary faces = faces used by exactly one tet, oriented away from that
  // tet's fourth vertex (outward).
  std::map<std::array<int, 3>, std::pair<int, std::array<int, 3>>> faceUse;
  static const int kFaceCorners[4][4] = {{1, 2, 3, 0}, {0, 2, 3, 1}, {0, 1, 3, 2}, {0, 1, 2, 3}};
  for (const std::array<int, 4> &tet : outTets) {
    for (const auto &corners : kFaceCorners) {
      std::array<int, 3> face = {tet[corners[0]], tet[corners[1]], tet[corners[2]]};
      std::array<int, 3> sorted = face;
      std::sort(sorted.begin(), sorted.end());
      auto &use = faceUse[sorted];
      if (use.first++ == 0) {
        const Eigen::Vector3d a = toEigenVec(outPositions[face[0]]);
        const Eigen::Vector3d normal = (toEigenVec(outPositions[face[1]]) - a).cross(toEigenVec(outPositions[face[2]]) - a);
        if (normal.dot(toEigenVec(outPositions[tet[corners[3]]]) - a) > 0.0) std::swap(face[1], face[2]);
        use.second = face;
      }
    }
  }
  std::vector<char> onSurface(outputCount, 0);
  std::vector<std::array<int, 3>> boundaryFaces;
  for (const auto &entry : faceUse) {
    if (entry.second.first != 1) continue;
    boundaryFaces.push_back(entry.second.second);
    for (int v : entry.second.second) onSurface[v] = 1;
  }

  // Renumber: input points (ids unchanged), new surface points, interior.
  std::vector<int> newId(outputCount, -1);
  int nextId = 0;
  for (int i = 0; i < inputCount; ++i) newId[i] = nextId++;
  for (int i = inputCount; i < outputCount; ++i)
    if (onSurface[i]) newId[i] = nextId++;
  const int surfaceCount = nextId;
  for (int i = inputCount; i < outputCount; ++i)
    if (!onSurface[i]) newId[i] = nextId++;

  // The face group of every TetGen boundary face (its facet marker - 1), by
  // sorted face key.
  std::unordered_map<std::array<int, 3>, int, TriangleKeyHash> groupOfFace;
  for (int f = 0; f < output.numberoftrifaces; ++f) {
    std::array<int, 3> key = {output.trifacelist[3 * f], output.trifacelist[3 * f + 1], output.trifacelist[3 * f + 2]};
    std::sort(key.begin(), key.end());
    const int group = output.trifacemarkerlist[f] - 1;
    if (group >= 0) groupOfFace[key] = group;
  }

  // Width group of every boundary face: its face group - 2 (-1 on a cap). A
  // new surface point is Bottom if it touches a Bottom face, Top if it
  // touches a Top face, Side otherwise (so a point on a floor-cap border or
  // inside the cap stays a Dirichlet Bottom vertex, one on the wall is free).
  // Input points keep their tags; interior points are free.
  zIntArray boundaryFaceWidthIds(boundaryFaces.size(), -1), boundaryRegionTags(boundaryFaces.size(), 0);
  std::vector<HarmonicField::TetFaceTag> boundaryTags(boundaryFaces.size(), HarmonicField::TetFaceTag::Interior);
  std::vector<char> touchesBottom(outputCount, 0), touchesTop(outputCount, 0);
  int unmarkedFaces = 0;
  for (size_t f = 0; f < boundaryFaces.size(); ++f) {
    std::array<int, 3> key = boundaryFaces[f];
    std::sort(key.begin(), key.end());
    const auto found = groupOfFace.find(key);
    if (found == groupOfFace.end()) {
      ++unmarkedFaces;
      continue;
    }
    const int group = found->second;
    boundaryTags[f] = group == 0 ? HarmonicField::TetFaceTag::Bottom : group == 1 ? HarmonicField::TetFaceTag::Top : HarmonicField::TetFaceTag::Side;
    if (group >= 2) { boundaryFaceWidthIds[f] = sideGroups.at(group-2).first; boundaryRegionTags[f] = sideGroups.at(group-2).second; }
    for (int v : boundaryFaces[f]) {
      if (group == 0) touchesBottom[v] = 1;
      else if (group == 1) touchesTop[v] = 1;
    }
  }
  if (unmarkedFaces > 0) {
    _Status = "Boundary faces have no TetGen facet marker; surface classification is incomplete.";
    printDebug("[buildTetVolume] " + _Status);
    return false;
  }

  _Tet.Vertices.assign(outputCount, zTetrahedron::TetVertex());
  _Tet.vertextag.assign(outputCount, static_cast<int>(VertexTag::Interior));
  for (int i = 0; i < outputCount; ++i) {
    const int id = newId[i];
    _Tet.Vertices[id].Position = outPositions[i];
    if (i < inputCount)
      _Tet.vertextag[id] = (i < static_cast<int>(_VertexTags.size())) ? _VertexTags[i]
                                                                     : static_cast<int>(VertexTag::Side);
    else if (onSurface[i])
      _Tet.vertextag[id] = static_cast<int>(touchesBottom[i] ? VertexTag::Bottom
                                          : touchesTop[i]  ? VertexTag::Top
                                                           : VertexTag::Side);
  }
  _Tet.Volumes.clear();
  _Tet.Volumes.reserve(outTets.size());
  for (const std::array<int, 4> &tet : outTets) {
    zTetrahedron::TetVolume volume;
    volume.VertexIds = {{newId[tet[0]], newId[tet[1]], newId[tet[2]], newId[tet[3]]}};
    _Tet.Volumes.push_back(volume);
  }

  // Rebuild the surface from the tet boundary (same ids as the tet mesh).
  TriangleSoup surface;
  surface.Positions.reserve(surfaceCount);
  for (int v = 0; v < surfaceCount; ++v) surface.Positions.push_back(_Tet.Vertices[v].Position);
  surface.Triangles.reserve(boundaryFaces.size());
  for (const std::array<int, 3> &face : boundaryFaces)
    surface.Triangles.push_back({newId[face[0]], newId[face[1]], newId[face[2]]});
  _Tet.Surface = std::move(surface);
  _Tet.SurfaceWidthIds = boundaryFaceWidthIds;
  _Tet.SurfaceRegionTags = boundaryRegionTags;
  _Tet.SurfaceTags = std::move(boundaryTags);

  return !_Tet.Volumes.empty() && buildTetTopology();
}

ZSPACE_INLINE bool HarmonicSlicer::buildTetTopology() {
  for (auto &vertex : _Tet.Vertices) {
    vertex.EdgeIds.clear();
    vertex.FaceIds.clear();
    vertex.VolumeIds.clear();
  }
  _Tet.Edges.clear();
  _Tet.Faces.clear();
  _Tet.FaceBasisGradients.clear();
  const float invalid = std::numeric_limits<float>::quiet_NaN();
  _Tet.UVW.assign(_Tet.Vertices.size(), zPoint(invalid, invalid, invalid));
  _Tet.UVWValid.assign(_Tet.Vertices.size(), false);
  _Tet.FieldValid = false;
  std::unordered_map<unsigned long long, int> edges;
  std::unordered_map<std::array<int, 3>, int, TriangleKeyHash> faces;
  const int edgeCorners[6][2] = {{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}};
  const int faceCorners[4][3] = {{1,2,3},{0,3,2},{0,1,3},{0,2,1}};
  for (int t = 0; t < static_cast<int>(_Tet.Volumes.size()); ++t) {
    const auto &tet = _Tet.Volumes[t].VertexIds;
    auto &volume = _Tet.Volumes[t];
    for (int v : tet) _Tet.Vertices[v].VolumeIds.push_back(t);
    for (int e = 0; e < 6; ++e) {
      const int a = tet[edgeCorners[e][0]], b = tet[edgeCorners[e][1]];
      const auto key = makeEdgeKey(a, b);
      auto inserted = edges.emplace(key, static_cast<int>(_Tet.Edges.size()));
      const int id = inserted.first->second;
      if (inserted.second) {
        zTetrahedron::TetEdge edge;
        edge.VertexIds = {{(std::min)(a,b), (std::max)(a,b)}};
        _Tet.Edges.push_back(edge);
        _Tet.Vertices[a].EdgeIds.push_back(id);
        _Tet.Vertices[b].EdgeIds.push_back(id);
      }
      _Tet.Edges[id].VolumeIds.push_back(t);
      volume.EdgeIds[e] = id;
    }
    for (int f = 0; f < 4; ++f) {
      std::array<int,3> tri = {{tet[faceCorners[f][0]],tet[faceCorners[f][1]],tet[faceCorners[f][2]]}};
      auto key = tri;
      std::sort(key.begin(), key.end());
      auto inserted = faces.emplace(key, static_cast<int>(_Tet.Faces.size()));
      const int id = inserted.first->second;
      if (inserted.second) {
        zTetrahedron::TetFace face;
        const auto a = toEigenVec(_Tet.Vertices[tri[0]].Position);
        if ((toEigenVec(_Tet.Vertices[tri[1]].Position)-a).cross(toEigenVec(_Tet.Vertices[tri[2]].Position)-a)
            .dot(toEigenVec(_Tet.Vertices[tet[f]].Position)-a) > 0) std::swap(tri[1],tri[2]);
        face.VertexIds = tri;
        face.VolumeIds[0] = t;
        Eigen::Matrix3d gradients = Eigen::Matrix3d::Zero();
        computeTriElement(_Tet.Vertices, tri, face, gradients);
        _Tet.FaceBasisGradients.push_back(gradients);
        for (int e = 0; e < 3; ++e) {
          face.EdgeIds[e] = edges.at(makeEdgeKey(tri[e],tri[(e+1)%3]));
          _Tet.Edges[face.EdgeIds[e]].FaceIds.push_back(id);
          _Tet.Vertices[tri[e]].FaceIds.push_back(id);
        }
        _Tet.Faces.push_back(face);
      } else {
        if (_Tet.Faces[id].VolumeIds[1] >= 0) {
          _Status = "A tetrahedral face has more than two incident volumes.";
          printDebug("[buildTetTopology] " + _Status);
          return false;
        }
        _Tet.Faces[id].VolumeIds[1] = t;
      }
      volume.FaceIds[f] = id;
    }
  }
  _Tet.FaceTags.assign(_Tet.Faces.size(), HarmonicField::TetFaceTag::Interior);
  _Tet.FaceWidthIds.assign(_Tet.Faces.size(), -1);
  _Tet.FaceRegionTags.assign(_Tet.Faces.size(), 0);
  _Tet.SurfaceFaceIds.resize(_Tet.Surface.Triangles.size());
  for (size_t s = 0; s < _Tet.Surface.Triangles.size(); ++s) {
    auto key = _Tet.Surface.Triangles[s];
    std::sort(key.begin(), key.end());
    const auto found = faces.find(key);
    if (found == faces.end() || _Tet.Faces[found->second].VolumeIds[1] >= 0) {
      _Status = "Surface triangle does not correspond to a tetrahedral boundary face.";
      printDebug("[buildTetTopology] " + _Status);
      return false;
    }
    const int faceId = found->second;
    _Tet.FaceTags[faceId] = _Tet.SurfaceTags[s];
    _Tet.FaceWidthIds[faceId] = _Tet.FaceTags[faceId] == HarmonicField::TetFaceTag::Side ? _Tet.SurfaceWidthIds[s] : -1;
    _Tet.FaceRegionTags[faceId] = _Tet.SurfaceRegionTags[s];
    _Tet.SurfaceFaceIds[s] = found->second;
  }
  _Tet.SurfaceWidthIds.clear();
  _Tet.SurfaceRegionTags.clear();
  _Tet.SurfaceTags.clear();
  return true;
}

ZSPACE_INLINE void HarmonicSlicer::updateFieldGradients() {
  _Tet.VolumeFieldGradients.assign(_Tet.Volumes.size(), Eigen::Vector3d::Zero());
  _Tet.FaceFieldGradients.assign(_Tet.Faces.size(), Eigen::Vector3d::Zero());
  if (!_Tet.FieldValid) return;
  for (size_t t = 0; t < _Tet.Volumes.size(); ++t) {
    if (!_Tet.Volumes[t].Valid) continue;
    Eigen::Vector4d values;
    for (int k = 0; k < 4; ++k) values[k] = _Tet.UVW[_Tet.Volumes[t].VertexIds[k]].z;
    _Tet.VolumeFieldGradients[t] = _Tet.VolumeBasisGradients[t] * values;
  }
  for (size_t f = 0; f < _Tet.Faces.size(); ++f) {
    const auto &face = _Tet.Faces[f];
    if (!face.Valid) continue;
    Eigen::Vector3d values;
    for (int k = 0; k < 3; ++k) values[k] = _Tet.UVW[face.VertexIds[k]].z;
    _Tet.FaceFieldGradients[f] = _Tet.FaceBasisGradients[f] * values;
  }
}

// Linear tet element used by the harmonic solve (Method 2, step 3) and the
// gradient-walking trace (step 4) - pure per-tet geometry, independent of
// any boundary values, so it's built once here alongside the tet volume
// itself rather than being recomputed by whichever step needs it.
ZSPACE_INLINE void HarmonicSlicer::buildTetElements(int &degenerateTets, bool useParallel) {
  for (auto &volume : _Tet.Volumes) { volume.Valid = false; volume.Volume = 0.0; }
  _Tet.VolumeBasisGradients.assign(_Tet.Volumes.size(), Eigen::Matrix<double, 3, 4>::Zero());
  _Tet.VolumeStiffness.assign(_Tet.Volumes.size(), Eigen::Matrix4d::Zero());
  _Tet.FieldValid = false;
  updateFieldGradients();
  int degenerateCount = 0;
  zTsDetail::ParallelErrors parallelErrors1(_Tet.Volumes.size());
#pragma omp parallel for schedule(static) if(useParallel) reduction(+:degenerateCount)
  for (int t = 0; t < static_cast<int>(_Tet.Volumes.size()); ++t) {
    try {
      const std::array<int, 4> &tet = _Tet.Volumes[t].VertexIds;
      zTetrahedron::TetVolume &element = _Tet.Volumes[t];
      Eigen::Matrix4d interpolation;
      Eigen::Matrix3d edgeMatrix;
      for (int local = 0; local < 4; ++local) {
        const zPoint &point = _Tet.Vertices[tet[local]].Position;
        interpolation.row(local) << 1.0, point.x, point.y, point.z;
        if (local > 0) {
          const zPoint &origin = _Tet.Vertices[tet[0]].Position;
          edgeMatrix.col(local - 1) << point.x - origin.x, point.y - origin.y, point.z - origin.z;
        }
      }
      const double determinant = edgeMatrix.determinant();
      element.Valid = std::fabs(determinant) >= 1e-5;
      if (!element.Valid) {
        ++degenerateCount;
        continue;
      }
      element.Volume = std::fabs(determinant) / 6.0;
      const Eigen::Matrix4d coefficients = interpolation.inverse();
      for (int local = 0; local < 4; ++local) _Tet.VolumeBasisGradients[t].col(local) = coefficients.block<3, 1>(1, local);
      _Tet.VolumeStiffness[t] = element.Volume * _Tet.VolumeBasisGradients[t].transpose() * _Tet.VolumeBasisGradients[t];
    } catch (...) { parallelErrors1.capture(t); }
  }
  parallelErrors1.rethrow();
  degenerateTets = degenerateCount;
}

//--------------------------
//---- STATE RESETS
//--------------------------

ZSPACE_INLINE void HarmonicSlicer::clearHeat() {
  _Tet.HeatValue.first.clear();
  _Tet.HeatValue.second.clear();
}

//--------------------------
//---- computeBasicStructure()
//--------------------------

// The input mesh as the triangle soup (polygons fanned), vertices in the
// input order - so _VertexTags tags them directly - with the width group of
// every triangle's input face (-1 for a Bottom/Top face).
ZSPACE_INLINE void HarmonicSlicer::buildSurfaceSoup() {
  zFnMesh fn(*i_BlockMesh);
  _Tet.Surface = TriangleSoup();
  fn.getVertexPositions(_Tet.Surface.Positions);

  zIntArray connects;
  zIntArray counts;
  fn.getPolygonData(connects, counts);
  _Tet.SurfaceWidthIds.clear();
  _Tet.SurfaceRegionTags.clear();
  _Tet.SurfaceTags.clear();
  int cursor = 0;
  for (int face = 0; face < static_cast<int>(counts.size()); ++face) {
    const bool side = face < static_cast<int>(i_FaceTags.size()) && isSideFaceTag(i_FaceTags[face]);
    const int widthId = !side ? -1 : (face < static_cast<int>(i_FaceWidthIds.size()) ? i_FaceWidthIds[face] : 0);
    const int count = counts[face];
    for (int i = 1; i + 1 < count; ++i) {
      const std::array<int, 3> tri = {connects[cursor], connects[cursor + i], connects[cursor + i + 1]};
      if (tri[0] != tri[1] && tri[1] != tri[2] && tri[2] != tri[0]) {
        _Tet.Surface.Triangles.push_back(tri);
        _Tet.SurfaceWidthIds.push_back(widthId);
        _Tet.SurfaceRegionTags.push_back(side ? i_FaceTags[face] : 0);
        _Tet.SurfaceTags.push_back(side ? HarmonicField::TetFaceTag::Side :
            i_FaceTags[face] == static_cast<int>(VertexTag::Bottom) ? HarmonicField::TetFaceTag::Bottom :
            i_FaceTags[face] == static_cast<int>(VertexTag::Top) ? HarmonicField::TetFaceTag::Top : HarmonicField::TetFaceTag::Interior);
      }
    }
    cursor += count;
  }
}

ZSPACE_INLINE void HarmonicSlicer::clearPrintPaths() {
  o_PrintPaths2D.clear();
  o_PrintPaths3D.clear();
  o_ResultantPrintPath = PrintPath();
  resetPrintStatistics();
}

ZSPACE_INLINE void HarmonicSlicer::clearFields() {
  o_FaceTagSources.clear();
  clearPrintPaths();
  o_PathFields.clear();
  o_WidthFields.clear();
  o_TrackFields.clear();
  o_BracingGraphs2D.clear();
  zFnMeshScalarField(o_ScalarField).clear();
}

ZSPACE_INLINE void HarmonicSlicer::clearUnroll() {
  clearFields();
  o_IsoMeshes2D.clear();
  o_IsoMeshBoundary2D.clear();
  o_IsoMeshBoundaryWidthId.clear();
}

ZSPACE_INLINE void HarmonicSlicer::clearContours() {
  clearUnroll();
  o_IsoMeshes.clear();
  o_IsoMeshBoundary.clear();
  o_IsoMeshUV.clear();
  o_VertexLayerHeights.clear();
  _HarmonicIsoEdgeWidthIds.clear();
}

ZSPACE_INLINE void HarmonicSlicer::clearParameterization() {
  clearContours();
  o_TraceGraphs.clear();
  _BottomUV.clear();
  _BottomUVScale = zPoint(1.0, 1.0, 0.0);
  _BottomUVOrigin = zPoint(0.0, 0.0, 0.0);
  _Tet.UVWValid.assign(_Tet.UVW.size(), false);
  const float invalid = std::numeric_limits<float>::quiet_NaN();
  for (zPoint &p : _Tet.UVW) { p.x = invalid; p.y = invalid; }
}

ZSPACE_INLINE void HarmonicSlicer::clearStructure() {
  clearParameterization();
  clearHeat();
  _Tet = HarmonicField();
  _VertexTags.clear();
}

ZSPACE_INLINE bool HarmonicSlicer::computeBasicStructure(bool useParallel) {
  const bool success = computeBasicStructureImpl(useParallel);
  if (success) zTsUtils::printSuccess("[HarmonicSlicer::computeBasicStructure] Success");
  else zTsUtils::printError("[HarmonicSlicer::computeBasicStructure] Failed");
  return success;
}

ZSPACE_INLINE bool HarmonicSlicer::computeBasicStructureImpl(bool useParallel) {
  // Drop every previous result up front: if this call fails (e.g. missing
  // tags), no later step may run on the last successful structure.
  clearStructure();

  if (!i_SlotPoint) {
    _Status = "failed: i_SlotPoint is required. Call setSlotAnchor() before computing.";
    printDebug("[computeBasicStructureImpl] " + _Status);
    return false;
  }

  if (!isLoaded()) {
    _Status = "computeBasicStructure failed: i_BlockMesh is null.";
    printDebug("[computeBasicStructureImpl] " + _Status);
    return false;
  }
  zFnMesh fnCheck(*i_BlockMesh);
  const int inputVertexCount = fnCheck.numVertices();
  if (inputVertexCount == 0) {
    _Status = "computeBasicStructure failed: input mesh has no vertices.";
    printDebug("[computeBasicStructureImpl] " + _Status);
    return false;
  }
  for (zItMeshHalfEdge he(*i_BlockMesh); !he.end(); he++) {
    if (he.onBoundary()) {
      _Status = "computeBasicStructure failed: input mesh is not closed - only closed meshes are supported.";
      printDebug("[computeBasicStructureImpl] " + _Status);
      return false;
    }
  }
  std::string tagError;
  if (!zTsUtils::vertexTagsFromFaceTags(*i_BlockMesh, i_FaceTags, _VertexTags, &tagError)) {
    _Status = "computeBasicStructure failed: face tags: " + tagError +
              " - pass one tag per face to setInputMesh(mesh, faceTags).";
    printDebug("[computeBasicStructureImpl] " + _Status);
    return false;
  }

  buildSurfaceSoup();

  int bottomCount = 0;
  int topCount = 0;
  for (int tag : _VertexTags) {
    if (tag == static_cast<int>(VertexTag::Bottom)) ++bottomCount;
    else if (tag == static_cast<int>(VertexTag::Top)) ++topCount;
  }

  if (!buildTetVolume()) {
    clearStructure(); // Keep the failure status but discard incomplete topology.
    return false;
  }
  int degenerateTets = 0;
  buildTetElements(degenerateTets, useParallel);

  // Layering basis: the shortest / longest Bottom->Top heat geodesic on the
  // (rebuilt, Steiner-including) surface - see computeHeatDistances().
  if (!computeHeatDistances()) return false; // _Status/printError already set

  std::ostringstream out;
  out << "computeBasicStructure: vertices=" << _Tet.Surface.Positions.size()
      << " triangles=" << _Tet.Surface.Triangles.size()
      << " bottom-tagged=" << bottomCount << " top-tagged=" << topCount
      << " | tet-vertices=" << _Tet.Vertices.size() << " tets=" << _Tet.Volumes.size()
      << " surface Steiner points=" << (static_cast<int>(_Tet.Surface.Positions.size()) - inputVertexCount)
      << " degenerate-tets=" << degenerateTets;
  _Status = out.str();

  return true;
}

} // namespace zSpace
