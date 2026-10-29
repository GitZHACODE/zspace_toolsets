#ifndef EIGEN_DONT_PARALLELIZE
#define EIGEN_DONT_PARALLELIZE
#endif

#include <zspace/zToolsets/fabrication/slicing/utils/Parallel.h>
#include <zspace/zToolsets/fabrication/slicing/zSlicer/HarmonicSlicer.h>

#ifndef ZSPACE_INLINE
#define ZSPACE_INLINE ZSPACE_TOOLSETS_INLINE
#endif

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <unordered_map>
#include <unordered_set>

namespace zSpace {

using namespace zTsUtils;

//--------------------------
//---- METHOD 2 EXTENSION : DIRECT-PRODUCT (u,v,w) PARAMETERIZATION
//--------------------------

// PCA-projects the tet volume's Bottom-tagged vertices onto their own
// best-fit plane and normalizes the result to [0,1] (_BottomUV, with
// _BottomUVScale / _BottomUVOrigin to go back to physical units) - the chart
// every trace ends on.
ZSPACE_INLINE bool HarmonicSlicer::buildBottomParameterization() {
  _BottomUV.assign(_Tet.Vertices.size(), zPoint());
  zIntArray bottomIds;
  for (size_t id = 0; id < _Tet.vertextag.size(); ++id)
    if (_Tet.vertextag[id] == static_cast<int>(VertexTag::Bottom))
      bottomIds.push_back(static_cast<int>(id));
  if (bottomIds.empty()) return false;

  Eigen::Vector3d centroid = Eigen::Vector3d::Zero();
  for (int id : bottomIds) centroid += toEigenVec(_Tet.Vertices[id].Position);
  centroid /= static_cast<double>(bottomIds.size());
  Eigen::Matrix3d covariance = Eigen::Matrix3d::Zero();
  for (int id : bottomIds) {
    const Eigen::Vector3d delta = toEigenVec(_Tet.Vertices[id].Position) - centroid;
    covariance += delta * delta.transpose();
  }
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(covariance);
  if (solver.info() != Eigen::Success || !solver.eigenvectors().allFinite()) return false;
  // Eigenvalues ascending: columns 2 and 1 span the best-fit plane.
  const Eigen::Vector3d axisU = solver.eigenvectors().col(2);
  const Eigen::Vector3d axisV = solver.eigenvectors().col(1);

  double uMin = std::numeric_limits<double>::max(), uMax = -std::numeric_limits<double>::max();
  double vMin = uMin, vMax = uMax;
  for (int id : bottomIds) {
    const Eigen::Vector3d delta = toEigenVec(_Tet.Vertices[id].Position) - centroid;
    _BottomUV[id] = zPoint(delta.dot(axisU), delta.dot(axisV), 0.0);
    uMin = std::min(uMin, static_cast<double>(_BottomUV[id].x));
    uMax = std::max(uMax, static_cast<double>(_BottomUV[id].x));
    vMin = std::min(vMin, static_cast<double>(_BottomUV[id].y));
    vMax = std::max(vMax, static_cast<double>(_BottomUV[id].y));
  }
  const double uRange = (uMax - uMin > 1e-6) ? (uMax - uMin) : 1.0;
  const double vRange = (vMax - vMin > 1e-6) ? (vMax - vMin) : 1.0;
  for (int id : bottomIds) {
    zPoint &p = _BottomUV[id];
    p.x = static_cast<float>((p.x - uMin) / uRange);
    p.y = static_cast<float>((p.y - vMin) / vRange);
  }
  _BottomUVScale = zPoint(static_cast<float>(uRange), static_cast<float>(vRange), 0.0);
  _BottomUVOrigin = zPoint(uMin, vMin, 0.0);

  return true;
}

// Walks tet-to-tet from vertexId's own barycentric one-hot coordinate,
// stepping along -grad(phi) inside each tet. The direction used at each
// point is interpolated from nodalGradient (a volume-weighted, per-vertex
// averaged gradient - see computeVolumetricParameterization()) rather than
// the current tet's own constant gradient: raw P1 element gradients are
// discontinuous across tets, so seeding/continuing a trace from a single
// element's constant gradient makes the walk sensitive to which of a
// shared vertex's many incident tets happened to be picked, producing
// local (u,v) jumps between neighboring, otherwise similar, traces. The
// nodal field is continuous (piecewise-linear), so the walk direction
// varies smoothly instead. A P1 tet's barycentric coordinates are still
// affine in the ACTUAL (unsmoothed) element gradients, so their directional
// derivative along a unit direction d is exactly Gradients.col(i).dot(d)
// regardless of which direction is chosen; the first coordinate to reach
// zero marks the face the walk exits through. Crossing into the neighbor
// tet just carries the shared-face weights across (same global vertex id ->
// same weight); no re-solving is needed. The walk terminates on a boundary
// face whose 3 vertices are all Bottom-tagged, and returns the
// barycentric-interpolated (u,v) there.
ZSPACE_INLINE bool HarmonicSlicer::traceToBottom(int vertexId, const VolumeTraceData &data,
    zPoint &outUV, zPointArray *outPath, std::string *failureReason) const {
  const auto &vertexTet = data.Seeds;
  const auto &faceAdjacency = data.Adjacency;
  const auto &nodalGradient = data.Gradient;
  if (failureReason) failureReason->clear();
  const auto fail = [&](const char *reason) {
    if (failureReason) *failureReason = reason;
    return false;
  };
  if (outPath) {
    outPath->clear();
    outPath->push_back(_Tet.Vertices[vertexId].Position);
  }
  int tetIndex = vertexTet[vertexId];
  if (tetIndex < 0) return fail("NoStartTet");

  std::array<int, 4> tet = _Tet.Volumes[tetIndex].VertexIds;
  Eigen::Vector4d b = Eigen::Vector4d::Zero();
  for (int local = 0; local < 4; ++local) if (tet[local] == vertexId) b[local] = 1.0;

  const int maxSteps = i_Settings._HTraceMaxSteps;
  for (int step = 0; step < maxSteps; ++step) {
    const zTetrahedron::TetVolume &element = _Tet.Volumes[tetIndex];
    const auto &gradients = _Tet.VolumeBasisGradients[tetIndex];
    if (!element.Valid) return fail("DegenerateElement");

    Eigen::Vector3d gradAtPoint = Eigen::Vector3d::Zero();
    for (int local = 0; local < 4; ++local) gradAtPoint += b[local] * nodalGradient[tet[local]];
    const double gradNorm = gradAtPoint.norm();
    if (!gradAtPoint.allFinite() || !std::isfinite(gradNorm)) return fail("NonFiniteGradient");
    if (gradNorm < 1e-15) return fail("ZeroGradient");
    const Eigen::Vector3d direction = -gradAtPoint / gradNorm;

    Eigen::Vector4d rate;
    for (int local = 0; local < 4; ++local) rate[local] = gradients.col(local).dot(direction);

    int exitLocal = -1;
    double exitT = std::numeric_limits<double>::max();
    for (int local = 0; local < 4; ++local) {
      if (rate[local] >= -1e-15) continue;
      const double t = b[local] / -rate[local];
      if (t < exitT) {
        exitT = t;
        exitLocal = local;
      }
    }
    if (exitLocal < 0) return fail("NoExit");

    // |grad(b_i)| is the reciprocal altitude opposite vertex i.
    // Limit each step locally, then re-evaluate the direction even when
    // the point remains inside the same simplex. Never step across a face.
    const double maxGradient = gradients.colwise().norm().maxCoeff();
    const double stepLimit = i_Settings._HTraceStepFraction / maxGradient;
    if (!std::isfinite(stepLimit) || !(stepLimit > 0.0) ||
        !std::isfinite(exitT) || exitT < 0.0) return fail("InvalidStep");
    const bool reachedBoundary = exitT <= stepLimit;
    const double stepLength = std::min(exitT, stepLimit);
    Eigen::Vector4d newB = b + stepLength * rate;
    if (reachedBoundary) newB[exitLocal] = 0.0;
    if (!newB.allFinite() || newB.minCoeff() < -1e-10) return fail("InvalidBarycentric");
    for (int local = 0; local < 4; ++local) newB[local] = std::max(0.0, newB[local]);
    const double weightSum = newB.sum();
    if (!(weightSum > 0.0)) return fail("InvalidBarycentric");
    newB /= weightSum;
    if (!reachedBoundary && (newB - b).squaredNorm() == 0.0) return fail("NoProgress");

    // Record every accepted substep, including the final boundary crossing.
    // Zero-distance transitions do not create zero-length graph edges.
    if (outPath) {
      Eigen::Vector3d crossing = Eigen::Vector3d::Zero();
      for (int local = 0; local < 4; ++local)
        crossing += newB[local] * toEigenVec(_Tet.Vertices[tet[local]].Position);
      if (crossing.allFinite() &&
          (crossing - toEigenVec(outPath->back())).squaredNorm() > 0.0)
        outPath->emplace_back(crossing.x(), crossing.y(), crossing.z());
    }

    if (!reachedBoundary) {
      b = newB;
      continue;
    }

    std::array<int, 3> face;
    int w = 0;
    for (int local = 0; local < 4; ++local) if (local != exitLocal) face[w++] = tet[local];
    std::array<int, 3> faceKey = face;
    std::sort(faceKey.begin(), faceKey.end());

    const auto found = faceAdjacency.find(faceKey);
    const std::array<int, 2> slots = (found != faceAdjacency.end()) ? found->second : std::array<int, 2>{-1, -1};
    const int neighborTet = (slots[0] == tetIndex) ? slots[1] : slots[0];

    if (neighborTet < 0) {
      const bool allBottom = _Tet.vertextag[face[0]] == static_cast<int>(VertexTag::Bottom) &&
                             _Tet.vertextag[face[1]] == static_cast<int>(VertexTag::Bottom) &&
                             _Tet.vertextag[face[2]] == static_cast<int>(VertexTag::Bottom);
      if (!allBottom) return fail("NonBottomBoundary");

      double u = 0.0, v = 0.0;
      for (int local = 0; local < 4; ++local) {
        if (local == exitLocal) continue;
        const zPoint &uv = _BottomUV[tet[local]];
        u += newB[local] * uv.x;
        v += newB[local] * uv.y;
      }
      outUV = zPoint(u, v, 0.0);
      return toEigenVec(outUV).allFinite() ? true : fail("InvalidBottomUV");
    }

    const std::array<int, 4> &neighborVerts = _Tet.Volumes[neighborTet].VertexIds;
    Eigen::Vector4d neighborB = Eigen::Vector4d::Zero();
    for (int nLocal = 0; nLocal < 4; ++nLocal) {
      for (int local = 0; local < 4; ++local) {
        if (local == exitLocal) continue;
        if (neighborVerts[nLocal] == tet[local]) neighborB[nLocal] = newB[local];
      }
    }

    tetIndex = neighborTet;
    tet = neighborVerts;
    b = neighborB;
  }
  return fail("StepLimit");
}

// Tangential P1 gradients of a surface triangle embedded in 3D, via the
// standard "rotate the opposite edge by 90 degrees" formula
// grad(b_i) = (normal x e_i) / (2*Area), e_i = the edge opposite vertex i.
// This is exactly the paper's rot() operator (Section IV.D) - it gives the
// in-plane gradient of a linear function on the triangle without leaving
// the surface.
ZSPACE_INLINE bool HarmonicSlicer::computeTriElement(const std::vector<zTetrahedron::TetVertex> &vertices, const std::array<int, 3> &tri, zTetrahedron::TetFace &element, Eigen::Matrix3d &gradients) {
  const Eigen::Vector3d p0 = toEigenVec(vertices[tri[0]].Position);
  const Eigen::Vector3d p1 = toEigenVec(vertices[tri[1]].Position);
  const Eigen::Vector3d p2 = toEigenVec(vertices[tri[2]].Position);

  const Eigen::Vector3d e0 = p2 - p1; // opposite vertex 0
  const Eigen::Vector3d e1 = p0 - p2; // opposite vertex 1
  const Eigen::Vector3d e2 = p1 - p0; // opposite vertex 2

  const Eigen::Vector3d normalRaw = e2.cross(-e1); // (p1-p0) x (p2-p0)
  const double doubleArea = normalRaw.norm();
  element.Valid = doubleArea >= 1e-8;
  if (!element.Valid) return false;

  element.Area = 0.5 * doubleArea;
  const Eigen::Vector3d normal = normalRaw / doubleArea;
  gradients.col(0) = normal.cross(e0) / doubleArea;
  gradients.col(1) = normal.cross(e1) / doubleArea;
  gradients.col(2) = normal.cross(e2) / doubleArea;
  return true;
}

// Walks triangle-to-triangle across the wall subset of _Tet.Surface
// from vertexId's own barycentric one-hot coordinate, stepping along
// -grad(phi) restricted to the surface (nodalGradient here is the surface's
// own area-weighted recovered gradient - see computeVolumetricParameterization()).
// This mirrors traceToBottom()'s tet-walking exactly, one dimension down: a
// triangle's barycentric coordinates are still affine in the true (tangential)
// element gradients, so the same "first coordinate to hit zero" rule finds
// the edge the walk exits through, and crossing into the neighbor triangle
// just carries the shared-edge weights across.
//
// The wall subset has open cap seams. Successful tracing terminates
// only when both endpoints of the crossed edge are Bottom-tagged,
// interpolating their _BottomUV at the crossing. An edge with just one
// Bottom endpoint is traversed normally instead of terminating. This keeps
// the whole trace on the boundary surface itself - it never steps into an
// interior tet - which is exactly the paper's wall-tracing case (Section
// IV.D: restrict the harmonic function to the wall, trace along its own
// gradient) instead of the general 3D volume case.
ZSPACE_INLINE bool HarmonicSlicer::traceWallToBottom(int vertexId, const WallTraceData &data,
    zPoint &outUV, zPointArray *outPath, std::string *failureReason) const {
  const auto &vertexTri = data.Seeds;
  const auto &edgeAdjacency = data.Adjacency;
  const auto &nodalGradient = data.Gradient;
  const auto &vertexWallTriangles = data.IncidentTriangles;
  if (failureReason) failureReason->clear();
  const auto fail = [&](const char *reason) {
    if (failureReason) *failureReason = reason;
    return false;
  };
  if (outPath) {
    outPath->clear();
    outPath->push_back(_Tet.Surface.Positions[vertexId]);
  }
  int triIndex = vertexTri[vertexId];
  if (triIndex < 0) return fail("NoStartTriangle");

  std::array<int, 3> tri = _Tet.Faces[_Tet.SurfaceFaceIds[triIndex]].VertexIds;
  Eigen::Vector3d b = Eigen::Vector3d::Zero();
  for (int local = 0; local < 3; ++local) if (tri[local] == vertexId) b[local] = 1.0;

  int valleyTarget = -1;
  bool resolveVertex = false;
  const int maxSteps = i_Settings._HTraceMaxSteps;
  for (int step = 0; step < maxSteps; ++step) {
    if (resolveVertex) {
      int atLocal = 0;
      b.maxCoeff(&atLocal);
      const int at = tri[atLocal];
      // Search the complete incident wall fan, not just the previous
      // neighbor - but via the precomputed vertex->wall-triangle adjacency
      // (every triangle incident to `at`), not a scan of every wall edge:
      // an edge incident to `at` is by definition carried by a triangle
      // that has `at` as one of its own 3 vertices, so this visits exactly
      // the same edges without touching unrelated ones.
      const std::vector<int> &incident = vertexWallTriangles[at];
      int bestEdgeTri = -1, bestEdgeTarget = -1;
      double bestSlope = 0.0;
      std::set<int> seenOther;
      for (int t : incident) {
        const auto &verts = _Tet.Faces[_Tet.SurfaceFaceIds[t]].VertexIds;
        for (int j = 0; j < 3; ++j) {
          const int other = verts[j];
          if (other == at || !seenOther.insert(other).second) continue;
          // A Bottom vertex is accepted only on an actual Bottom-Bottom wall edge.
          if (_Tet.vertextag[at] == static_cast<int>(VertexTag::Bottom) &&
              _Tet.vertextag[other] == static_cast<int>(VertexTag::Bottom)) {
            outUV = _BottomUV[at];
            return toEigenVec(outUV).allFinite() ? true : fail("InvalidBottomUV");
          }
          const double length = (toEigenVec(_Tet.Surface.Positions[other]) -
                                 toEigenVec(_Tet.Surface.Positions[at])).norm();
          if (!(length > 0.0)) continue;
          const double slope = (static_cast<double>(_Tet.UVW[at].z) - _Tet.UVW[other].z) / length;
          if (std::isfinite(slope) && slope > bestSlope) {
            bestSlope = slope;
            bestEdgeTarget = other;
            bestEdgeTri = t;
          }
        }
      }
      int chosen = -1;
      double bestMargin = 0.0;
      for (int t : incident) {
        const zTetrahedron::TetFace &candidate = _Tet.Faces[_Tet.SurfaceFaceIds[t]];
    const auto &candidateGradients = _Tet.FaceBasisGradients[_Tet.SurfaceFaceIds[t]];
        if (!candidate.Valid) continue;
        const auto &verts = _Tet.Faces[_Tet.SurfaceFaceIds[t]].VertexIds;
        const Eigen::Vector3d normal = candidateGradients.col(0).cross(candidateGradients.col(1)).normalized();
        Eigen::Vector3d d = -nodalGradient[at];
        d -= normal * d.dot(normal);
        if (!d.allFinite() || d.norm() < 1e-15) continue;
        d.normalize();
        Eigen::Vector3d values(_Tet.UVW[verts[0]].z, _Tet.UVW[verts[1]].z, _Tet.UVW[verts[2]].z);
        if ((candidateGradients * values).dot(d) >= 0.0) continue;
        double margin = 1.0;
        for (int j = 0; j < 3; ++j) if (verts[j] != at)
          margin = std::min(margin, candidateGradients.col(j).dot(d) / candidateGradients.col(j).norm());
        if (margin > bestMargin) { bestMargin = margin; chosen = t; }
      }
      if (chosen < 0) {
        if (bestEdgeTri < 0) return fail("ValleyStalled");
        chosen = bestEdgeTri;
        valleyTarget = bestEdgeTarget;
      }
      triIndex = chosen;
      tri = _Tet.Faces[_Tet.SurfaceFaceIds[chosen]].VertexIds;
      b.setZero();
      for (int j = 0; j < 3; ++j) if (tri[j] == at) b[j] = 1.0;
      resolveVertex = false;
    }

    const zTetrahedron::TetFace &element = _Tet.Faces[_Tet.SurfaceFaceIds[triIndex]];
    const auto &gradients = _Tet.FaceBasisGradients[_Tet.SurfaceFaceIds[triIndex]];
    if (!element.Valid) return fail("DegenerateElement");

    // In edge mode advance along the linear edge field in bounded substeps.
    if (valleyTarget >= 0) {
      Eigen::Vector3d x = Eigen::Vector3d::Zero();
      int targetLocal = -1;
      for (int j = 0; j < 3; ++j) {
        x += b[j] * toEigenVec(_Tet.Surface.Positions[tri[j]]);
        if (tri[j] == valleyTarget) targetLocal = j;
      }
      if (targetLocal < 0) return fail("ValleyInvalidEdge");
      const Eigen::Vector3d target = toEigenVec(_Tet.Surface.Positions[valleyTarget]);
      const double distance = (target - x).norm();
      const double limit = i_Settings._HTraceStepFraction / gradients.colwise().norm().maxCoeff();
      if (!std::isfinite(distance) || !std::isfinite(limit) || !(limit > 0.0)) return fail("InvalidStep");
      const double fraction = distance <= limit ? 1.0 : limit / distance;
      const Eigen::Vector3d oldB = b;
      b *= 1.0 - fraction;
      b[targetLocal] += fraction;
      if (fraction < 1.0 && (b - oldB).squaredNorm() == 0.0) return fail("NoProgress");
      const Eigen::Vector3d point = (1.0 - fraction) * x + fraction * target;
      if (outPath && (point - toEigenVec(outPath->back())).squaredNorm() > 0.0)
        outPath->emplace_back(point.x(), point.y(), point.z());
      if (fraction == 1.0) { valleyTarget = -1; resolveVertex = true; }
      continue;
    }

    Eigen::Vector3d gradAtPoint = b[0] * nodalGradient[tri[0]] + b[1] * nodalGradient[tri[1]] + b[2] * nodalGradient[tri[2]];
    const Eigen::Vector3d normal = gradients.col(0).cross(gradients.col(1)).normalized();
    gradAtPoint -= normal * gradAtPoint.dot(normal);
    const double gradNorm = gradAtPoint.norm();
    if (!gradAtPoint.allFinite() || !std::isfinite(gradNorm)) return fail("NonFiniteGradient");
    if (gradNorm < 1e-15) return fail("ZeroGradient");
    const Eigen::Vector3d direction = -gradAtPoint / gradNorm;

    Eigen::Vector3d rate;
    for (int local = 0; local < 3; ++local) rate[local] = gradients.col(local).dot(direction);

    int exitLocal = -1;
    double exitT = std::numeric_limits<double>::max();
    for (int local = 0; local < 3; ++local) {
      if (rate[local] >= -1e-15) continue;
      const double t = b[local] / -rate[local];
      if (t < exitT) {
        exitT = t;
        exitLocal = local;
      }
    }
    if (exitLocal < 0) {
      // No tangential direction decreases any coordinate of this triangle -
      // exactly what happens when the current position is a mesh vertex
      // sitting on a fold/ridge: the recovered -grad(phi) there lands
      // almost exactly along the fold edge, which neither adjacent triangle
      // sees as an interior-entry direction. If we're sitting at a vertex
      // (b is one-hot), fall into the same apex fan-search/valley-following
      // resolveVertex() uses mid-walk instead of failing outright - it scans
      // every wall triangle actually incident to this vertex, not just the
      // one we happened to be in.
      if (b.maxCoeff() >= 1.0 - 1e-9) {
        resolveVertex = true;
        continue;
      }
      return fail("NoExit");
    }

    // |grad(b_i)| is the reciprocal altitude opposite vertex i.
    // Limit each step locally, then re-evaluate the direction even when
    // the point remains inside the same simplex. Never step across a face.
    const double maxGradient = gradients.colwise().norm().maxCoeff();
    const double stepLimit = i_Settings._HTraceStepFraction / maxGradient;
    if (!std::isfinite(stepLimit) || !(stepLimit > 0.0) ||
        !std::isfinite(exitT) || exitT < 0.0) return fail("InvalidStep");
    const bool reachedBoundary = exitT <= stepLimit;
    const double stepLength = std::min(exitT, stepLimit);
    Eigen::Vector3d newB = b + stepLength * rate;
    if (reachedBoundary) newB[exitLocal] = 0.0;
    if (!newB.allFinite() || newB.minCoeff() < -1e-10) return fail("InvalidBarycentric");
    for (int local = 0; local < 3; ++local) newB[local] = std::max(0.0, newB[local]);
    const double weightSum = newB.sum();
    if (!(weightSum > 0.0)) return fail("InvalidBarycentric");
    newB /= weightSum;
    if (!reachedBoundary && (newB - b).squaredNorm() == 0.0) return fail("NoProgress");

    if (outPath) {
      Eigen::Vector3d crossing = Eigen::Vector3d::Zero();
      for (int local = 0; local < 3; ++local)
        crossing += newB[local] * toEigenVec(_Tet.Surface.Positions[tri[local]]);
      if (crossing.allFinite() &&
          (crossing - toEigenVec(outPath->back())).squaredNorm() > 0.0)
        outPath->emplace_back(crossing.x(), crossing.y(), crossing.z());
    }

    if (!reachedBoundary) {
      b = newB;
      continue;
    }

    const int localA = (exitLocal + 1) % 3;
    const int localB = (exitLocal + 2) % 3;
    const int edgeA = tri[localA];
    const int edgeB = tri[localB];

    const bool aBottom = _Tet.vertextag[edgeA] == static_cast<int>(VertexTag::Bottom);
    const bool bBottom = _Tet.vertextag[edgeB] == static_cast<int>(VertexTag::Bottom);
    if (aBottom && bBottom) {
      const zPoint &uvA = _BottomUV[edgeA];
      const zPoint &uvB = _BottomUV[edgeB];
      const double sum = newB[localA] + newB[localB];
      if (!(sum > 1e-12)) return fail("InvalidBottomWeights");
      const double wa = newB[localA] / sum;
      const double wb = 1.0 - wa;
      outUV = zPoint(uvA.x * wa + uvB.x * wb, uvA.y * wa + uvB.y * wb, 0.0);
      return toEigenVec(outUV).allFinite() ? true : fail("InvalidBottomUV");
    }

    const unsigned long long edgeKey = makeEdgeKey(edgeA, edgeB);
    const auto found = edgeAdjacency.find(edgeKey);
    const std::array<int, 2> slots = (found != edgeAdjacency.end()) ? found->second : std::array<int, 2>{-1, -1};
    const int neighborTri = (slots[0] == triIndex) ? slots[1] : slots[0];
    if (neighborTri < 0) return fail("NonBottomSeamOrMissingNeighbor"); // left the wall through a non-bottom seam or missing adjacency

    const std::array<int, 3> &neighborVerts = _Tet.Faces[_Tet.SurfaceFaceIds[neighborTri]].VertexIds;
    Eigen::Vector3d neighborB = Eigen::Vector3d::Zero();
    for (int nLocal = 0; nLocal < 3; ++nLocal) {
      if (neighborVerts[nLocal] == edgeA) neighborB[nLocal] = newB[localA];
      else if (neighborVerts[nLocal] == edgeB) neighborB[nLocal] = newB[localB];
    }

    // Test the neighbor direction before crossing. If both sides drive
    // toward the edge, follow the edge's own decreasing scalar field.
    const zTetrahedron::TetFace &neighborElement = _Tet.Faces[_Tet.SurfaceFaceIds[neighborTri]];
    const auto &neighborGradients = _Tet.FaceBasisGradients[_Tet.SurfaceFaceIds[neighborTri]];
    if (!neighborElement.Valid) return fail("DegenerateElement");
    Eigen::Vector3d ng = Eigen::Vector3d::Zero();
    int opposite = -1;
    for (int j = 0; j < 3; ++j) {
      ng += neighborB[j] * nodalGradient[neighborVerts[j]];
      if (neighborVerts[j] != edgeA && neighborVerts[j] != edgeB) opposite = j;
    }
    if (!ng.allFinite() || opposite < 0) return fail("NonFiniteGradient");
    const double entryRate = neighborGradients.col(opposite).dot(-ng);
    if (entryRate < 0.0) {
      Eigen::Vector3d currentGradient = Eigen::Vector3d::Zero();
      for (int j = 0; j < 3; ++j) currentGradient += newB[j] * nodalGradient[tri[j]];
      const double currentInward = gradients.col(exitLocal).dot(-currentGradient);
      b = newB;
      // A finite step may have overshot a turn: stay in this face if the
      // refreshed direction already enters it. Only opposed directions slide.
      if (currentInward > 0.0) continue;
      if (b.maxCoeff() >= 1.0 - 1e-12) {
        int j = 0; b.maxCoeff(&j); b.setZero(); b[j] = 1.0;
        resolveVertex = true;
        continue;
      }
      const double edgePhi = newB[localA] * _Tet.UVW[edgeA].z + newB[localB] * _Tet.UVW[edgeB].z;
      const int lower = _Tet.UVW[edgeA].z < _Tet.UVW[edgeB].z ? edgeA : edgeB;
      if (!std::isfinite(edgePhi) || !(static_cast<double>(_Tet.UVW[lower].z) < edgePhi))
        return fail("ValleyFlatEdge");
      valleyTarget = lower;
      continue;
    }

    triIndex = neighborTri;
    tri = neighborVerts;
    b = neighborB;
  }
  return fail("StepLimit");
}

ZSPACE_INLINE bool HarmonicSlicer::computeVolumetricParameterization(bool useParallel) {
  zTsUtils::ScopedSlicingTiming methodTiming("HarmonicSlicer.computeVolumetricParameterization");
  const bool success = computeVolumetricParameterizationImpl(useParallel);
  if (success) zTsUtils::printSuccess("[HarmonicSlicer::computeVolumetricParameterization] Success");
  else zTsUtils::printError("[HarmonicSlicer::computeVolumetricParameterization] Failed");
  return success;
}

ZSPACE_INLINE bool HarmonicSlicer::computeVolumetricParameterizationImpl(bool useParallel) {
  clearParameterization();
  if (!_Tet.FieldValid || _Tet.UVW.empty() || _Tet.UVW.size() != _Tet.Vertices.size()) {
    _Status = "computeVolumetricParameterization failed: no harmonic field. Call computeHarmonicField() first.";
    printDebug("[computeVolumetricParameterizationImpl] " + _Status);
    return false;
  }

  const float invalidUV = std::numeric_limits<float>::quiet_NaN();
  if (!std::isfinite(i_Settings._HTraceStepFraction) ||
      i_Settings._HTraceStepFraction <= 0.0 || i_Settings._HTraceStepFraction > 1.0 ||
      i_Settings._HTraceMaxSteps <= 0) {
    _Status = "computeVolumetricParameterization failed: _HTraceStepFraction must be in (0,1] and _HTraceMaxSteps must be positive.";
    printDebug("[computeVolumetricParameterizationImpl] " + _Status);
    return false;
  }
  if (!buildBottomParameterization()) {
    _Status = "computeVolumetricParameterization failed: could not build a finite bottom chart; failed="
        + std::to_string(_Tet.Vertices.size()) + ".";
    printDebug("[computeVolumetricParameterizationImpl] " + _Status);
    return false;
  }

  // Seed ids are selected after recovering gradients, not by mesh order.
  VolumeTraceData volume;
  auto &vertexTet = volume.Seeds;
  vertexTet.assign(_Tet.Vertices.size(), -1);

  // Volume-weighted average of grad(phi) over every tet incident to each
  // vertex - a continuous (piecewise-linear) recovered gradient field, used
  // instead of the raw per-tet constant gradient so traceToBottom()'s walk
  // direction doesn't depend on which single incident tet happened to seed
  // or continue it (see the comment on traceToBottom()).
  auto &nodalGradient = volume.Gradient;
  nodalGradient.assign(_Tet.Vertices.size(), Eigen::Vector3d::Zero());
  std::vector<double> nodalWeight(_Tet.Vertices.size(), 0.0);
  for (size_t t = 0; t < _Tet.Volumes.size(); ++t) {
    const zTetrahedron::TetVolume &element = _Tet.Volumes[t];
    const auto &gradients = _Tet.VolumeBasisGradients[t];
    if (!element.Valid) continue;
    const std::array<int, 4> &tet = _Tet.Volumes[t].VertexIds;
    const Eigen::Vector3d &gradPhi = _Tet.VolumeFieldGradients[t];
    for (int local = 0; local < 4; ++local) {
      nodalGradient[tet[local]] += element.Volume * gradPhi;
      nodalWeight[tet[local]] += element.Volume;
    }
  }
  for (size_t i = 0; i < nodalGradient.size(); ++i) {
    if (nodalWeight[i] > 1e-5) nodalGradient[i] /= nodalWeight[i];
  }

  // Test the tangent cone at each incident vertex. Normalize rates by their
  // maximum magnitude so the margin is independent of model scale. A positive
  // minimum margin prefers strict entry over directions on faces/edges.
  // Negative rates are conservatively rejected, even if tiny: accepting them
  // would immediately create a zero-distance outward step in the walker.
  const auto entryMargin = [](const auto &gradients, int start, const Eigen::Vector3d &direction) {
    const auto rates = (gradients.transpose() * direction).eval();
    if (!rates.allFinite()) return -1.0;
    const double scale = rates.cwiseAbs().maxCoeff();
    if (!(scale > 0.0) || !(rates[start] < 0.0)) return -1.0;
    double margin = 1.0;
    for (int j = 0; j < rates.size(); ++j) {
      if (j == start) continue;
      if (rates[j] < 0.0) return -1.0;
      margin = std::min(margin, rates[j] / scale);
    }
    return margin;
  };
  std::vector<double> tetMargin(_Tet.Vertices.size(), -1.0);
  for (size_t t = 0; t < _Tet.Volumes.size(); ++t) {
    const auto &element = _Tet.Volumes[t];
    const auto &gradients = _Tet.VolumeBasisGradients[t];
    if (!element.Valid) continue;
    for (int local = 0; local < 4; ++local) {
      const int id = _Tet.Volumes[t].VertexIds[local];
      const auto &g = nodalGradient[id];
      if (!g.allFinite() || g.norm() < 1e-15) continue;
      const Eigen::Vector3d direction = -g.normalized();
      const double margin = entryMargin(gradients, local, direction);
      if (margin > tetMargin[id]) {
        tetMargin[id] = margin;
        vertexTet[id] = static_cast<int>(t);
      }
    }
  }

  // Sorted-face-key -> the one or two tets sharing it, so traceToBottom() can
  // walk tet-to-tet across shared faces.
  auto &faceAdjacency = volume.Adjacency;
  for (const auto &face : _Tet.Faces) {
    auto key = face.VertexIds;
    std::sort(key.begin(), key.end());
    faceAdjacency.emplace(key, face.VolumeIds);
  }

  // Build adjacency and recovered gradients on wall triangles only. A cap
  // triangle is classified by its source facet marker, independently of vertex tags.
  // Top vertices incident to a wall triangle must also use surface tracing.
  // Surface vertex ids retain the existing TetGen input-index correspondence.
  const size_t surfaceVertexCount = _Tet.Surface.Positions.size();
  WallTraceData wall;
  auto &vertexTri = wall.Seeds;
  vertexTri.assign(surfaceVertexCount, -1);
  std::vector<bool> onWall(surfaceVertexCount, false);
  std::vector<int> wallTriangles;
  auto &edgeAdjacency = wall.Adjacency;
  // Every wall triangle incident to each vertex, built once alongside
  // edgeAdjacency - traceWallToBottom()'s resolveVertex() ridge handling
  // uses this instead of scanning the entire edgeAdjacency map per call.
  auto &vertexWallTriangles = wall.IncidentTriangles;
  vertexWallTriangles.resize(surfaceVertexCount);
  // Surface indices map to shared faces; vertices and basis columns use face order.
  auto &surfaceNodalGradient = wall.Gradient;
  surfaceNodalGradient.assign(surfaceVertexCount, Eigen::Vector3d::Zero());
  std::vector<double> surfaceNodalWeight(surfaceVertexCount, 0.0);
  for (size_t t = 0; t < _Tet.Surface.Triangles.size(); ++t) {
    const std::array<int, 3> &tri = _Tet.Faces[_Tet.SurfaceFaceIds[t]].VertexIds;
    const auto &face = _Tet.Faces[_Tet.SurfaceFaceIds[t]];
    if (_Tet.FaceTags[_Tet.SurfaceFaceIds[t]] != HarmonicField::TetFaceTag::Side) continue;
    wallTriangles.push_back(static_cast<int>(t));
    for (int local = 0; local < 3; ++local) {
      onWall[tri[local]] = true;
      vertexWallTriangles[tri[local]].push_back(static_cast<int>(t));
      const unsigned long long key = makeEdgeKey(tri[local], tri[(local + 1) % 3]);
      auto &slots = edgeAdjacency.try_emplace(key, std::array<int, 2>{-1, -1}).first->second;
      if (slots[0] < 0) slots[0] = static_cast<int>(t);
      else slots[1] = static_cast<int>(t);
    }

    const auto &element = face;
    if (!face.Valid) continue;
    const Eigen::Vector3d &gradPhi = _Tet.FaceFieldGradients[_Tet.SurfaceFaceIds[t]];
    for (int local = 0; local < 3; ++local) {
      surfaceNodalGradient[tri[local]] += element.Area * gradPhi;
      surfaceNodalWeight[tri[local]] += element.Area;
    }
  }
  for (size_t i = 0; i < surfaceNodalGradient.size(); ++i) {
    if (surfaceNodalWeight[i] > 1e-8) surfaceNodalGradient[i] /= surfaceNodalWeight[i];
  }

  std::vector<double> triMargin(surfaceVertexCount, -1.0);
  for (int t : wallTriangles) {
    const zTetrahedron::TetFace &element = _Tet.Faces[_Tet.SurfaceFaceIds[t]];
    const auto &gradients = _Tet.FaceBasisGradients[_Tet.SurfaceFaceIds[t]];
    if (!element.Valid) continue;
    const auto &tri = _Tet.Faces[_Tet.SurfaceFaceIds[t]].VertexIds;
    const Eigen::Vector3d normal = gradients.col(0).cross(gradients.col(1)).normalized();
    for (int local = 0; local < 3; ++local) {
      const int id = tri[local];
      Eigen::Vector3d direction = -surfaceNodalGradient[id];
      direction -= normal * direction.dot(normal);
      if (!direction.allFinite() || direction.norm() < 1e-15) continue;
      direction.normalize();
      const double margin = entryMargin(gradients, local, direction);
      if (margin > triMargin[id]) {
        triMargin[id] = margin;
        vertexTri[id] = t;
      }
    }
  }

  // Fold/ridge vertices can fail the strict margin test on every incident
  // wall triangle: -grad(phi) there lands almost exactly along the fold
  // edge itself, so neither adjacent triangle sees a strict interior-entry
  // direction (the "other" coordinate dips slightly negative on both
  // sides), leaving vertexTri[id] at its -1 sentinel ("NoStartTriangle")
  // without ever reaching traceWallToBottom()'s own resolveVertex()/
  // valleyTarget ridge handling. Give every wall vertex a best-effort seed
  // triangle regardless of margin so that machinery actually gets to run;
  // this only fires where the strict pass above found nothing.
  for (int t : wallTriangles) {
    const auto &tri = _Tet.Faces[_Tet.SurfaceFaceIds[t]].VertexIds;
    for (int local = 0; local < 3; ++local) {
      const int id = tri[local];
      if (vertexTri[id] < 0) vertexTri[id] = t;
    }
  }

  // Vertex-indexed diagnostic graphs: retain failed paths only. Bottom
  // vertices and successful traces have empty graphs; validity = _Tet.UVWValid.
  o_TraceGraphs.assign(_Tet.Vertices.size(), zObjGraph());
  struct TraceResult {
    zPoint UV;
    bool Valid = false;
    bool Wall = false;
    std::string FailureReason;
    zPointArray FailedPath;
  };
  std::vector<TraceResult> results(_Tet.Vertices.size());
  // Keep the shared UVW field read-only until every trace has finished.
  zTsDetail::ParallelErrors parallelErrors1(static_cast<int>(_Tet.Vertices.size()));
#pragma omp parallel for schedule(dynamic) if(useParallel)
  for (int i = 0; i < static_cast<int>(_Tet.Vertices.size()); ++i) {
    try {
      TraceResult &result = results[i];
      result.UV = zPoint(invalidUV, invalidUV, 0.0);
      zPointArray path;
      if (_Tet.vertextag[i] == static_cast<int>(VertexTag::Bottom)) {
        result.UV = _BottomUV[i];
        result.Valid = true;
      } else {
        result.Wall = static_cast<size_t>(i) < surfaceVertexCount &&
            (_Tet.vertextag[i] == static_cast<int>(VertexTag::Side) ||
             (_Tet.vertextag[i] == static_cast<int>(VertexTag::Top) && onWall[i]));
        result.Valid = result.Wall
            ? traceWallToBottom(i, wall, result.UV, &path, &result.FailureReason)
            : traceToBottom(i, volume, result.UV, &path, &result.FailureReason);
        if (!result.Valid) result.FailedPath.swap(path);
      }
    } catch (...) { parallelErrors1.capture(i); }
  }
  parallelErrors1.rethrow();

  // Commit in vertex order: vector<bool>, statistics and graph creation
  // stay serial, with the same diagnostic ordering as before.
  std::map<std::string, std::array<size_t, 2>> failureCounts;
  size_t failed2D = 0, failed3D = 0;
  for (size_t i = 0; i < results.size(); ++i) {
    const TraceResult &result = results[i];
    _Tet.UVWValid[i] = result.Valid;
    _Tet.UVW[i].x = result.Valid ? result.UV.x : invalidUV;
    _Tet.UVW[i].y = result.Valid ? result.UV.y : invalidUV;
    if (result.Valid) continue;
    if (result.Wall) ++failed2D;
    else ++failed3D;
    ++failureCounts[result.FailureReason][result.Wall ? 0 : 1];

    zPointArray &path = results[i].FailedPath;
    if (path.size() < 2) continue;
    zIntArray edges;
    edges.reserve(2 * (path.size() - 1));
    for (size_t j = 1; j < path.size(); ++j) {
      edges.push_back(static_cast<int>(j - 1));
      edges.push_back(static_cast<int>(j));
    }
    zFnGraph(o_TraceGraphs[i]).create(path, edges);
  }

  // o_IsoMeshUV/o_IsoMeshes themselves are (re)built by
  // computeContours(), which runs after this step and reads the real
  // _Tet.UVW set above directly - no placeholder/refresh pass needed.

  std::ostringstream out;
  out << "computeVolumetricParameterization: vertices=" << _Tet.UVW.size()
      << " valid=" << (_Tet.UVW.size() - failed2D - failed3D)
      << " failed=" << (failed2D + failed3D)
      << " (2D=" << failed2D << ", 3D=" << failed3D << ").";
  for (const auto &entry : failureCounts) {
    out << "\n  " << entry.first << ": 2D=" << entry.second[0]
        << " 3D=" << entry.second[1]
        << " total=" << (entry.second[0] + entry.second[1]);
  }
  _Status = out.str();
  for (const auto &entry : failureCounts) {
    if (entry.second[0] + entry.second[1] > 0)
      printDebug("[computeVolumetricParameterizationImpl] UV tracing failed: " + entry.first);
  }
  return true;
}

} // namespace zSpace
