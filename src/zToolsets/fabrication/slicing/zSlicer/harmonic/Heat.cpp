#ifndef EIGEN_DONT_PARALLELIZE
#define EIGEN_DONT_PARALLELIZE
#endif

#include <zspace/zToolsets/fabrication/slicing/zSlicer/HarmonicSlicer.h>

#ifndef ZSPACE_INLINE
#define ZSPACE_INLINE ZSPACE_TOOLSETS_INLINE
#endif

#include <Eigen/SparseCholesky>

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>

namespace zSpace {

using namespace zTsUtils;

//--------------------------
//---- METHOD 1 : HEAT GEODESICS
//--------------------------

ZSPACE_INLINE bool HarmonicSlicer::buildSurfaceOperators(const TriangleSoup &soup, SurfaceOperators &ops, int &degenerateTriangles) {
  const int n = static_cast<int>(soup.Positions.size());
  ops.Mass = Eigen::VectorXd::Zero(n);
  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(soup.Triangles.size() * 12);
  degenerateTriangles = 0;

  for (const auto &tri : soup.Triangles) {
    const Eigen::Vector3d p[3] = {toEigenVec(soup.Positions[tri[0]]), toEigenVec(soup.Positions[tri[1]]), toEigenVec(soup.Positions[tri[2]])};
    const double area2 = (p[1] - p[0]).cross(p[2] - p[0]).norm();
    if (area2 < 1e-8) {
      ++degenerateTriangles;
      continue;
    }
    for (int c = 0; c < 3; ++c) {
      ops.Mass[tri[c]] += area2 / 6.0;

      // Cotangent at corner c weights the opposite edge (c+1, c+2).
      const int i = tri[(c + 1) % 3];
      const int j = tri[(c + 2) % 3];
      const double cot = (p[(c + 1) % 3] - p[c]).dot(p[(c + 2) % 3] - p[c]) / area2;
      const double w = 0.5 * cot;
      triplets.emplace_back(i, j, -w);
      triplets.emplace_back(j, i, -w);
      triplets.emplace_back(i, i, w);
      triplets.emplace_back(j, j, w);
    }
  }

  ops.MeanEdgeLength = meanEdgeLength(soup.Positions, soup.Triangles);
  if (!(ops.MeanEdgeLength > 0.0)) return false;
  ops.Laplacian.resize(n, n);
  ops.Laplacian.setFromTriplets(triplets.begin(), triplets.end());
  return true;
}

// Heat method (Crane, Weischedel, Wardetzky): diffuse from the sources,
// normalise the negated gradient, then recover distance from its divergence.
ZSPACE_INLINE bool HarmonicSlicer::computeHeatGeodesics(const TriangleSoup &soup, const SurfaceOperators &ops,
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> &heatSolver,
    const zIntArray &sourceIds, zFloatArray &outDistance, std::string &error) {
  const int n = static_cast<int>(soup.Positions.size());
  if (n == 0 || sourceIds.empty()) {
    error = "empty mesh or source set";
    return false;
  }

  Eigen::VectorXd delta = Eigen::VectorXd::Zero(n);
  for (int id : sourceIds) delta[id] = 1.0;

  const Eigen::VectorXd u = heatSolver.solve(delta);
  if (heatSolver.info() != Eigen::Success || !u.allFinite()) {
    error = "heat solve failed";
    return false;
  }

  Eigen::VectorXd divergence = Eigen::VectorXd::Zero(n);
  for (const auto &tri : soup.Triangles) {
    const Eigen::Vector3d p[3] = {toEigenVec(soup.Positions[tri[0]]), toEigenVec(soup.Positions[tri[1]]), toEigenVec(soup.Positions[tri[2]])};
    const Eigen::Vector3d normalScaled = (p[1] - p[0]).cross(p[2] - p[0]);
    const double area2 = normalScaled.norm();
    if (area2 < 1e-8) continue;
    const Eigen::Vector3d normal = normalScaled / area2;

    // grad u = (1 / 2A) * sum_i u_i (N x e_i), e_i opposite vertex i.
    Eigen::Vector3d gradient = Eigen::Vector3d::Zero();
    for (int c = 0; c < 3; ++c) {
      const Eigen::Vector3d opposite = p[(c + 2) % 3] - p[(c + 1) % 3];
      gradient += u[tri[c]] * normal.cross(opposite);
    }
    gradient /= area2;
    const double gradientNorm = gradient.norm();
    if (gradientNorm < 1e-14) continue;
    const Eigen::Vector3d X = -gradient / gradientNorm;

    for (int c = 0; c < 3; ++c) {
      const Eigen::Vector3d &pi = p[c];
      const Eigen::Vector3d &pj = p[(c + 1) % 3];
      const Eigen::Vector3d &pk = p[(c + 2) % 3];
      const Eigen::Vector3d e1 = pj - pi;
      const Eigen::Vector3d e2 = pk - pi;
      const double cotAtK = (pi - pk).dot(pj - pk) / area2;
      const double cotAtJ = (pi - pj).dot(pk - pj) / area2;
      divergence[tri[c]] += 0.5 * (cotAtK * e1.dot(X) + cotAtJ * e2.dot(X));
    }
  }

  // The cotangent matrix here is the positive semi-definite weak form of
  // -Laplace, so Laplace(phi) = div X becomes L phi = -div. The source
  // loop is a Dirichlet boundary (phi = 0) so every source vertex sits
  // exactly on the zero level, matching the harmonic cap constraints.
  std::vector<int> freeIndex(n, -1);
  int freeCount = 0;
  for (int id : sourceIds) freeIndex[id] = -2;
  for (int i = 0; i < n; ++i) {
    if (freeIndex[i] == -1) freeIndex[i] = freeCount++;
  }

  Eigen::VectorXd phi = Eigen::VectorXd::Zero(n);
  if (freeCount > 0) {
    std::vector<Eigen::Triplet<double>> triplets;
    Eigen::VectorXd rhs(freeCount);
    for (int k = 0; k < ops.Laplacian.outerSize(); ++k) {
      for (Eigen::SparseMatrix<double>::InnerIterator it(ops.Laplacian, k); it; ++it) {
        const int row = freeIndex[static_cast<int>(it.row())];
        const int column = freeIndex[static_cast<int>(it.col())];
        if (row >= 0 && column >= 0) triplets.emplace_back(row, column, it.value());
      }
    }
    for (int i = 0; i < n; ++i) {
      if (freeIndex[i] >= 0) rhs[freeIndex[i]] = -divergence[i];
    }

    Eigen::SparseMatrix<double> poissonMatrix(freeCount, freeCount);
    poissonMatrix.setFromTriplets(triplets.begin(), triplets.end());
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> poissonSolver;
    poissonSolver.compute(poissonMatrix);
    if (poissonSolver.info() != Eigen::Success) {
      error = "poisson factorization failed";
      return false;
    }
    const Eigen::VectorXd solution = poissonSolver.solve(rhs);
    if (poissonSolver.info() != Eigen::Success || !solution.allFinite()) {
      error = "poisson solve failed";
      return false;
    }
    for (int i = 0; i < n; ++i) {
      if (freeIndex[i] >= 0) phi[i] = solution[freeIndex[i]];
    }
  }

  if (phi.mean() < 0.0) phi = -phi;

  outDistance.assign(n, 0.0f);
  for (int i = 0; i < n; ++i) outDistance[i] = static_cast<float>(std::max(0.0, phi[i]));
  return true;
}

// d_b + d_t at a vertex is the shortest Bottom->Top surface path through
// it; over the side (wall) vertices its minimum is the shortest wall
// geodesic and its maximum the longest. Cap vertices are left out: a path
// across a cap is not a printed column. False before computeHeatDistances().
ZSPACE_INLINE bool HarmonicSlicer::columnLengths(double &outShortest, double &outLongest) const {
  const zFloatArray &toBottom = _Tet.HeatValue.first;
  const zFloatArray &toTop = _Tet.HeatValue.second;
  const size_t n = _Tet.Surface.Positions.size();
  if (n == 0 || toBottom.size() != n || toTop.size() != n || _Tet.vertextag.size() < n) return false;
  double shortest = (std::numeric_limits<double>::max)(), longest = 0.0;
  int sideCount = 0;
  for (size_t i = 0; i < n; ++i) {
    const int tag = _Tet.vertextag[i];
    if (tag == static_cast<int>(VertexTag::Bottom) || tag == static_cast<int>(VertexTag::Top)) continue;
    const double length = static_cast<double>(toBottom[i]) + static_cast<double>(toTop[i]);
    shortest = (std::min)(shortest, length);
    longest = (std::max)(longest, length);
    ++sideCount;
  }
  if (sideCount == 0 || shortest <= 0.0) return false;
  outShortest = shortest;
  outLongest = longest;
  return true;
}

ZSPACE_INLINE bool HarmonicSlicer::computeHeatDistances() {
  zTsUtils::ScopedSlicingTiming methodTiming("HarmonicSlicer.computeHeatDistances");
  const bool success = computeHeatDistancesImpl();
  if (success) zTsUtils::printSuccess("[HarmonicSlicer::computeHeatDistances] Success");
  else zTsUtils::printError("[HarmonicSlicer::computeHeatDistances] Failed");
  return success;
}

ZSPACE_INLINE bool HarmonicSlicer::computeHeatDistancesImpl() {
  clearContours();
  clearHeat();
  const size_t n = _Tet.Surface.Positions.size();
  if (_Tet.Surface.Triangles.empty() || _Tet.vertextag.size() < n) {
    _Status = "computeHeatDistances failed: no tagged surface. Call computeBasicStructure() first.";
    printDebug("[computeHeatDistancesImpl] " + _Status);
    return false;
  }

  // Surface vertex i is tet vertex i, so _Tet.vertextag tags it.
  zIntArray bottomIds;
  zIntArray topIds;
  for (int i = 0; i < static_cast<int>(n); ++i) {
    if (_Tet.vertextag[i] == static_cast<int>(VertexTag::Bottom)) bottomIds.push_back(i);
    else if (_Tet.vertextag[i] == static_cast<int>(VertexTag::Top)) topIds.push_back(i);
  }

  SurfaceOperators ops;
  int degenerateTriangles = 0;
  if (!buildSurfaceOperators(_Tet.Surface, ops, degenerateTriangles)) {
    _Status = "computeHeatDistances failed: could not build surface operators.";
    printDebug("[computeHeatDistancesImpl] " + _Status);
    return false;
  }

  // Both source sets use the same diffusion matrix; only the Poisson
  // Dirichlet systems below depend on the chosen cap.
  Eigen::SparseMatrix<double> heatMatrix = ops.Laplacian * (ops.MeanEdgeLength * ops.MeanEdgeLength);
  for (int i = 0; i < static_cast<int>(n); ++i) heatMatrix.coeffRef(i, i) += ops.Mass[i];
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> heatSolver;
  heatSolver.compute(heatMatrix);
  if (heatSolver.info() != Eigen::Success) {
    _Status = "computeHeatDistances failed: heat factorization failed";
    printDebug("[computeHeatDistancesImpl] " + _Status);
    return false;
  }
  std::string error;
  if (!computeHeatGeodesics(_Tet.Surface, ops, heatSolver, bottomIds, _Tet.HeatValue.first, error)) {
    _Status = "computeHeatDistances failed (bottom): " + error;
    printDebug("[computeHeatDistancesImpl] " + _Status);
    return false;
  }
  if (!computeHeatGeodesics(_Tet.Surface, ops, heatSolver, topIds, _Tet.HeatValue.second, error)) {
    clearHeat();
    _Status = "computeHeatDistances failed (top): " + error;
    printDebug("[computeHeatDistancesImpl] " + _Status);
    return false;
  }

  double shortest = 0.0, longest = 0.0;
  if (!columnLengths(shortest, longest)) {
    clearHeat();
    _Status = "computeHeatDistances failed: no side vertex between the Bottom and Top tags to measure a "
              "Bottom->Top distance on.";
    printDebug("[computeHeatDistancesImpl] " + _Status);
    return false;
  }

  const int layers = zTsUtils::layersFromSpacing(shortest, i_Settings.LayerSpacing);
  std::ostringstream out;
  out << "computeHeatDistances: Bottom->Top wall geodesic shortest=" << shortest << " longest=" << longest;
  if (layers >= 2)
    out << " -> " << layers << " layers at LayerSpacing=" << i_Settings.LayerSpacing << ", layer height "
        << shortest / (layers - 1) << " .. " << longest / (layers - 1);
  _Status = out.str();

  return true;
}

} // namespace zSpace
