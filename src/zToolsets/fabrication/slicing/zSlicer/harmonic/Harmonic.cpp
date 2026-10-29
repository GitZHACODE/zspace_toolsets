#ifndef EIGEN_DONT_PARALLELIZE
#define EIGEN_DONT_PARALLELIZE
#endif

#include <zspace/zToolsets/fabrication/slicing/zSlicer/HarmonicSlicer.h>

#ifndef ZSPACE_INLINE
#define ZSPACE_INLINE ZSPACE_TOOLSETS_INLINE
#endif

#include <Eigen/SparseCholesky>

#include <algorithm>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace zSpace {

using namespace zTsUtils;

//--------------------------
//---- METHOD 2 : HARMONIC FIELD (w)
//--------------------------

// Dirichlet solve K x = b on the tet volume. fixedValue < 0 marks a free
// vertex. Without targetGradients this is the Laplace (harmonic) problem;
// with them b also carries sum_t vol_t grad(N_i) . F_t, the normal equations
// of min sum_t vol_t |grad x - F_t|^2 (gradient equalisation). The solved
// field is clamped to [0, 1] (bottom=0/top=1).
ZSPACE_INLINE bool HarmonicSlicer::prepareTetDirichlet(const std::vector<double> &fixedValue,
    TetDirichletSystem &system, std::string &error) const {
  const int vertexCount = static_cast<int>(_Tet.Vertices.size());
  system.FixedField.assign(vertexCount, 0.0f);
  system.FreeIndex.assign(vertexCount, -1);
  int freeCount = 0;
  for (int i = 0; i < vertexCount; ++i) {
    if (fixedValue[i] < 0.0) system.FreeIndex[i] = freeCount++;
    else system.FixedField[i] = static_cast<float>(fixedValue[i]);
  }
  system.BoundaryRhs = Eigen::VectorXd::Zero(freeCount);
  std::vector<Eigen::Triplet<double>> triplets;
  for (size_t t = 0; t < _Tet.Volumes.size(); ++t) {
    const zTetrahedron::TetVolume &element = _Tet.Volumes[t];
    if (!element.Valid) continue;
    const auto &tet = _Tet.Volumes[t].VertexIds;
    for (int i = 0; i < 4; ++i) {
      const int row = system.FreeIndex[tet[i]];
      if (row < 0) continue;
      for (int j = 0; j < 4; ++j) {
        const int column = system.FreeIndex[tet[j]];
        const double value = _Tet.VolumeStiffness[t](i, j);
        if (column >= 0) triplets.emplace_back(row, column, value);
        else system.BoundaryRhs[row] -= value * fixedValue[tet[j]];
      }
    }
  }
  if (freeCount == 0) return true;
  Eigen::SparseMatrix<double> stiffness(freeCount, freeCount);
  stiffness.setFromTriplets(triplets.begin(), triplets.end());
  system.Solver.compute(stiffness);
  if (system.Solver.info() != Eigen::Success) {
    error = "factorization failed";
    return false;
  }
  return true;
}

ZSPACE_INLINE bool HarmonicSlicer::solveTetDirichlet(TetDirichletSystem &system,
    zFloatArray &outField, std::string &error,
    const std::vector<Eigen::Vector3d> *targetGradients) const {
  outField = system.FixedField;
  if (system.BoundaryRhs.size() == 0) return true;
  Eigen::VectorXd rhs = system.BoundaryRhs;
  if (targetGradients) {
    for (size_t t = 0; t < _Tet.Volumes.size(); ++t) {
      const zTetrahedron::TetVolume &element = _Tet.Volumes[t];
      if (!element.Valid) continue;
      const Eigen::Vector4d source = element.Volume * (_Tet.VolumeBasisGradients[t].transpose() * (*targetGradients)[t]);
      for (int i = 0; i < 4; ++i) {
        const int row = system.FreeIndex[_Tet.Volumes[t].VertexIds[i]];
        if (row >= 0) rhs[row] += source[i];
      }
    }
  }
  const Eigen::VectorXd solution = system.Solver.solve(rhs);
  if (system.Solver.info() != Eigen::Success || !solution.allFinite()) {
    error = "returned an invalid field";
    return false;
  }
  for (size_t i = 0; i < system.FreeIndex.size(); ++i) {
    const int row = system.FreeIndex[i];
    if (row >= 0) outField[i] = std::clamp(static_cast<float>(solution[row]), 0.0f, 1.0f);
  }
  return true;
}

ZSPACE_INLINE Eigen::Vector3d HarmonicSlicer::tetGradient(size_t tet, const zFloatArray &field) const {
  Eigen::Vector3d gradient = Eigen::Vector3d::Zero();
  const zTetrahedron::TetVolume &element = _Tet.Volumes[tet];
  if (!element.Valid) return gradient;
  for (int local = 0; local < 4; ++local)
    gradient += static_cast<double>(field[_Tet.Volumes[tet].VertexIds[local]]) * _Tet.VolumeBasisGradients[tet].col(local);
  return gradient;
}

// Gradient equalisation (footing comparison sketch, Method 3 stage A): the
// target magnitude g is the volume-weighted mean |grad| of the input field
// (the constant its [0, 1] Dirichlet range can best accommodate); each
// iteration re-solves min sum_t vol_t |grad psi - g d_t|^2 with the same caps,
// d_t the unit gradient direction of the previous iterate (the harmonic
// field's on the first). Same field range and caps as the harmonic field, so
// the layering (layer count, iso values) is unchanged - only where the
// evenly spaced iso surfaces sit.
ZSPACE_INLINE bool HarmonicSlicer::equaliseField(TetDirichletSystem &system, zFloatArray &ioField, std::string &outReport) const {
  const size_t tetCount = _Tet.Volumes.size();
  const auto gradientSpread = [&](const zFloatArray &field, double &outMean, double &outSpread) {
    double minG = (std::numeric_limits<double>::max)(), maxG = 0.0, sum = 0.0, volume = 0.0;
    for (size_t t = 0; t < tetCount; ++t) {
      if (!_Tet.Volumes[t].Valid) continue;
      const double g = tetGradient(t, field).norm();
      minG = (std::min)(minG, g);
      maxG = (std::max)(maxG, g);
      sum += _Tet.Volumes[t].Volume * g;
      volume += _Tet.Volumes[t].Volume;
    }
    outMean = volume > 0.0 ? sum / volume : 0.0;
    outSpread = outMean > 0.0 ? (maxG - minG) / outMean : 0.0;
  };

  double target = 0.0, spreadBefore = 0.0;
  gradientSpread(ioField, target, spreadBefore);
  if (target <= 0.0) {
    outReport = "no usable gradient";
    return false;
  }

  std::vector<Eigen::Vector3d> targets(tetCount, Eigen::Vector3d::Zero());
  for (int iteration = 0; iteration < i_Settings._HEqualiseIterations; ++iteration) {
    for (size_t t = 0; t < tetCount; ++t) {
      const Eigen::Vector3d gradient = tetGradient(t, ioField);
      const double magnitude = gradient.norm();
      targets[t] = magnitude > 1e-15 ? Eigen::Vector3d(target * gradient / magnitude) : Eigen::Vector3d::Zero();
    }
    zFloatArray psi;
    std::string error;
    if (!solveTetDirichlet(system, psi, error, &targets)) {
      outReport = "iteration " + std::to_string(iteration) + " " + error;
      return false;
    }
    ioField.swap(psi);
  }

  double meanAfter = 0.0, spreadAfter = 0.0;
  gradientSpread(ioField, meanAfter, spreadAfter);
  std::ostringstream out;
  out << "equalised " << i_Settings._HEqualiseIterations << "x, target |grad|=" << target
      << ", spread (max-min)/mean " << spreadBefore << " -> " << spreadAfter;
  outReport = out.str();
  return true;
}

ZSPACE_INLINE std::vector<double> HarmonicSlicer::capDirichletValues() const {
  std::vector<double> fixedValue(_Tet.Vertices.size(), -1.0);
  for (size_t id = 0; id < _Tet.vertextag.size(); ++id) {
    if (_Tet.vertextag[id] == static_cast<int>(VertexTag::Bottom)) fixedValue[id] = 0.0;
    else if (_Tet.vertextag[id] == static_cast<int>(VertexTag::Top)) fixedValue[id] = 1.0;
  }
  return fixedValue;
}

// The harmonic field w (Bottom = 0, Top = 1, harmonic inside), optionally
// gradient-equalised (_HEqualiseIterations). w is written straight into
// _Tet.UVW's z instead of its own array: _Tet.UVW is allocated here with
// (u,v) left as NaN placeholders until computeVolumetricParameterization()
// fills them in.
ZSPACE_INLINE bool HarmonicSlicer::computeHarmonicField() {
  zTsUtils::ScopedSlicingTiming methodTiming("HarmonicSlicer.computeHarmonicField");
  const bool success = computeHarmonicFieldImpl();
  if (success) zTsUtils::printSuccess("[HarmonicSlicer::computeHarmonicField] Success");
  else zTsUtils::printError("[HarmonicSlicer::computeHarmonicField] Failed");
  return success;
}

ZSPACE_INLINE bool HarmonicSlicer::computeHarmonicFieldImpl() {
  clearParameterization();
  const float invalid = std::numeric_limits<float>::quiet_NaN();
  _Tet.UVW.assign(_Tet.Vertices.size(), zPoint(invalid, invalid, invalid));
  _Tet.UVWValid.assign(_Tet.Vertices.size(), false);
  _Tet.FieldValid = false;
  updateFieldGradients();
  if (_Tet.Volumes.empty()) {
    _Status = "computeHarmonicField failed: no tet volume. Call computeBasicStructure() first.";
    printDebug("[computeHarmonicFieldImpl] " + _Status);
    return false;
  }

  zFloatArray phi;
  TetDirichletSystem system;
  std::string error;
  if (!prepareTetDirichlet(capDirichletValues(), system, error) ||
      !solveTetDirichlet(system, phi, error)) {
    _Status = "computeHarmonicField failed: harmonic solve " + error + ".";
    printDebug("[computeHarmonicFieldImpl] " + _Status);
    return false;
  }

  // Optional gradient equalisation: the slicing field (w) becomes psi; the
  // tracing, contours and layering downstream read it unchanged.
  std::string equaliseReport;
  if (i_Settings._HEqualiseIterations > 0 && !equaliseField(system, phi, equaliseReport)) {
    _Status = "computeHarmonicField failed: gradient equalisation " + equaliseReport + ".";
    printDebug("[computeHarmonicFieldImpl] " + _Status);
    return false;
  }
  if (phi.empty()) return false;

  const float invalidUV = std::numeric_limits<float>::quiet_NaN();
  _Tet.UVW.assign(phi.size(), zPoint(invalidUV, invalidUV, 0.0));
  for (size_t i = 0; i < phi.size(); ++i) _Tet.UVW[i].z = phi[i];
  _Tet.UVWValid.assign(phi.size(), false);
  _Tet.FieldValid = true;
  updateFieldGradients();

  std::ostringstream out;
  out << "computeHarmonicField: tet vertices=" << _Tet.Vertices.size() << " tets=" << _Tet.Volumes.size()
      << (equaliseReport.empty() ? std::string(" | plain harmonic field") : " | " + equaliseReport);
  _Status = out.str();

  return true;
}

} // namespace zSpace
