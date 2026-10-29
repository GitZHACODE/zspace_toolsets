#pragma once

#include <zspace/zToolsets/fabrication/slicing/utils/Utils.h>
#include <functional>
#include <unordered_map>

namespace zSpace { namespace zTsRelax {

// All geometry is in model units. The reference is immutable and includes the
// caps; WidthIds is -1 on caps and the input side-width group on side triangles.
struct Reference {
  zPointArray Positions;
  std::vector<std::array<int, 3>> Triangles;
  zIntArray WidthIds;
};

// Connectivity stays fixed after optional quality preconditioning. Each boundary edge retains its
// input width ID and is constrained to a connected planar input-side patch.
// Vertices at patch junctions slide along their shared reference edges.
struct Layer {
  zPointArray Positions, UV;
  std::vector<std::array<int, 3>> Triangles;
  std::unordered_map<unsigned long long, int> BoundaryWidths;
  // TCP clearance along -Normals to the preceding section (S1..Slast).
  // Unprinted S0 stores an optional next-gap diagnostic, zero if unavailable.
  zFloatArray Heights;
  zVectorArray Normals; // recovered TCP normals, including fixed top; S0 is unprinted
};

struct Settings {
  int Iterations = 8;
  double Strength = 0.25;
  // 0 selects the separation-based limit. An explicit value is additionally
  // bounded by 20% of the initial distance to either neighbouring layer.
  double MaxDeviation = 0.0;
  // Optional targeted quality repair before fairing/unrolling. Zero preserves
  // connectivity. Caps and boundary vertices are never collapsed.
  double MinEdgeLength = 0.0;
  bool Parallel = true;
};

struct Report {
  int AcceptedIterations = 0;
  int MovedBoundaryVertices = 0;
  int PlaneHeightEstimates = 0;
  int UnavailableAuxiliaryGaps = 0;
  int RecoveredBoundaryTags = 0;
  int UVConstrainedVertices = 0;
  int SliverConstrainedVertices = 0;
  int GeometryConstrainedVertices = 0;
  std::string LastRejectedReason;
  double MaxDisplacement = 0.0;
  double MaxBoundaryError = 0.0;
  double MinimumInitialGap = 0.0;
  // Boundary/setup, gap bounds, sparse targets, trials, normals/heights, total.
  std::array<double,6> StageMilliseconds{};
};

struct HeightReport {
  int AcceptedIterations = 0, ViolatingVertices = 0;
  int Factorizations = 0, OuterIterations = 0;
  int DegreesOfFreedom = 0, InnerIterations = 0, LineSearchTrials = 0;
  int InitialSliverVertices = 0, UncontrollableSamples = 0, FixedVertices = 0;
  int UnavailableAuxiliaryGaps = 0;
  int InitialSections = 0, FinalSections = 0;
  bool Converged = false;
  double Minimum = 0, Maximum = 0, WorstViolation = 0;
  double InitialMinimum = 0, InitialMaximum = 0, Milliseconds = 0;
  double MaxDisplacement = 0, MaxBoundaryError = 0;
  double InitialBendingEnergy = 0, FinalBendingEnergy = 0;
  double MinimumDirectedGap = 0, MaximumDirectedGap = 0;
  std::string LastRejectedReason;
};

// A null UV pointer asks only whether the point is in the input volume. A
// non-null pointer must receive finite transported (u,v) at that position.
using SampleVolume = std::function<bool(const zPoint &, zPoint *)>;

ZSPACE_TOOLSETS bool relax(const Reference &reference, std::vector<Layer> &layers,
                          const Settings &settings, const SampleVolume &sample,
                          Report &report, std::string &error);

// Ordered low-dimensional initializer. Endpoints remain 0/1; the geometric
// height evaluator, not gradient values or clamped output attributes, drives it.
ZSPACE_TOOLSETS bool redistributeIsoValues(std::vector<Layer> layers,
    const PrintHeightOptimizationSettings &settings, std::vector<double> &levels,
    bool parallel, std::string &error);

// Coupled smooth correction of relaxed meshes. A valid geometry may retain
// violations when constrained: inspect Converged and validate final paths.
ZSPACE_TOOLSETS bool optimizeHeightDomain(const Reference &reference,
    std::vector<Layer> &layers, const PrintHeightOptimizationSettings &settings,
    const SampleVolume &sample, bool parallel, HeightReport &report, std::string &error,
    const std::vector<zPointArray> *printTargets = nullptr,
    double fairingStrength = 0.25);

// Exact directed TCP rays. targets/normals p belong to section p+1 and cast
// along -normal to section p; finite nonzero normals are normalized internally.
// A miss is unsupported, never plane-extrapolated.
ZSPACE_TOOLSETS bool measureTargetHeights(const std::vector<Layer> &layers,
    const std::vector<zPointArray> &targets,
    const std::vector<zVectorArray> &targetNormals,
    std::vector<zFloatArray> &outHeights,std::string &error,bool parallel=true);

} } // namespace zSpace::zTsRelax
