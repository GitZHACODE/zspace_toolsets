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

// Connectivity stays fixed during relaxation. Each boundary edge retains its
// input width ID and is constrained to a connected planar input-side patch.
// Vertices at patch junctions slide along their shared reference edges.
struct Layer {
  zPointArray Positions, UV;
  std::vector<std::array<int, 3>> Triangles;
  std::unordered_map<unsigned long long, int> BoundaryWidths;
  zFloatArray Heights;
  zVectorArray Normals; // empty on fixed caps; preserve their face normals
};

struct Settings {
  int Iterations = 8;
  double Strength = 0.25;
  // 0 selects the separation-based limit. An explicit value is additionally
  // bounded by 20% of the initial distance to either neighbouring layer.
  double MaxDeviation = 0.0;
  bool Parallel = true;
};

struct Report {
  int AcceptedIterations = 0;
  int MovedBoundaryVertices = 0;
  int PlaneHeightEstimates = 0;
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

// A null UV pointer asks only whether the point is in the input volume. A
// non-null pointer must receive finite transported (u,v) at that position.
using SampleVolume = std::function<bool(const zPoint &, zPoint *)>;

ZSPACE_TOOLSETS bool relax(const Reference &reference, std::vector<Layer> &layers,
                          const Settings &settings, const SampleVolume &sample,
                          Report &report, std::string &error);

} } // namespace zSpace::zTsRelax
