#include <zspace/slicing.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

using namespace zSpace;
namespace {
void require(bool value, const std::string &message) {
  if (!value) throw std::runtime_error(message);
}
zPointArray positions(const PrintPath &path) {
  zPointArray p;
  zFnGraph(const_cast<zObjGraph &>(path.graph)).getVertexPositions(p);
  return p;
}
PrintPath path(zPointArray p, bool closed = false) {
  PrintPath out;
  zIntArray edges;
  for (int i = 1; i < int(p.size()); ++i) { edges.push_back(i-1); edges.push_back(i); }
  if (closed) { edges.push_back(int(p.size())-1); edges.push_back(0); }
  zFnGraph(out.graph).create(p, edges);
  for (size_t i = 0; i < p.size(); ++i) {
    out.widths.push_back(0.2f + 0.1f * float(i));
    out.heights.push_back(1.0f + float(i));
    out.normals.push_back(zVector(0,0,1));
    out.trackIds.push_back(i == 0 ? 0.0f : 1.0f);
  }
  return out;
}
size_t findPoint(const PrintPath &p, const zPoint &target, double tolerance = 1e-7) {
  const auto points = positions(p);
  for (size_t i = 0; i < points.size(); ++i)
    if (zTsUtils::vecDistSq(points[i], target) <= tolerance*tolerance) return i;
  throw std::runtime_error("Missing required point (" + std::to_string(target.x) + ", " +
      std::to_string(target.y) + ", " + std::to_string(target.z) + ")");
}
void aligned(const PrintPath &p) {
  const size_t n = positions(p).size();
  require(p.widths.size() == n && p.heights.size() == n && p.normals.size() == n &&
      p.trackIds.size() == n && p.featureFlags.size() == n, "Per-point attributes lost alignment");
  zPointArray polyline;
  zScalarArray arc;
  if (!zTsUtils::graphtopolyline(p.graph, polyline, arc)) {
    zIntArray edges;
    zFnGraph(const_cast<zObjGraph &>(p.graph)).getEdgeData(edges);
    const auto points = positions(p);
    std::string detail = "Invalid path connectivity: " + std::to_string(n) + " vertices";
    for (size_t e = 0; e < edges.size(); e += 2)
      if (zTsUtils::vecDistSq(points[edges[e]], points[edges[e+1]]) == 0)
        detail += "; zero-length edge " + std::to_string(edges[e]) + "-" + std::to_string(edges[e+1]);
    throw std::runtime_error(detail);
  }
}
PrintPath resample(PrintPath input, float spacing, float angle) {
  PrintPathArray layers{{input}};
  PrintSynSetting settings;
  settings.ResampleLength = spacing; settings.FeatureAngleTol = angle;
  std::string status;
  std::vector<int> bad;
  require(zPrintSyn::resamplePrintPaths(layers, settings, status, &bad, false), status);
  require(bad.empty(), "Resampling opened a closed curve");
  aligned(layers[0][0]);
  return layers[0][0];
}
void resampling() {
  // Strict threshold semantics: equality is intentionally excluded.
  const auto corner = path({{0,0,0},{1.005,0,0},{1.005,1,0}});
  require(zPrintSyn::cornerVertexIds(corner, 89).empty(), "Below threshold selected a corner");
  require(zPrintSyn::cornerVertexIds(corner, 90).empty(), "Threshold equality changed");
  require(zPrintSyn::cornerVertexIds(corner, 91) == zIntArray{1}, "Above threshold missed a corner");
  auto sampled = resample(corner, 0.2f, 91);
  size_t id = findPoint(sampled, {1.005,0,0});
  require(sampled.featureFlags[id] & PrintPath::AngleFeature, "Lost angle identity");
  require(std::abs(sampled.widths[id] - corner.widths[1]) < 1e-7 &&
      sampled.heights[id] == corner.heights[1] && sampled.trackIds[id] == corner.trackIds[1],
      "Nearby regular sample replaced a corner's attributes");

  // Features separated by <0.01 mm and <0.1 mm across closure are distinct.
  for (bool clockwise : {false, true}) {
    zPointArray p{{0,0,0},{1,0,0},{1,0.005,0},{0,0.005,0}};
    if (clockwise) std::reverse(p.begin(), p.end());
    auto result = resample(path(p, true), 0.2f, 91);
    for (const auto &point : p)
      require(result.featureFlags[findPoint(result, point)] & PrintPath::AngleFeature,
          "Thin-loop or closure feature was removed/reordered incorrectly");
    require(zTsUtils::isClosedGraph(result.graph), "Closed loop was opened");
  }
  // Short edges cannot disable angle detection; fine sampling must remain fine.
  zPointArray tiny{{0,0,0},{0.0004,0,0},{0.0004,0.0003,0},{0,0.0003,0}};
  sampled = resample(path(tiny, true), 0.0001f, 91);
  for (const auto &point : tiny)
    require(sampled.featureFlags[findPoint(sampled, point, 1e-10)] & PrintPath::AngleFeature,
        "Short feature-to-feature edge lost a corner");
  require(positions(sampled).size() > 4, "Fine regular samples were filtered away");

  // Consecutive coincident feature vertices coalesce, combining their masks.
  auto coincident = path({{0,0,0},{1,0,0},{1,0,0},{2,0,0}});
  coincident.featureFlags = {0, PrintPath::AngleFeature, PrintPath::IntersectionFeature, 0};
  sampled = resample(coincident, 0.3f, 0);
  require(sampled.featureFlags[findPoint(sampled, {1,0,0})] ==
      (PrintPath::AngleFeature | PrintPath::IntersectionFeature), "Coincident feature flags lost");

  // Reconstruction intersects the horizontal/vertical supports at (2,0).
  auto bevel = path({{0,0,0},{1,0,0},{1.8,0,0},{2,0.2,0},{2,1,0},{2,2,0}});
  require(zPrintSyn::reconstructCorners(bevel, 150), "Corner reconstruction failed");
  id = findPoint(bevel, {2,0,0});
  require(bevel.featureFlags[id] & PrintPath::IntersectionFeature, "Intersection not marked");
  const float cornerWidth = bevel.widths[id], cornerHeight = bevel.heights[id];
  // A lower threshold deliberately does not re-detect this 90-degree corner.
  sampled = resample(bevel, 0.3f, 80);
  id = findPoint(sampled, {2,0,0});
  require(sampled.featureFlags[id] & PrintPath::IntersectionFeature, "Intersection required re-detection");
  require(sampled.widths[id] == cornerWidth && sampled.heights[id] == cornerHeight,
      "Reconstructed intersection attributes changed");
  sampled = resample(sampled, 0.23f, 150);
  require(sampled.featureFlags[findPoint(sampled, {2,0,0})] & PrintPath::IntersectionFeature,
      "Repeated reconstruction moved an existing intersection anchor");
  // Endpoint survival even when spacing exceeds the whole open path length.
  sampled = resample(path({{0,0,0},{0.005,0,0}}), 1, 90);
  require(positions(sampled).size() == 2, "Short open path lost its endpoint");

  PrintPathArray invalid{{corner}};
  invalid[0][0].featureFlags = {1};
  std::string status;
  require(!zPrintSyn::resamplePrintPaths(invalid, PrintSynSetting{}, status, nullptr, false),
      "Malformed feature metadata accepted");
  PrintSynSetting badSettings;
  badSettings.FeatureAngleTol = std::numeric_limits<float>::quiet_NaN();
  invalid = {{corner}};
  require(!zPrintSyn::resamplePrintPaths(invalid, badSettings, status, nullptr, false),
      "Nonfinite feature threshold accepted");
}

zObjMesh chart() {
  zObjMesh mesh;
  zPointArray p{{0,0,0},{2,0,0},{2,2,0},{0,2,0}};
  zIntArray counts{3,3}, faces{0,1,2,0,2,3};
  zFnMesh(mesh).create(p, counts, faces);
  return mesh;
}
// Independent analytic piecewise-affine surface: stretched XY, folded across
// the diagonal, rigidly rotated about Y, then translated per layer.
zPoint expected(const zPoint &p, double layer) {
  const double x = 2*p.x, y = 0.5*p.y, z = 0.4*(std::max)(0.0, double(p.y)-p.x);
  return zPoint(10 + std::cos(0.5)*x + std::sin(0.5)*z, -3+y,
      layer - std::sin(0.5)*x + std::cos(0.5)*z);
}
zObjMesh section(double layer) {
  auto mesh = chart();
  zFnMesh fn(mesh);
  auto *p = fn.getRawVertexPositions();
  for (int i = 0; i < fn.numVertices(); ++i) p[i] = expected(p[i], layer);
  return mesh;
}
void mapping() {
  auto input = resample(path({{0.2,0.2,0},{1.005,0.2,0},{1.005,1.6,0},{0.2,1.6,0}}, true), 0.2f, 91);
  input.trackIds.assign(input.trackIds.size(), 0.0f);
  const auto source = positions(input);
  zObjMeshArray charts{chart(), chart()}, sections{section(0), section(1)};
  PrintPathArray mapped;
  std::string status;
  std::vector<zFloatArray> heights{{0.2f,0.4f,0.6f,0.8f},{1,1,1,1}};
  require(zPrintSyn::projectPrintPaths({{input},{input}}, charts, sections, nullptr,
      mapped, status, nullptr, false, &heights), status);
  double maxError = 0;
  for (int l = 0; l < 2; ++l) {
    aligned(mapped[l][0]);
    const auto points = positions(mapped[l][0]);
    require(points.size() == source.size() && mapped[l][0].featureFlags == input.featureFlags &&
        mapped[l][0].widths == input.widths && mapped[l][0].trackIds == input.trackIds,
        "Projection changed target identity/attributes");
    for (size_t i = 0; i < source.size(); ++i) {
      const double error = std::sqrt(zTsUtils::vecDistSq(points[i], expected(source[i], l)));
      maxError = (std::max)(maxError, error);
      require(error < 3e-6, "Incorrect analytic 2D-to-3D point mapping");
      const auto &p = source[i];
      const double h = l ? 1.0 : (p.x >= p.y ? 0.2+0.1*p.x+0.1*p.y : 0.2-0.1*p.x+0.3*p.y);
      require(std::abs(mapped[l][0].heights[i]-h) < 1e-6, "Height interpolation mismatch");
      require(std::abs(mapped[l][0].normals[i].length()-1) < 1e-6, "Nonunit mapped normal");
    }
  }
  std::cout << "Analytic folded/rotated mapping max error: " << maxError << " model units\n";
  PrintPath combined;
  require(zPrintSyn::linkPrintPaths(mapped, combined, status), status);
  require(combined.featureFlags.size() == source.size()*2, "Linking lost feature masks");
  for (int l = 0; l < 2; ++l) for (size_t i = 0; i < source.size(); ++i)
    require(combined.featureFlags[l*source.size()+i] == input.featureFlags[i], "Linked mask order changed");
  require(zPrintSyn::linkPrintPaths({{mapped[0][0], mapped[1][0]}}, combined, status), status);
  std::vector<zTsUtils::GraphPolyline> components;
  require(zTsUtils::graphPolylines(combined.graph, components) && components.size() == 2 &&
      combined.featureFlags.size() == source.size()*2, "Feature propagation linked disjoint curves");

  // Both seam mapping strategies must preserve features outside explicit gaps.
  for (auto strategy : {zPrintSyn::SeamMapping::SharedTopology, zPrintSyn::SeamMapping::WorldNearest}) {
    PrintPathArray paths{{input},{input}};
    PrintSynSetting settings;
    settings.ResampleLength = 0.17f; settings.FeatureAngleTol = 0;
    settings.SeamOffset = 0.11f;
    require(zPrintSyn::computePipeline(charts, sections, nullptr, strategy, settings,
        paths, mapped, combined, status, false, &heights), status);
    for (int l = 0; l < 2; ++l) {
      aligned(paths[l][0]); aligned(mapped[l][0]);
      for (size_t i = 0; i < source.size(); ++i) if (input.featureFlags[i]) {
        const auto id = findPoint(mapped[l][0], expected(source[i], l), 3e-6);
        require(mapped[l][0].featureFlags[id] == input.featureFlags[i], "Seam/projection lost feature");
        const auto linked = findPoint(combined, expected(source[i], l), 3e-6);
        require(combined.featureFlags[linked] == input.featureFlags[i], "Final linked target lost feature");
      }
    }
  }
  // An explicit seam may remove only the cut arc; surviving corner metadata
  // and the new endpoint must stay aligned. The last square corner is cut.
  auto square = path({{0.2,0.2,0},{1.8,0.2,0},{1.8,1.8,0},{0.2,1.8,0}}, true);
  square.featureFlags.assign(4, PrintPath::AngleFeature);
  square.trackIds.assign(4, 0);
  PrintPathArray gap{{square}};
  std::vector<int> bad;
  require(zPrintSyn::alignStartPoints(gap, zObjMeshArray{chart()}, 0, 1.8,
      status, &bad, false) && bad.empty(), status);
  aligned(gap[0][0]);
  require(!zTsUtils::isClosedGraph(gap[0][0].graph) &&
      std::count_if(gap[0][0].featureFlags.begin(), gap[0][0].featureFlags.end(),
          [](int flag) { return flag != 0; }) == 3, "Explicit seam cut wrong feature set");
  for (const zPoint &p : zPointArray{{0.2,0.2,0},{1.8,0.2,0},{1.8,1.8,0}})
    require(gap[0][0].featureFlags[findPoint(gap[0][0], p)] & PrintPath::AngleFeature,
        "Seam removed a feature outside the cut arc");
  require(gap[0][0].featureFlags[findPoint(gap[0][0], {0.4,1.8,0}, 1e-6)] == 0,
      "Interpolated seam endpoint acquired a false feature mask");
  gap = {{square}}; gap[0][0].featureFlags.pop_back();
  require(!zPrintSyn::alignStartPoints(gap, zObjMeshArray{chart()}, 0, 0,
      status, nullptr, false), "Seam accepted malformed feature metadata");

  auto invalid = input;
  auto p = positions(invalid); p[0].x = 5;
  zIntArray edges;
  zFnGraph(invalid.graph).getEdgeData(edges);
  zFnGraph(invalid.graph).create(p, edges);
  require(!zPrintSyn::projectPrintPaths({{invalid}}, {chart()}, {section(0)}, nullptr,
      mapped, status, &bad, false), "Out-of-chart feature was silently snapped");
  require(bad == std::vector<int>{0} && status.find("layer 0, path 0: vertex 0") != std::string::npos,
      "Mapping diagnostic omitted layer/path/vertex: " + status);
  require(positions(mapped[0][0]).empty(), "Mapping failure retained an invalid projected path");
  require(!zPrintSyn::projectPrintPaths({{invalid}}, {chart()}, {section(0)}, nullptr,
      mapped, status, nullptr, false), "Omitting bad-layer output hid mapping failure");
  auto mismatched = section(0);
  zPointArray sectionPoints;
  zFnMesh(mismatched).getVertexPositions(sectionPoints);
  zIntArray quadCounts{4}, quadFace{0,1,2,3};
  zFnMesh(mismatched).create(sectionPoints, quadCounts, quadFace);
  require(!zPrintSyn::projectPrintPaths({{input}}, {chart()}, {mismatched}, nullptr,
      mapped, status, nullptr, false) && status.find("TopologyMismatch") != std::string::npos,
      "Mismatched 2D/3D face connectivity accepted");

  auto overlapping = chart();
  zPointArray duplicate{{0,0,0},{2,0,0},{0,2,0},{0,0,0},{2,0,0},{0,2,0}};
  zIntArray counts{3,3}, faces{0,1,2,3,4,5};
  zFnMesh(overlapping).create(duplicate, counts, faces);
  auto probe = path({{0.2,0.2,0},{0.4,0.2,0}});
  std::string error;
  require(!zTsUtils::remapGraph(probe.graph, overlapping, overlapping, nullptr, nullptr, nullptr, &error) &&
      error.find("AmbiguousChartOverlap") != std::string::npos, "Overlapping chart accepted: " + error);
}
void slicerIntegration(zTsSlicerMode mode) {
  zObjMesh mesh;
  zPointArray p{{0,0,0},{4,0,0},{4,4,0},{0,4,0},
      {0,0,1},{4,0,1},{4,4,1},{0,4,1},{0,0,2},{4,0,2},{4,4,2},{0,4,2}};
  zIntArray counts(10,4), faces{0,3,2,1,8,9,10,11,0,1,5,4,1,2,6,5,2,3,7,6,3,0,4,7,
      4,5,9,8,5,6,10,9,6,7,11,10,7,4,8,11}, tags(10,Side);
  tags[0] = Bottom; tags[1] = Top;
  zFnMesh(mesh).create(p, counts, faces);
  zIntArray serialFlags;
  for (bool parallel : {false,true}) {
    zTsSlicing session;
    session.setInputMesh(mesh, tags); session.setSlotAnchor(zPoint(0,0,0)); session.setSlicerMode(mode);
    session.getSettings().LayerSpacing = 0.5f;
    session.getSettings().BoundaryWidth = {0.4f}; session.getSettings().FieldCellSize = 0.1f;
    session.getPrintSynSetting().FeatureAngleTol = 140;
    session.getPrintSynSetting().ResampleLength = 0.2f;
    require(session.computeSlicer(parallel), session.getStatus());
    require(session.computePostProcess(parallel), session.getStatus());
    const auto &two = session.get2DPrintPaths(), &three = session.get3DPrintPaths();
    const auto &combined = session.getResultantPrintPath();
    aligned(combined);
    int features = 0;
    for (size_t l = 0; l < two.size(); ++l) for (size_t r = 0; r < two[l].size(); ++r) {
      aligned(three[l][r]);
      const auto points = positions(three[l][r]);
      require(two[l][r].featureFlags == three[l][r].featureFlags, "Slicer projection changed feature identity");
      for (size_t i = 0; i < points.size(); ++i) if (three[l][r].featureFlags[i]) {
        ++features;
        const size_t id = findPoint(combined, points[i], 1e-6);
        require(combined.featureFlags[id] == three[l][r].featureFlags[i], "Slicer linking lost feature identity");
      }
    }
    require(features > 0, "Slicer integration fixture produced no protected features");
    if (!parallel) serialFlags = combined.featureFlags;
    else require(serialFlags == combined.featureFlags, "Serial/parallel feature masks differ");
  }
}
}
int main() {
  try {
    resampling(); mapping(); slicerIntegration(zTsSlicerMode::Topology);
#ifdef ZSPACE_TOOLSETS_WITH_TETGEN
    slicerIntegration(zTsSlicerMode::Harmonic);
#endif
    return 0;
  }
  catch (const std::exception &e) { std::cerr << e.what() << '\n'; return 1; }
}
