#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <zspace/slicing.h>
using namespace zSpace;
static void require(bool ok, const char *message) {
  if (!ok)
    throw std::runtime_error(message);
}
static PrintPath path(zPointArray points, zIntArray edges) {
  PrintPath p;
  zFnGraph(p.graph).create(points, edges);
  p.widths.assign(points.size(), 2);
  p.heights.assign(points.size(), 4);
  p.normals.assign(points.size(),
                   zVector(0, 0, 2)); // normalization must preserve depth
  return p;
}
static void checkClosed(zObjMesh &mesh, size_t vertices, size_t faces) {
  zPointArray points;
  zIntArray ids, counts;
  zFnMesh fn(mesh);
  fn.getVertexPositions(points);
  fn.getPolygonData(ids, counts);
  require(points.size() == vertices && counts.size() == faces,
          "Unexpected loft topology");
  for (const auto &p : points)
    require(std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z),
            "Invalid mesh coordinate");
  std::map<std::pair<int, int>, int> directed;
  size_t offset = 0;
  for (int count : counts) {
    require(count == 4, "Expected quads");
    for (int i = 0; i < count; ++i) {
      const int a = ids[offset + i], b = ids[offset + (i + 1) % count];
      require(a >= 0 && b >= 0 && size_t(a) < vertices && size_t(b) < vertices,
              "Invalid mesh index");
      ++directed[{a, b}];
    }
    offset += count;
  }
  for (const auto &e : directed)
    require(e.second == 1 && directed[{e.first.second, e.first.first}] == 1,
            "Mesh edge winding/closure mismatch");
}
int main() {
  try {
    auto straight = path({{0, 0, 0}, {10, 0, 0}}, {0, 1});
    straight.widths[1] = 6;
    straight.heights[1] = 2;
    zObjMesh mesh;
    std::string error;
    require(zTsUtils::buildPrintMesh({{straight}}, mesh, &error),
            error.c_str());
    checkClosed(mesh, 8, 6);
    zPointArray points;
    zFnMesh(mesh).getVertexPositions(points);
    const zPointArray expected = {{0, -1, 0}, {0, 1, 0},   {0, 1, 4},
                                  {0, -1, 4}, {10, -3, 0}, {10, 3, 0},
                                  {10, 3, 2}, {10, -3, 2}};
    for (size_t i = 0; i < points.size(); ++i)
      require(zTsUtils::vecDistSq(points[i], expected[i]) < 1e-12,
              "Rectangle width centering/one-sided height mismatch");
    // Width and height axes rotate with the supplied normal, rather than world
    // Z.
    auto tilted = straight;
    tilted.normals.assign(2, zVector(0, 1, 0));
    require(zTsUtils::buildPrintMesh({{tilted}}, mesh, &error), error.c_str());
    zFnMesh(mesh).getVertexPositions(points);
    require(zTsUtils::vecDistSq(points[0], zPoint(0, 0, 1)) < 1e-12 &&
                zTsUtils::vecDistSq(points[2], zPoint(0, 4, -1)) < 1e-12,
            "Normal was ignored");
    // A sloping path requires projecting the normal into the section plane.
    auto slope = path({{0, 0, 0}, {10, 0, 10}}, {0, 1});
    require(zTsUtils::buildPrintMesh({{slope}}, mesh, &error), error.c_str());
    zFnMesh(mesh).getVertexPositions(points);
    const double rise = 2 * std::sqrt(2.0);
    require(zTsUtils::vecDistSq(points[0], zPoint(0, -1, 0)) < 1e-10 &&
                zTsUtils::vecDistSq(points[2], zPoint(-rise, 1, rise)) < 1e-10,
            "Section is not perpendicular to the path with projected height");
    checkClosed(mesh, 8, 6);
    // A 3:4 corner follows the neighbor chord, not the unit-edge bisector.
    auto corner = path({{0, 0, 0}, {3, 0, 0}, {3, 4, 0}},
                       {0, 1, 1, 2});
    require(zTsUtils::buildPrintMesh({{corner}}, mesh, &error), error.c_str());
    zFnMesh(mesh).getVertexPositions(points);
    require(zTsUtils::vecDistSq(points[4], zPoint(3.8, -0.6, 0)) < 1e-10 &&
                zTsUtils::vecDistSq(points[5], zPoint(2.2, 0.6, 0)) < 1e-10 &&
                zTsUtils::vecDistSq(points[6], zPoint(2.2, 0.6, 4)) < 1e-10,
            "Corner chord weighting or forward section ordering mismatch");
    checkClosed(mesh, 12, 10);
    auto loop = path({{0, 0, 0}, {10, 0, 0}, {10, 10, 0}, {0, 10, 0}},
                     {0, 1, 1, 2, 2, 3, 3, 0});
    require(
        zTsUtils::buildPrintMesh({{loop}, {straight, straight}}, mesh, &error),
        error.c_str());
    checkClosed(mesh, 32,
                28); // loop has 16 side quads; each independent open path has 6
    zIntArray ids, counts;
    zFnMesh(mesh).getPolygonData(ids, counts);
    for (size_t f = 0; f < counts.size(); ++f) {
      const int component = ids[4 * f] / 16;
      for (int j = 1; j < 4; ++j)
        require(ids[4 * f + j] / 16 == component,
                "Loft bridged separate layers");
      if (f >= 16)
        for (int j = 1; j < 4; ++j)
          require((ids[4 * f + j] - 16) / 8 == (ids[4 * f] - 16) / 8,
                  "Loft bridged separate paths");
    }
    straight.heights.assign(2, 0);
    require(zTsUtils::buildPrintMesh({{straight}}, mesh, &error),
            error.c_str());
    zFnMesh(mesh).getVertexPositions(points);
    for (const auto &p : points)
      require(p.z == 0, "Zero height gained artificial thickness");
    straight.widths.assign(2, 0);
    require(zTsUtils::buildPrintMesh({{straight}}, mesh, &error),
            error.c_str());
    zFnMesh(mesh).getVertexPositions(points);
    for (const auto &p : points)
      require(p.y == 0, "Zero width gained artificial thickness");
    auto invalid = straight;
    invalid.widths.pop_back();
    require(!zTsUtils::buildPrintMesh({{invalid}}, mesh, &error) &&
                !error.empty() && zFnMesh(mesh).numVertices() == 0,
            "Missing attributes retained stale mesh");
    invalid = straight;
    invalid.heights[0] = -1;
    require(!zTsUtils::buildPrintMesh({{invalid}}, mesh, &error),
            "Negative dimension accepted");
    invalid = straight;
    invalid.normals[0] = zVector(0, 0, 0);
    require(!zTsUtils::buildPrintMesh({{invalid}}, mesh, &error),
            "Zero normal accepted");
    invalid = straight;
    invalid.widths[0] = std::numeric_limits<float>::infinity();
    require(!zTsUtils::buildPrintMesh({{invalid}}, mesh, &error),
            "Nonfinite dimension accepted");
    invalid = straight;
    invalid.normals.assign(2, zVector(1, 0, 0));
    require(!zTsUtils::buildPrintMesh({{invalid}}, mesh, &error),
            "Parallel normal/tangent accepted");
    invalid = path({{0, 0, 0}, {1, 0, 1}, {0, 0, 2}}, {0, 1, 1, 2});
    require(!zTsUtils::buildPrintMesh({{invalid}}, mesh, &error) &&
                zFnMesh(mesh).numVertices() == 0,
            "Normal parallel to corner chord retained an invalid mesh");
    invalid =
        path({{0, 0, 0}, {1, 0, 0}, {2, 0, 0}, {1, 1, 0}}, {0, 1, 1, 2, 1, 3});
    require(!zTsUtils::buildPrintMesh({{invalid}}, mesh, &error),
            "Branched graph accepted");
    invalid = path({{0, 0, 0}, {0, 0, 0}}, {0, 1});
    require(!zTsUtils::buildPrintMesh({{invalid}}, mesh, &error),
            "Zero segment accepted");
    require(!zTsUtils::buildPrintMesh({}, mesh, &error),
            "Empty paths accepted");
    std::cout << "PASS print mesh one-sided rectangles, projected normals, "
                 "weighted corner chords, closures, separate "
                 "paths and invalid inputs\n";
    return 0;
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
