#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <zspace/slicing.h>
using namespace zSpace;
constexpr int ProfileSides = 16;
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
    require(count == 4 || count == ProfileSides, "Expected side quads or ellipse end caps");
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
static void checkOutward(zObjMesh &mesh) {
  zPointArray points; zIntArray ids, counts; zFnMesh fn(mesh);
  fn.getVertexPositions(points); fn.getPolygonData(ids, counts);
  double volume = 0; size_t offset = 0;
  for (int count : counts) {
    const auto a = zTsUtils::toEigenVec(points[ids[offset]]);
    for (int j = 1; j + 1 < count; ++j) {
      const auto b = zTsUtils::toEigenVec(points[ids[offset+j]]);
      const auto c = zTsUtils::toEigenVec(points[ids[offset+j+1]]);
      volume += a.dot(b.cross(c)) / 6;
    }
    offset += count;
  }
  require(volume > 0, "Ellipse loft is not outward wound");
}
static void checkPoint(const zPoint& actual, const zPoint& expected, const char* message) {
  require(zTsUtils::vecDistSq(actual, expected) < 1e-10, message);
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
    checkClosed(mesh, 32, 18);
    checkOutward(mesh);
    zPointArray points;
    zFnMesh(mesh).getVertexPositions(points);
    for (int ring = 0; ring < 2; ++ring) {
      const int base = ring * ProfileSides;
      const double x = ring * 10, width = ring ? 3 : 1, height = ring ? 2 : 4;
      checkPoint(points[base], zPoint(x, 0, 0), "TCP is not the exact ellipse top");
      checkPoint(points[base+4], zPoint(x, -width, -height/2), "Ellipse left extent/center mismatch");
      checkPoint(points[base+8], zPoint(x, 0, -height), "Ellipse bottom is not TCP minus height");
      checkPoint(points[base+12], zPoint(x, width, -height/2), "Ellipse right extent/center mismatch");
      for (int j = 0; j < ProfileSides; ++j) {
        const auto& q = points[base+j];
        const double u = q.y/width, v = (q.z+height/2)/(height/2);
        require(std::abs(u*u+v*v-1) < 2e-6, "A section vertex is not on the ellipse");
        require(q.z <= 0 && q.z >= -height, "Ellipse exceeds the TCP/bottom interval");
      }
      require(points[base+2].y > -width && points[base+2].y < 0 &&
              points[base+2].z > -height/2 && points[base+2].z < 0,
              "Intermediate ellipse profile is still a rectangle");
    }
    // Width and height axes rotate with the supplied normal, rather than world
    // Z.
    auto tilted = straight;
    tilted.normals.assign(2, zVector(0, 1, 0));
    require(zTsUtils::buildPrintMesh({{tilted}}, mesh, &error), error.c_str());
    zFnMesh(mesh).getVertexPositions(points);
    checkPoint(points[0], zPoint(0, 0, 0), "Rotated ellipse moved the TCP");
    checkPoint(points[4], zPoint(0, -2, 1), "Rotated ellipse width/center mismatch");
    checkPoint(points[8], zPoint(0, -4, 0), "Rotated ellipse ignored the normal");
    checkClosed(mesh, 32, 18); checkOutward(mesh);
    auto downward = straight;
    downward.normals.assign(2, zVector(0, 0, -2));
    require(zTsUtils::buildPrintMesh({{downward}}, mesh, &error), error.c_str());
    zFnMesh(mesh).getVertexPositions(points);
    checkPoint(points[0], zPoint(0, 0, 0), "Negative growth normal moved TCP");
    checkPoint(points[8], zPoint(0, 0, 4), "Negative growth normal was flipped");
    checkClosed(mesh, 32, 18); checkOutward(mesh);
    // A sloping path requires projecting the normal into the section plane.
    auto slope = path({{0, 0, 0}, {10, 0, 10}}, {0, 1});
    require(zTsUtils::buildPrintMesh({{slope}}, mesh, &error), error.c_str());
    zFnMesh(mesh).getVertexPositions(points);
    const double rise = 2 * std::sqrt(2.0);
    checkPoint(points[0], zPoint(0, 0, 0), "Sloping ellipse moved TCP");
    checkPoint(points[4], zPoint(rise/2, -1, -rise/2),
               "Sloping ellipse width/height center mismatch");
    checkPoint(points[8], zPoint(rise, 0, -rise),
               "Section is not perpendicular to the path with projected height");
    checkClosed(mesh, 32, 18); checkOutward(mesh);
    // A 3:4 corner follows the neighbor chord, not the unit-edge bisector.
    auto corner = path({{0, 0, 0}, {3, 0, 0}, {3, 4, 0}},
                       {0, 1, 1, 2});
    require(zTsUtils::buildPrintMesh({{corner}}, mesh, &error), error.c_str());
    zFnMesh(mesh).getVertexPositions(points);
    checkPoint(points[16], zPoint(3, 0, 0), "Corner ellipse moved TCP");
    checkPoint(points[20], zPoint(3.8, -0.6, -2), "Corner chord left weighting mismatch");
    checkPoint(points[28], zPoint(2.2, 0.6, -2), "Corner chord right weighting mismatch");
    checkPoint(points[24], zPoint(3, 0, -4), "Corner height axis mismatch");
    checkClosed(mesh, 48, 34); checkOutward(mesh);
    auto loop = path({{0, 0, 0}, {10, 0, 0}, {10, 10, 0}, {0, 10, 0}},
                     {0, 1, 1, 2, 2, 3, 3, 0});
    require(
        zTsUtils::buildPrintMesh({{loop}, {straight, straight}}, mesh, &error),
        error.c_str());
    checkClosed(mesh, 128, 100); // loop:64 quads; independent open paths:18 faces each
    zIntArray ids, counts;
    zFnMesh(mesh).getPolygonData(ids, counts);
    size_t offset = 0;
    for (int count : counts) {
      const auto component = [](int id) { return id < 64 ? 0 : id < 96 ? 1 : 2; };
      for (int j = 1; j < count; ++j)
        require(component(ids[offset+j]) == component(ids[offset]),
                "Loft bridged separate components or layers");
      offset += count;
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
    std::cout << "PASS print mesh TCP-top ellipses, projected normals, "
                 "weighted corner chords, closures, separate "
                 "paths and invalid inputs\n";
    return 0;
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
