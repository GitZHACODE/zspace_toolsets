#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <zspace/zToolsets/fabrication/slicing/utils/Utils.h>

namespace zSpace {
namespace zTsUtils {

bool buildPrintMesh(const PrintPathArray &layers, zObjMesh &outMesh,
                    std::string *outError) {
  zFnMesh(outMesh).clear();
  if (outError)
    outError->clear();
  std::string location;
  try {
    zPointArray vertices;
    zIntArray counts, connects;
    const auto finite = [](const zVector &v) {
      return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
    };
    const auto unit = [](const zVector &v, const char *context) {
      const double length = std::hypot(double(v.x), double(v.y), double(v.z));
      if (!std::isfinite(length) || length <= 1e-12)
        throw std::runtime_error(std::string("Zero or invalid section direction: ")+context);
      return zVector(v.x / length, v.y / length, v.z / length);
    };
    const auto quad = [&](int a, int b, int c, int d) {
      counts.push_back(4);
      connects.insert(connects.end(), {a, b, c, d});
    };
    for (size_t l = 0; l < layers.size(); ++l)
      for (size_t p = 0; p < layers[l].size(); ++p) {
        location = "Layer " + std::to_string(l + 1) + ", path " +
                   std::to_string(p + 1) + ": ";
        const auto &path = layers[l][p];
        zPointArray points;
        zFnGraph(const_cast<zObjGraph &>(path.graph))
            .getVertexPositions(points);
        if (points.size() != path.widths.size() ||
            points.size() != path.heights.size() ||
            points.size() != path.normals.size())
          throw std::runtime_error("Missing per-point dimensions/normals.");
        for (size_t i = 0; i < points.size(); ++i)
          if (!finite(points[i]) || !finite(path.normals[i]) ||
              !std::isfinite(path.widths[i]) ||
              !std::isfinite(path.heights[i]) || path.widths[i] < 0 ||
              path.heights[i] < 0)
            throw std::runtime_error("Coordinates/normals must be finite and "
                                     "dimensions nonnegative.");
        std::vector<GraphPolyline> curves;
        if (!graphPolylines(path.graph, curves))
          throw std::runtime_error("Expected unbranched print curves.");
        size_t visited = 0;
        for (const auto &curve : curves) {
          const auto &ids = curve.Vertices;
          const size_t n = ids.size();
          if (n < (curve.Closed ? 3u : 2u))
            throw std::runtime_error("Insufficient curve vertices.");
          if (n >
              (size_t(std::numeric_limits<int>::max()) - vertices.size()) / 4)
            throw std::runtime_error("Print mesh exceeds the index range.");
          visited += n;
          const int base = static_cast<int>(vertices.size());
          const size_t segments = curve.Closed ? n : n - 1;
          std::vector<zVector> directions;
          for (size_t i = 0; i < segments; ++i)
            directions.push_back(
                unit(vecSub(points[ids[(i + 1) % n]], points[ids[i]]),"edge"));
          for (size_t i = 0; i < n; ++i) {
            const int id = ids[i];
            const zVector normal = unit(path.normals[id],"normal");
            const zVector incoming =
                directions[i == 0 ? (curve.Closed ? segments - 1 : 0) : i - 1];
            const zVector outgoing =
                directions[i < segments ? i : segments - 1];
            // Average unit edge directions at corners. A reversal uses the
            // outgoing edge; reject a path wholly parallel to its print normal.
            zVector tangent(incoming.x + outgoing.x, incoming.y + outgoing.y,
                            incoming.z + outgoing.z);
            zVector width = vecCross(normal, tangent);
            if (vecDot(width, width) <= 1e-24)
              width = vecCross(normal, outgoing);
            if (vecDot(width, width) <= 1e-24)
              width = vecCross(normal, incoming);
            width = unit(width,"cross-section width");
            const double w = path.widths[id] * 0.5, h = path.heights[id] * 0.5;
            for (const auto &sign :
                 {std::array<int, 2>{-1, -1}, {1, -1}, {1, 1}, {-1, 1}}) {
              const auto &center = points[id];
              zPoint corner(
                  center.x + sign[0] * w * width.x + sign[1] * h * normal.x,
                  center.y + sign[0] * w * width.y + sign[1] * h * normal.y,
                  center.z + sign[0] * w * width.z + sign[1] * h * normal.z);
              if (!finite(corner))
                throw std::runtime_error("Section coordinate overflow.");
              vertices.push_back(corner);
            }
          }
          for (size_t i = 0; i < segments; ++i) {
            const int a = base + int(4 * i), b = base + int(4 * ((i + 1) % n));
            for (int j = 0; j < 4; ++j)
              quad(a + j, a + (j + 1) % 4, b + (j + 1) % 4, b + j);
          }
          if (!curve.Closed) {
            quad(base + 3, base + 2, base + 1, base);
            const int end = base + int(4 * (n - 1));
            quad(end, end + 1, end + 2, end + 3);
          }
        }
        if (visited != points.size())
          throw std::runtime_error("Unused print vertices.");
      }
    if (vertices.empty())
      throw std::runtime_error("No print paths to loft.");
    zFnMesh(outMesh).create(vertices, counts, connects);
    return true;
  } catch (const std::exception &error) {
    zFnMesh(outMesh).clear();
    if (outError)
      *outError = location + error.what();
    return false;
  }
}
} // namespace zTsUtils
} // namespace zSpace
