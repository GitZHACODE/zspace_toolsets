#include <zspace/zToolsets/fabrication/slicing/utils/IO.h>
// Toolsets owns the slicer schema; Core owns generic USD operations.
#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <zspace/zIO/zUSDapi.h>
#include <zspace/zToolsets/fabrication/slicing/utils/Utils.h>
namespace zSpace::zTsIO {
namespace {
inline std::string formatPrintPathUSDA(const std::string &text) {
  constexpr size_t maxLineLength = 500;
  std::istringstream input(text);
  std::ostringstream output;
  std::string line;
  while (std::getline(input, line)) {
    const auto first = line.find_first_not_of(" \t");
    const std::string declaration =
        first == std::string::npos ? "" : line.substr(first);
    const bool statisticArray =
        declaration.find("double[] printheightdomain = ") == 0 ||
        declaration.find("double[] printwidthdomain = ") == 0;
    const bool statisticScalar =
        declaration.find("double printlength = ") == 0 ||
        declaration.find("double printvolume = ") == 0;
    if (statisticArray || statisticScalar) {
      const auto equals = line.find("= ");
      std::string values = line.substr(equals + 2);
      for (char &c : values)
        if (c == '[' || c == ']' || c == ',')
          c = ' ';
      std::istringstream numbers(values);
      double value;
      output << line.substr(0, equals + 2) << std::fixed
             << std::setprecision(2);
      if (statisticArray)
        output << '[';
      bool firstValue = true;
      while (numbers >> value) {
        if (!firstValue)
          output << ", ";
        output << value;
        firstValue = false;
      }
      if (statisticArray)
        output << ']';
      output << '\n';
      continue;
    }
    const bool vectors =
        declaration.find("custom point3f[] points = [") == 0 ||
        declaration.find("custom normal3f[] printNormals = [") == 0;
    const bool scalars =
        declaration.find("custom float[] printWidths = [") == 0 ||
        declaration.find("custom float[] printHeight = [") == 0;
    const auto begin = line.find("= [");
    const auto end = line.rfind(']');
    if ((!vectors && !scalars) || begin == std::string::npos ||
        end == std::string::npos || end <= begin + 3) {
      output << line << '\n';
      continue;
    }
    const std::string indent(first, ' ');
    const std::string itemIndent = indent + "    ";
    output << line.substr(0, begin + 3) << '\n' << itemIndent;
    size_t column = itemIndent.size();
    size_t itemStart = begin + 3;
    int depth = 0;
    for (size_t i = itemStart; i <= end; ++i) {
      const char c = line[i];
      if (c == '(')
        ++depth;
      if (c == ')')
        --depth;
      if (i != end && (c != ',' || depth != 0))
        continue;
      std::string item = line.substr(itemStart, i - itemStart);
      itemStart = i + 1;
      const auto start = item.find_first_not_of(" \t\r\n");
      if (start == std::string::npos)
        continue;
      item = item.substr(start, item.find_last_not_of(" \t\r\n") - start + 1);
      if (i != end)
        item += ',';
      if (column > itemIndent.size()) {
        if (column + 1 + item.size() > maxLineLength) {
          output << '\n' << itemIndent;
          column = itemIndent.size();
        } else {
          output << ' ';
          ++column;
        }
      }
      output << item;
      column += item.size();
    }
    output << '\n' << indent << line.substr(end) << '\n';
  }
  return output.str();
}

inline void requireType(USDdoc &doc, const std::string &path,
                        const char *type) {
  if (doc.getPrimAtPath(path).getTypeName() != type)
    throw std::runtime_error("Expected " + std::string(type) + " at " + path);
}
// Scale lengths; ratios and directions use scale=1.
double format(double value, double scale = 1.0) {
  const double result = value * scale;
  if (!std::isfinite(result))
    throw std::runtime_error("Numeric conversion overflow.");
  return result;
}
// Import Y-up to Z-up; export applies the inverse rotation.
std::array<double, 3> format(const std::array<double, 3> &value, double scale,
                             bool yUp, bool toInternal = true) {
  std::array<double, 3> result{format(value[0], scale), format(value[1], scale),
                               format(value[2], scale)};
  if (yUp) {
    const double y = result[1], z = result[2];
    result[1] = toInternal ? -z : z;
    result[2] = toInternal ? y : -y;
  }
  return result;
}
void readSettings(USDdoc &doc, double scale, zFloatArray &boundary,
                  zFloatArray &bracing, SlicerUSDInput &out) {
  out.Settings = defaultSlicerSettings();
  out.PrintSettings = defaultPrintSynSettings();
  auto &s = out.Settings;
  auto &p = out.PrintSettings;
  // Defaults are already in mm; only authored values are converted.
  boundary = s.BoundaryWidth;
  bracing = {s.BracingWidth};
  auto settingsPrim = doc.getPrimAtPath("/Slicer_setting");
  for (const auto &name : doc.readPrim("/Slicer_setting").attributes) {
    auto attribute = settingsPrim.getAttribute(name);
    if (name == "slicer:BoundaryWidth") {
      auto values = attribute.get<std::vector<float>>();
      if (!values.empty()) {
        for (auto &v : values)
          v = static_cast<float>(format(v, scale));
        boundary = std::move(values);
      }
    } else if (name == "slicer:BracingWidth") {
      auto values = attribute.get<std::vector<float>>();
      if (!values.empty()) {
        for (auto &v : values)
          v = static_cast<float>(format(v, scale));
        bracing = std::move(values);
      }
    } else if (name == "slicer:DoubleBoundary")
      s.DoubleBoundary = attribute.get<bool>();
    else if (name == "slicer:OverlapWidth")
      s.OverlapWidth = format(attribute.get<double>(), scale);
    else if (name == "slicer:SlotOffset")
      s.SlotOffset = format(attribute.get<double>());
    else if (name == "slicer:BracingSlot")
      s.BracingSlot = format(attribute.get<double>());
    else if (name == "slicer:FieldScaleUp")
      s.FieldScaleUp = format(attribute.get<double>());
    else if (name == "slicer:FieldCellSize")
      s.FieldCellSize = format(attribute.get<double>(), scale);
    else if (name == "slicer:LayerSpacing")
      s.LayerSpacing = format(attribute.get<double>(), scale);
    else if (name == "slicer:_HTraceStepFraction")
      s._HTraceStepFraction = format(attribute.get<double>());
    else if (name == "slicer:_HTraceMaxSteps")
      s._HTraceMaxSteps = attribute.get<int>();
    else if (name == "slicer:_HEqualiseIterations")
      s._HEqualiseIterations = attribute.get<int>();
    else if (name == "slicer:_HArapIterations")
      s._HArapIterations = attribute.get<int>();
    else if (name == "slicer:_HUnrollStart")
      s._HUnrollStart = attribute.get<int>();
    else if (name == "printSyn:ResampleLength")
      p.ResampleLength = format(attribute.get<double>(), scale);
    else if (name == "printSyn:FeatureAngleTol")
      p.FeatureAngleTol = format(attribute.get<double>());
    else if (name == "printSyn:SeamOffset")
      p.SeamOffset = format(attribute.get<double>(), scale);
    else if (name == "printSyn:SeamSize")
      p.SeamSize = format(attribute.get<double>(), scale);
    else if (name.rfind("slicer:", 0) == 0 || name.rfind("printSyn:", 0) == 0)
      out.UnknownSettings.push_back(name);
  }
  s.BoundaryWidth = boundary;
  s.BracingWidth = bracing.front();
}
inline void readBlock(USDdoc &doc, const std::string &blockPath, double scale,
                      bool yUp, const zFloatArray &boundary,
                      SlicerUSDInput &out) {
  const auto coordinates = doc.readWorldPoints(blockPath);
  zPointArray points;
  points.reserve(coordinates.size());
  for (const auto &point : coordinates) {
    const auto p = format(point, scale, yUp);
    points.emplace_back(p[0], p[1], p[2]);
  }
  auto counts = doc.getPrimAtPath(blockPath).getAttribute("faceVertexCounts").get<std::vector<int>>(),
       indices = doc.getPrimAtPath(blockPath).getAttribute("faceVertexIndices").get<std::vector<int>>();
  if (counts.empty())
    throw std::runtime_error("Block contains no faces.");
  size_t total = 0;
  for (int c : counts) {
    if (c < 3)
      throw std::runtime_error("Invalid polygon size.");
    total += c;
  }
  if (total != indices.size())
    throw std::runtime_error("Face counts do not match indices.");
  for (int id : indices)
    if (id < 0 || size_t(id) >= points.size())
      throw std::runtime_error("Invalid vertex index.");
  out.FaceTags = doc.getPrimAtPath("/Slicer_data").getAttribute("block:FaceTag").get<std::vector<int>>();
  out.FaceWidthIds = doc.getPrimAtPath("/Slicer_data").getAttribute("block:widthid").get<std::vector<int>>();
  if (out.FaceTags.size() != counts.size() ||
      out.FaceWidthIds.size() != counts.size())
    throw std::runtime_error("Face attribute count mismatch.");
  for (size_t i = 0; i < counts.size(); ++i) {
    int t = out.FaceTags[i], w = out.FaceWidthIds[i];
    if (t < 1 || t > 6 ||
        (isSideFaceTag(t) ? (w < 0 || size_t(w) >= boundary.size()) : w != -1))
      throw std::runtime_error("Invalid face tag or widthid.");
  }
  const bool reverse = (doc.getPrimAtPath(blockPath).getAttribute("orientation").get<std::string>() ==
                        "leftHanded") != doc.readPrim(blockPath).mirrored;
  if (reverse) {
    size_t start = 0;
    for (int c : counts) {
      std::reverse(indices.begin() + start, indices.begin() + start + c);
      start += c;
    }
  }
  zFnMesh(out.Block).create(points, counts, indices);
  if (zFnMesh(out.Block).numVertices() != points.size() ||
      zFnMesh(out.Block).numPolygons() != counts.size())
    throw std::runtime_error("Mesh construction changed indexing.");
}
inline void readBracing(USDdoc &doc, const std::string &bracingPath,
                        double scale, bool yUp, const zFloatArray &bracing,
                        SlicerUSDInput &out) {
  if (doc.getPrimAtPath(bracingPath).getAttribute("type").get<std::string>() != "linear")
    throw std::runtime_error("Bracing requires linear curves.");
  const auto wrap = doc.getPrimAtPath(bracingPath).getAttribute("wrap").get<std::string>();
  if (wrap != "periodic" && wrap != "nonperiodic")
    throw std::runtime_error("Invalid curve wrap.");
  const bool closed = wrap == "periodic";
  const auto coordinates = doc.readWorldPoints(bracingPath);
  zPointArray points;
  points.reserve(coordinates.size());
  for (const auto &point : coordinates) {
    const auto p = format(point, scale, yUp);
    points.emplace_back(p[0], p[1], p[2]);
  }
  auto counts = doc.getPrimAtPath(bracingPath).getAttribute("curveVertexCounts").get<std::vector<int>>();
  out.BracingWidthIds = doc.getPrimAtPath("/Slicer_data").getAttribute("bracing:widthid").get<std::vector<int>>();
  size_t total = 0;
  for (int c : counts) {
    if (c < (closed ? 3 : 2))
      throw std::runtime_error("Curve has too few points.");
    total += c;
  }
  if (counts.empty() || total != points.size() ||
      out.BracingWidthIds.size() != counts.size())
    throw std::runtime_error("Curve array count mismatch.");
  out.BracingCurveCount = counts.size();
  size_t start = 0;
  zIntArray edges;
  for (size_t i = 0; i < counts.size(); ++i) {
    int n = counts[i], id = out.BracingWidthIds[i];
    if (id < 0 || size_t(id) >= bracing.size())
      throw std::runtime_error("Invalid bracing widthid.");

    for (int j = 0; j < (closed ? n : n - 1); ++j) {
      int next = (j + 1) % n;
      if (zTsUtils::vecDistSq(points[start + j], points[start + next]) <= 0)
        throw std::runtime_error("Zero-length bracing segment.");
      edges.push_back(static_cast<int>(start) + j);
      edges.push_back(static_cast<int>(start) + next);
    }
    out.BracingPrintWidths.push_back(bracing[id]);
    start += n;
  }
  zFnGraph(out.Bracing).create(points, edges);
}

inline void writePrintPath(USDdoc &doc, const PrintPath &path, double scale,
                                    bool yUp, int layerCount, int blockid) {
  if (layerCount <= 0)
    throw std::runtime_error("Layer count must be positive.");
  std::vector<zTsUtils::GraphPolyline> curves;
  if (!zTsUtils::graphPolylines(path.graph, curves) || curves.size() != 1)
    throw std::runtime_error(
        "Print path must be one connected, unbranched curve.");
  zPointArray positions;
  zFnGraph(const_cast<zObjGraph &>(path.graph)).getVertexPositions(positions);
  if (path.normals.size() != positions.size() ||
      path.widths.size() != positions.size() ||
      path.heights.size() != positions.size() ||
      (!path.faceTags.empty() && path.faceTags.size() != positions.size()))
    throw std::runtime_error("Path attribute count mismatch.");
  auto ids = curves[0].Vertices;
  if (ids.size() > 1 && ids.front() == ids.back())
    ids.pop_back();
  if (ids.size() != positions.size() ||
      ids.size() < (curves[0].Closed ? 3u : 2u))
    throw std::runtime_error("Invalid path vertices.");
  std::vector<float> points, normals, widths, heights;
  std::vector<int> faceTags;
  points.reserve(ids.size() * 3);
  normals.reserve(ids.size() * 3);
  widths.reserve(ids.size());
  heights.reserve(ids.size());
  const auto f = [](double v) {
    if (!std::isfinite(v) || std::abs(v) > std::numeric_limits<float>::max())
      throw std::runtime_error("Export float overflow.");
    return static_cast<float>(v);
  };
  double length = 0, volume = 0;
  for (size_t i = 0; i < ids.size(); ++i) {
    int id = ids[i];
    if (id < 0 || size_t(id) >= positions.size())
      throw std::runtime_error("Invalid path index.");
    const auto &p = positions[id];
    const auto &n = path.normals[id];
    for (double v :
         {double(p.x), double(p.y), double(p.z), double(n.x), double(n.y),
          double(n.z), double(path.widths[id]), double(path.heights[id])})
      if (!std::isfinite(v))
        throw std::runtime_error("Non-finite path data.");
    if (path.widths[id] < 0 || path.heights[id] < 0)
      throw std::runtime_error("Negative print width/height.");
    const auto point =
        format(std::array<double, 3>{p.x, p.y, p.z}, scale, yUp, false);
    const auto normal =
        format(std::array<double, 3>{n.x, n.y, n.z}, 1.0, yUp, false);
    for (double v : point)
      points.push_back(f(v));
    for (double v : normal)
      normals.push_back(f(v));
    widths.push_back(f(format(path.widths[id], scale)));
    heights.push_back(f(format(path.heights[id], scale)));
    if (!path.faceTags.empty()) faceTags.push_back(path.faceTags[id]);
    if (i + 1 < ids.size() || curves[0].Closed) {
      double d = std::sqrt(zTsUtils::vecDistSq(
                     p, positions[ids[(i + 1) % ids.size()]])) *
                 scale;
      length += d;
      volume += path.widths[id] * scale * path.heights[id] * scale * d;
    }
  }
  const auto round = [](double v) {
    if (!std::isfinite(v))
      throw std::runtime_error("Export overflow.");
    return std::round(v * 100) / 100;
  };
  auto domain = [&](const auto &a) {
    double lo = format(double(a[0]), scale), hi = lo;
    for (const auto &v : a) {
      lo = std::min(lo, format(double(v), scale));
      hi = std::max(hi, format(double(v), scale));
    }
    return std::vector<double>{round(lo), round(hi)};
  };
  for (const auto *name : {"/TDP_gcode", "/Gcode", "/printpath"}) doc.removePrim(name);
  auto prim = doc.definePrim("/TDP_gcode", "Scope");
  prim.createAttribute("points", "point3f[]").set(points);
  prim.createAttribute("printNormals", "normal3f[]").set(normals);
  prim.createAttribute("printWidths", "float[]").set(widths);
  prim.createAttribute("printHeight", "float[]").set(heights);
  if (!faceTags.empty()) prim.createAttribute("printFaceTags", "int[]").set(faceTags);
  auto metadata = doc.readCustomLayerData();
  metadata.erase("printstatisticsvalid");
  metadata["blockid"] = {"int", blockid};
  metadata["layercount"] = {"int", layerCount};
  metadata["printheightdomain"] = {"double[]", domain(path.heights)};
  metadata["printwidthdomain"] = {"double[]", domain(path.widths)};
  metadata["printlength"] = {"double", round(length)};
  metadata["printvolume"] = {"double", round(volume)};
  doc.setCustomLayerData(metadata);
}

// Geometry follows the shared print-mesh builder; USD authoring stays in Core IO.
void writePrintGeometry(USDdoc &doc, const PrintPathArray &layers,
                        double scale, bool yUp) {
  zObjMesh mesh;
  std::string error;
  if (!zTsUtils::buildPrintMesh(layers, mesh, &error))
    throw std::runtime_error("Print mesh: " + error);
  const auto scalar = [](double v) {
    if (!std::isfinite(v) || std::abs(v) > std::numeric_limits<float>::max())
      throw std::runtime_error("Print geometry conversion overflow.");
    return static_cast<float>(v);
  };
  auto append = [&](std::vector<float> &out, const zPoint &p) {
    const auto v = format(std::array<double,3>{p.x,p.y,p.z}, scale, yUp, false);
    for (double x : v) out.push_back(scalar(x));
  };
  zPointArray positions;
  zIntArray counts, connects;
  zFnMesh fn(mesh);
  fn.getVertexPositions(positions);
  fn.getPolygonData(connects, counts);
  std::vector<float> points;
  for (const auto &p : positions) append(points, p);
  for (const char *name : {"/printMesh", "/TDP_printmesh", "/TDP_printpath"})
    doc.removePrim(name);
  auto prim = doc.definePrim("/TDP_printmesh", "Mesh");
  prim.createAttribute("points", "point3f[]", false).set(points);
  prim.createAttribute("faceVertexCounts", "int[]", false).set(counts);
  prim.createAttribute("faceVertexIndices", "int[]", false).set(connects);
  prim.createAttribute("subdivisionScheme", "token", false).set(std::string("none"));
  points.clear(); counts.clear();
  std::vector<float> widths;
  for (const auto &layer : layers) for (const auto &path : layer) {
    positions.clear();
    zFnGraph(const_cast<zObjGraph &>(path.graph)).getVertexPositions(positions);
    std::vector<zTsUtils::GraphPolyline> components;
    if (!zTsUtils::graphPolylines(path.graph, components))
      throw std::runtime_error("Invalid renderable print path.");
    for (const auto &component : components) {
      auto ids = component.Vertices;
      if (component.Closed) ids.push_back(ids.front());
      counts.push_back(static_cast<int>(ids.size()));
      for (int id : ids) {
        append(points, positions[id]);
        widths.push_back(scalar(path.widths[id] * scale));
      }
    }
  }
  auto curves = doc.definePrim("/TDP_printpath", "BasisCurves");
  curves.createAttribute("type", "token", false).set(std::string("linear"));
  curves.createAttribute("wrap", "token", false).set(std::string("nonperiodic"));
  curves.createAttribute("points", "point3f[]", false).set(points);
  curves.createAttribute("curveVertexCounts", "int[]", false).set(counts);
  zUSDEdits edits;
  zUSDPrimEdit edit;
  edit.path = "/TDP_printpath"; edit.type = "BasisCurves";
  USDdoc::writeFloatArray(edit, "widths", widths, "vertex", false);
  edits.prims.push_back(std::move(edit));
  doc.applyEdits(edits);
}

} // namespace

bool importUSD(const std::string &filePath, SlicerUSDInput &out,
                        std::string *outError, const std::string &blockPath,
                        const std::string &bracingPath,
                        const std::string &startPath) {
  out = SlicerUSDInput();
  if (outError)
    outError->clear();
  try {
    USDdoc doc;
    doc.open(filePath);
    const auto metadata = doc.readMetadata();
    const double scale = format(metadata.metersPerUnit, 1000.0);
    if (scale <= 0 || (metadata.upAxis != "Y" && metadata.upAxis != "Z"))
      throw std::runtime_error("Invalid USD units or up axis.");
    const bool yUp = metadata.upAxis == "Y";
    out.MetersPerUnit = metadata.metersPerUnit;
    out.UpAxis = metadata.upAxis;
    requireType(doc, blockPath, "Mesh");
    requireType(doc, bracingPath, "BasisCurves");
    requireType(doc, startPath, "Points");
    requireType(doc, "/Slicer_setting", "Scope");
    requireType(doc, "/Slicer_data", "Scope");
    zFloatArray boundaryWidths, bracingWidths;
    readSettings(doc, scale, boundaryWidths, bracingWidths, out);
    readBlock(doc, blockPath, scale, yUp, boundaryWidths, out);
    readBracing(doc, bracingPath, scale, yUp, bracingWidths, out);
    const auto coordinates = doc.readWorldPoints(startPath);
    if (coordinates.size() != 1)
      throw std::runtime_error("Slot requires exactly one point.");
    const auto anchor = format(coordinates.front(), scale, yUp);
    out.StartPoint = zPoint(anchor[0], anchor[1], anchor[2]);
    return true;
  } catch (const std::exception &e) {
    out = SlicerUSDInput();
    if (outError)
      *outError = e.what();
    return false;
  }
}

bool exportUSD(const PrintPath &path, const std::string &input,
                        const std::string &output, int layerCount,
                        std::string *error, int blockid, const PrintPathArray *layerPaths) {
  const char *phase = "open source document";
  if (error)
    error->clear();
  try {
    USDdoc doc;
    doc.open(input);
    const auto metadata = doc.readMetadata();
    if (!std::isfinite(metadata.metersPerUnit) || metadata.metersPerUnit <= 0 ||
        (metadata.upAxis != "Y" && metadata.upAxis != "Z"))
      throw std::runtime_error("Invalid USD units or up axis.");
    const double scale = format(0.001 / metadata.metersPerUnit);
    if (scale <= 0)
      throw std::runtime_error("Invalid USD scale.");
    const bool yUp = metadata.upAxis == "Y";
    requireType(doc, "/AG_block", "Mesh");
    requireType(doc, "/Slicer_bracing", "BasisCurves");
    requireType(doc, "/Slicer_slot", "Points");
    phase = "prepare print path";
    writePrintPath(doc, path, scale, yUp, layerCount, blockid);
    phase = "prepare print geometry";
    if (layerPaths && layerPaths->size() != size_t(layerCount))
      throw std::runtime_error("Print mesh layers must match the actual layer count.");
    PrintPathArray fallback;
    if (!layerPaths) fallback = {{path}};
    writePrintGeometry(doc, layerPaths ? *layerPaths : fallback, scale, yUp);
    phase = "save output document";
    const auto dot = output.find_last_of('.');
    auto suffix = dot == std::string::npos ? std::string() : output.substr(dot);
    std::transform(suffix.begin(), suffix.end(), suffix.begin(),
                   [](unsigned char c) { return char(std::tolower(c)); });
    if (suffix == ".usda") {
      // Complete serialization before truncating the destination (which may be
      // the input).
      const auto text = formatPrintPathUSDA(doc.text());
      std::ofstream file(output, std::ios::binary | std::ios::trunc);
      file << text;
      file.close();
      if (!file)
        throw std::runtime_error("Cannot save USDA.");
    } else if (suffix == ".usd" || suffix == ".usdc")
      doc.save(output);
    else
      throw std::runtime_error("Expected .usd, .usda or .usdc.");
    return true;
  } catch (const std::exception &e) {
    if (error)
      *error = std::string(phase) + ": " + e.what();
    return false;
  }
}

} // namespace zSpace::zTsIO
