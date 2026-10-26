// This file is part of zspace, a simple C++ collection of geometry
// data-structures & algorithms, data analysis & visualization framework.
//
// Copyright (C) 2019 ZSPACE
//
// This Source Code Form is subject to the terms of the MIT License
// If a copy of the MIT License was not distributed with this file, You can
// obtain one at https://opensource.org/licenses/MIT.
//
// Author : Weijie Tao <weijie.tao@zspace.com>
//

#ifndef EIGEN_DONT_PARALLELIZE
#define EIGEN_DONT_PARALLELIZE
#endif

// Parse USD/Windows SDK headers before legacy zSpace global namespace imports.
#if defined(ZSPACE_TS3DP_OPENUSD)
#include <pxr/usd/usd/stage.h>
#include <pxr/usd/sdf/layer.h>
#include <pxr/usd/usdGeom/tokens.h>
#include <pxr/usd/usdGeom/mesh.h>
#include <pxr/usd/usdGeom/basisCurves.h>
#include <pxr/usd/usdGeom/points.h>
#include <pxr/usd/usdGeom/scope.h>
#include <pxr/usd/usdGeom/primvarsAPI.h>
#include <pxr/usd/usdGeom/xformCache.h>
#include <pxr/usd/usdGeom/metrics.h>
#include <pxr/base/vt/dictionary.h>
#include <pxr/base/vt/array.h>
#include <nlohmann/json.hpp>
#include <cmath>
#include <stdexcept>
#endif

#include <zspace/zToolsets/fabrication/slicing/utils/Utils.h>

#ifndef ZSPACE_INLINE
#define ZSPACE_INLINE ZSPACE_TOOLSETS_INLINE
#endif

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>



namespace zSpace {
namespace zTsIO {

#if defined(ZSPACE_TS3DP_OPENUSD)
namespace {
// Only reflow the numeric arrays authored by exportPrintPathUSD. Other USD
// syntax, including strings, comments and attribute metadata, stays untouched.
std::string formatPrintPathUSDA(const std::string &text) {
  constexpr size_t maxLineLength = 500;
  std::istringstream input(text);
  std::ostringstream output;
  std::string line;
  while (std::getline(input, line)) {
    const auto first = line.find_first_not_of(" \t");
    const std::string declaration = first == std::string::npos ? "" : line.substr(first);
    const bool statisticArray = declaration.find("double[] printheightdomain = ") == 0 ||
                                declaration.find("double[] printwidthdomain = ") == 0;
    const bool statisticScalar = declaration.find("double printlength = ") == 0 ||
                                 declaration.find("double printvolume = ") == 0;
    if (statisticArray || statisticScalar) {
      const auto equals = line.find("= ");
      std::string values = line.substr(equals + 2);
      for (char &c : values) if (c == '[' || c == ']' || c == ',') c = ' ';
      std::istringstream numbers(values);
      double value;
      output << line.substr(0, equals + 2);
      // Preserve reusable domain precision, including metre-based stages.
      if(declaration.find("double[] printheightdomain = ")==0)output<<std::defaultfloat<<std::setprecision(17);
      else output<<std::fixed<<std::setprecision(2);
      if (statisticArray) output << '[';
      bool firstValue = true;
      while (numbers >> value) {
        if (!firstValue) output << ", ";
        output << value;
        firstValue = false;
      }
      if (statisticArray) output << ']';
      output << '\n';
      continue;
    }
    const bool vectors = declaration.find("custom point3f[] points = [") == 0 ||
                         declaration.find("custom normal3f[] printNormals = [") == 0;
    const bool scalars = declaration.find("custom float[] printWidths = [") == 0 ||
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
      if (c == '(') ++depth;
      if (c == ')') --depth;
      if (i != end && (c != ',' || depth != 0)) continue;
      std::string item = line.substr(itemStart, i - itemStart);
      itemStart = i + 1;
      const auto start = item.find_first_not_of(" \t\r\n");
      if (start == std::string::npos) continue;
      item = item.substr(start, item.find_last_not_of(" \t\r\n") - start + 1);
      if (i != end) item += ',';
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
}
#endif

ZSPACE_INLINE bool exportPrintPathUSD(const PrintPath &path,
    const std::string &inputFilePath, const std::string &filePath, int layerCount, std::string *outError, int blockid,
    const PrintPathArray *layerPaths) {
  if (outError) outError->clear();
  const auto fail = [&](const std::string &message) {
    if (outError) *outError = message;
    return false;
  };
#if !defined(ZSPACE_TS3DP_OPENUSD)
  return fail("USD export requires ZSPACE_TS3DP_OPENUSD and linked OpenUSD libraries.");
#else
  try {
    if (layerCount <= 0) return fail("Export requires a positive actual layer count.");
    if (layerPaths && layerPaths->size() != size_t(layerCount))
      return fail("Print mesh layers must match the actual layer count.");
    zObjMesh printMesh;
    std::string meshError;
    PrintPathArray fallback;
    if (!layerPaths) fallback = {{path}};
    const auto &renderPaths = layerPaths ? *layerPaths : fallback;
    if (!zTsUtils::buildPrintMesh(renderPaths, printMesh, &meshError))
      return fail("Print mesh: " + meshError);
    const auto dot = filePath.find_last_of('.');
    std::string extension = dot == std::string::npos ? "" : filePath.substr(dot);
    std::transform(extension.begin(), extension.end(), extension.begin(),
                   [](unsigned char c) { return char(std::tolower(c)); });
    if (extension != ".usd" && extension != ".usda" && extension != ".usdc")
      return fail("Use a .usd, .usda or .usdc output file.");
    std::vector<zTsUtils::GraphPolyline> curves;
    if (!zTsUtils::graphPolylines(path.graph, curves) || curves.size() != 1)
      return fail("Print path must be one connected, unbranched curve.");
    zPointArray positions;
    zFnGraph(const_cast<zObjGraph &>(path.graph)).getVertexPositions(positions);
    if (path.normals.size() != positions.size() || path.widths.size() != positions.size() ||
        path.heights.size() != positions.size())
      return fail("Normals, widths and heights must have one value per graph vertex.");
    auto ids = curves.front().Vertices;
    if (ids.size() > 1 && ids.front() == ids.back()) ids.pop_back();
    if (ids.size() < (curves.front().Closed ? 3u : 2u) || ids.size() != positions.size())
      return fail("Print path contains insufficient or unused vertices.");
    pxr::VtArray<pxr::GfVec3f> points, normals;
    pxr::VtArray<float> widths, heights;
    for (int id : ids) {
      if (id < 0 || size_t(id) >= positions.size()) return fail("Invalid path vertex index.");
      const auto &p = positions[id];
      const auto &n = path.normals[id];
      if (!std::isfinite(p.x) || !std::isfinite(p.y) || !std::isfinite(p.z) ||
          !std::isfinite(n.x) || !std::isfinite(n.y) || !std::isfinite(n.z) ||
          !std::isfinite(path.widths[id]) || !std::isfinite(path.heights[id]) ||
          path.widths[id] < 0 || path.heights[id] < 0)
        return fail("Path attributes must be finite; widths and heights must be nonnegative.");
      points.push_back(pxr::GfVec3f(p.x, p.y, p.z));
      normals.push_back(pxr::GfVec3f(n.x, n.y, n.z));
      widths.push_back(path.widths[id]);
      heights.push_back(path.heights[id]);
    }
    const auto heightRange = std::minmax_element(heights.begin(), heights.end());
    const auto widthRange = std::minmax_element(widths.begin(), widths.end());
    double printLength = 0.0, printVolume = 0.0;
    const size_t segmentCount = curves.front().Closed ? points.size() : points.size() - 1;
    for (size_t i = 0; i < segmentCount; ++i) {
      const auto &a = points[i];
      const auto &b = points[(i + 1) % points.size()];
      const double distance = std::hypot(double(b[0]) - a[0],
          double(b[1]) - a[1], double(b[2]) - a[2]);
      printLength += distance;
      printVolume += double(heights[i]) * double(widths[i]) * distance;
    }
    if (!std::isfinite(printLength) || !std::isfinite(printVolume))
      return fail("Print statistics overflow.");
    // Load a fresh private layer so exporting over the input does not mutate
    // another live stage or reuse stale cached file contents.
    auto layer = pxr::SdfLayer::OpenAsAnonymous(inputFilePath);
    if (!layer) return fail("Could not open input USD file.");
    auto stage = pxr::UsdStage::Open(layer);
    if (!stage) return fail("Could not open input USD stage.");
    // Preserve the source stage units while exporting internally millimetric paths.
    const double sourceUnit = pxr::UsdGeomGetStageMetersPerUnit(stage);
    if (!std::isfinite(sourceUnit) || sourceUnit <= 0) return fail("Invalid stage metersPerUnit.");
    const double fromMillimetres = 0.001 / sourceUnit;
    for (auto &p : points) p *= static_cast<float>(fromMillimetres);
    for (auto &w : widths) w *= static_cast<float>(fromMillimetres);
    for (auto &h : heights) h *= static_cast<float>(fromMillimetres);
    printLength *= fromMillimetres;
    printVolume *= fromMillimetres * fromMillimetres * fromMillimetres;
    for (const auto &p : points)
      if (!std::isfinite(p[0]) || !std::isfinite(p[1]) || !std::isfinite(p[2])) return fail("Export coordinate conversion overflow.");
    for (float w : widths) if (!std::isfinite(w)) return fail("Export width conversion overflow.");
    for (float h : heights) if (!std::isfinite(h)) return fail("Export height conversion overflow.");
    if (!std::isfinite(printLength) || !std::isfinite(printVolume)) return fail("Export statistics conversion overflow.");
    if (!stage->GetPrimAtPath(pxr::SdfPath("/AG_block")) ||
        !stage->GetPrimAtPath(pxr::SdfPath("/Slicer_bracing")) ||
        !stage->GetPrimAtPath(pxr::SdfPath("/Slicer_slot")))
      return fail("Input USD must contain /AG_block, /Slicer_bracing and /Slicer_slot.");
    // Replace current output and remove the legacy curve output on re-export.
    for (const auto *name : {"/TDP_gcode", "/Gcode", "/printpath", "/printMesh", "/TDP_printmesh", "/TDP_printpath"}) {
      const pxr::SdfPath outputPath(name);
      if (stage->GetPrimAtPath(outputPath) &&
          (!stage->RemovePrim(outputPath) || stage->GetPrimAtPath(outputPath)))
        return fail(std::string("Cannot replace ") + name + " composed from another USD layer.");
    }
    auto scope = pxr::UsdGeomScope::Define(stage, pxr::SdfPath("/TDP_gcode"));
    auto prim = scope.GetPrim();
    auto metadata = stage->GetRootLayer()->GetCustomLayerData();
    const auto round2 = [](double value) { return std::round(value * 100.0) / 100.0; };
    pxr::VtArray<double> heightDomain, widthDomain;
    heightDomain.push_back(*heightRange.first);
    heightDomain.push_back(*heightRange.second);
    widthDomain.push_back(round2(*widthRange.first));
    widthDomain.push_back(round2(*widthRange.second));
    metadata["blockid"] = pxr::VtValue(blockid);
    metadata["printheightdomain"] = pxr::VtValue(heightDomain);
    metadata["printwidthdomain"] = pxr::VtValue(widthDomain);
    metadata["printlength"] = pxr::VtValue(round2(printLength));
    metadata["printvolume"] = pxr::VtValue(round2(printVolume));
    metadata.erase("printstatisticsvalid");
    metadata["layercount"] = pxr::VtValue(layerCount);
    stage->GetRootLayer()->SetCustomLayerData(metadata);
    // Parallel arrays: element i always describes the same ordered print point.
    const bool ok = prim.CreateAttribute(pxr::TfToken("points"), pxr::SdfValueTypeNames->Point3fArray, true).Set(points) &&
        prim.CreateAttribute(pxr::TfToken("printNormals"), pxr::SdfValueTypeNames->Normal3fArray, true).Set(normals) &&
        prim.CreateAttribute(pxr::TfToken("printWidths"), pxr::SdfValueTypeNames->FloatArray, true).Set(widths) &&
        prim.CreateAttribute(pxr::TfToken("printHeight"), pxr::SdfValueTypeNames->FloatArray, true).Set(heights);
    if (!ok) return fail("Could not author USD Gcode attributes.");
    zPointArray meshPositions;
    zIntArray meshCounts, meshConnects;
    zFnMesh meshFn(printMesh);
    meshFn.getVertexPositions(meshPositions);
    meshFn.getPolygonData(meshConnects, meshCounts);
    pxr::VtArray<pxr::GfVec3f> meshPoints;
    for (const auto &p : meshPositions) {
      const pxr::GfVec3f point(p.x * fromMillimetres, p.y * fromMillimetres, p.z * fromMillimetres);
      if (!std::isfinite(point[0]) || !std::isfinite(point[1]) || !std::isfinite(point[2]))
        return fail("Print mesh coordinate conversion overflow.");
      meshPoints.push_back(point);
    }
    auto mesh = pxr::UsdGeomMesh::Define(stage, pxr::SdfPath("/TDP_printmesh"));
    mesh.GetPrim().SetDisplayName("print mesh");
    if (!mesh.CreatePointsAttr().Set(meshPoints) ||
        !mesh.CreateFaceVertexCountsAttr().Set(pxr::VtArray<int>(meshCounts.begin(), meshCounts.end())) ||
        !mesh.CreateFaceVertexIndicesAttr().Set(pxr::VtArray<int>(meshConnects.begin(), meshConnects.end())) ||
        !mesh.CreateSubdivisionSchemeAttr().Set(pxr::UsdGeomTokens->none))
      return fail("Could not author USD print mesh geometry.");
    // Use the same per-layer curves as the loft. Nonperiodic curves with an
    // explicit repeated endpoint allow closed and open paths in one prim.
    pxr::VtArray<pxr::GfVec3f> curvePoints;
    pxr::VtArray<int> curveCounts;
    pxr::VtArray<float> curveWidths;
    for (const auto &paths : renderPaths) for (const auto &renderPath : paths) {
      zPointArray positions;
      zFnGraph(const_cast<zObjGraph &>(renderPath.graph)).getVertexPositions(positions);
      std::vector<zTsUtils::GraphPolyline> components;
      if (!zTsUtils::graphPolylines(renderPath.graph, components))
        return fail("Invalid renderable print path.");
      for (const auto &component : components) {
        auto ids = component.Vertices;
        if (component.Closed) ids.push_back(ids.front());
        curveCounts.push_back(static_cast<int>(ids.size()));
        for (int id : ids) {
          const auto &p = positions[id];
          const pxr::GfVec3f point(p.x*fromMillimetres, p.y*fromMillimetres, p.z*fromMillimetres);
          const float width = static_cast<float>(renderPath.widths[id]*fromMillimetres);
          if (!std::isfinite(point[0]) || !std::isfinite(point[1]) || !std::isfinite(point[2]) || !std::isfinite(width))
            return fail("Print curve coordinate/width conversion overflow.");
          curvePoints.push_back(point);
          curveWidths.push_back(width);
        }
      }
    }
    auto curvesPrim = pxr::UsdGeomBasisCurves::Define(stage, pxr::SdfPath("/TDP_printpath"));
    curvesPrim.GetPrim().SetDisplayName("print path");
    if (!curvesPrim.CreateTypeAttr().Set(pxr::UsdGeomTokens->linear) ||
        !curvesPrim.CreateWrapAttr().Set(pxr::UsdGeomTokens->nonperiodic) ||
        !curvesPrim.CreateCurveVertexCountsAttr().Set(curveCounts) ||
        !curvesPrim.CreatePointsAttr().Set(curvePoints) ||
        !curvesPrim.CreateWidthsAttr().Set(curveWidths) ||
        !curvesPrim.SetWidthsInterpolation(pxr::UsdGeomTokens->vertex))
      return fail("Could not author USD print path curves.");
    if (extension == ".usda") {
      std::string text;
      if (!stage->GetRootLayer()->ExportToString(&text))
        return fail("Could not serialize USDA output.");
      const auto formatted = formatPrintPathUSDA(text);
      std::ofstream output(filePath, std::ios::binary | std::ios::trunc);
      if (!output) return fail("Could not open USDA output file.");
      output.write(formatted.data(), static_cast<std::streamsize>(formatted.size()));
      output.close();
      if (!output) return fail("Could not write USDA output file.");
    } else if (!stage->GetRootLayer()->Export(filePath)) {
      return fail("Could not write USD output file.");
    }
    return true;
  } catch (const std::exception &e) {
    return fail(e.what());
  }
#endif
}

#if defined(ZSPACE_TS3DP_OPENUSD)
namespace {
nlohmann::json usdSettingsJSON(const pxr::VtValue &value) {
  if (value.IsHolding<pxr::VtDictionary>()) {
    auto result = nlohmann::json::object();
    for (const auto &entry : value.UncheckedGet<pxr::VtDictionary>())
      result[entry.first] = usdSettingsJSON(entry.second);
    return result;
  }
  if (value.IsHolding<bool>()) return value.UncheckedGet<bool>();
  if (value.IsHolding<int>()) return value.UncheckedGet<int>();
  if (value.IsHolding<float>()) {
    const auto v = value.UncheckedGet<float>();
    if (!std::isfinite(v)) throw std::runtime_error("Settings contain a non-finite number.");
    return v;
  }
  if (value.IsHolding<double>()) {
    const auto v = value.UncheckedGet<double>();
    if (!std::isfinite(v)) throw std::runtime_error("Settings contain a non-finite number.");
    return v;
  }
  if (value.IsHolding<std::string>()) return value.UncheckedGet<std::string>();
  auto result = nlohmann::json::array();
  if (value.IsHolding<pxr::VtArray<float>>()) {
    for (auto v : value.UncheckedGet<pxr::VtArray<float>>()) result.push_back(usdSettingsJSON(pxr::VtValue(v)));
    return result;
  }
  if (value.IsHolding<pxr::VtArray<double>>()) {
    for (auto v : value.UncheckedGet<pxr::VtArray<double>>()) result.push_back(usdSettingsJSON(pxr::VtValue(v)));
    return result;
  }
  if (value.IsHolding<pxr::VtArray<int>>()) {
    for (auto v : value.UncheckedGet<pxr::VtArray<int>>()) result.push_back(v);
    return result;
  }
  throw std::runtime_error("Unsupported stage setting type; use JSON-compatible scalar or numeric array values.");
}

pxr::UsdPrim usdGeometry(const pxr::UsdStageRefPtr &stage, const std::string &path,
                         const std::string &kind) {
  const auto root = stage->GetPrimAtPath(pxr::SdfPath(path));
  if (!root) throw std::runtime_error("Missing " + kind + " object: " + path);
  if ((kind == "block" && !root.IsA<pxr::UsdGeomMesh>()) || (kind == "bracing" && !root.IsA<pxr::UsdGeomBasisCurves>()) || (kind == "slotpoint" && !root.IsA<pxr::UsdGeomPoints>()))
    throw std::runtime_error("Incorrect USD geometry type at " + path);
  return root;
}

zFloatArray usdWidthTable(const nlohmann::json &settings, const std::string &key) {
  const auto &table = settings.at("slicer").at(key);
  if (!table.is_array() || table.empty())
    throw std::runtime_error(key + " must be a nonempty width array.");
  zFloatArray result;
  for (const auto &entry : table) {
    if (!entry.is_number()) throw std::runtime_error(key + " contains a nonnumeric width.");
    const float width = entry.get<float>();
    if (!std::isfinite(width) || width <= 0)
      throw std::runtime_error(key + " requires finite positive widths.");
    result.push_back(width);
  }
  return result;
}

pxr::VtArray<int> usdDataIds(const pxr::UsdPrim &prim, const char *name, size_t count, bool optional = false) {
  const auto pv = prim.GetAttribute(pxr::TfToken(name));
  pxr::VtArray<int> values;
  if (!pv && optional) { values.assign(count, 0); return values; }
  if (!pv || !pv.Get(&values) || values.size() != count)
    throw std::runtime_error(prim.GetPath().GetString() + ": " + name + " must be an int[] with one value per corresponding face/curve.");
  return values;
}
zPointArray usdWorldPoints(const pxr::UsdGeomPointBased &geometry, pxr::UsdGeomXformCache &cache, double toMillimetres) {
  pxr::VtArray<pxr::GfVec3f> points;
  if (!geometry.GetPointsAttr().Get(&points)) throw std::runtime_error("Missing points on " + geometry.GetPath().GetString());
  const auto matrix = cache.GetLocalToWorldTransform(geometry.GetPrim());
  zPointArray result;
  result.reserve(points.size());
  for (const auto &point : points) {
    const auto p = matrix.Transform(pxr::GfVec3d(point)) * toMillimetres;
    if (!std::isfinite(p[0]) || !std::isfinite(p[1]) || !std::isfinite(p[2]))
      throw std::runtime_error("Geometry contains non-finite transformed positions.");
    result.push_back(zPoint(p[0], p[1], p[2]));
  }
  return result;
}
}
#endif

ZSPACE_INLINE bool readSlicerInputUSD(const std::string &filePath, SlicerUSDInput &out,
    std::string *outError, const std::string &blockPath,
    const std::string &bracingPath, const std::string &startPath) {
  out = SlicerUSDInput();
  if (outError) outError->clear();
#if !defined(ZSPACE_TS3DP_OPENUSD)
  if (outError) *outError = "USD input requires ZSPACE_TS3DP_OPENUSD and linked OpenUSD libraries.";
  return false;
#else
  try {
    auto stage = pxr::UsdStage::Open(filePath);
    if (!stage) throw std::runtime_error("Could not open USD stage: " + filePath);
    const auto settingsPrim = stage->GetPrimAtPath(pxr::SdfPath("/Slicer_setting"));
    if (!settingsPrim || !settingsPrim.IsA<pxr::UsdGeomScope>())
      throw std::runtime_error("Missing Scope /Slicer_setting.");
    nlohmann::json settings = nlohmann::json::object();
    for (const std::string group : {"slicer", "printSyn"}) {
      settings[group] = nlohmann::json::object();
      const std::string prefix = group + ":";
      for (const auto &attr : settingsPrim.GetAttributes()) {
        const auto name = attr.GetName().GetString();
        if (name.compare(0, prefix.size(), prefix) != 0) continue;
        pxr::VtValue value;
        if (!attr.Get(&value)) throw std::runtime_error("Cannot read setting " + name);
        settings[group][name.substr(prefix.size())] = usdSettingsJSON(value);
      }
      if (settings[group].empty())
        throw std::runtime_error("Slicer_setting requires attributes in namespace " + group);
    }
    out.MetersPerUnit = pxr::UsdGeomGetStageMetersPerUnit(stage);
    out.UpAxis = pxr::UsdGeomGetStageUpAxis(stage).GetString();
    if (!(out.MetersPerUnit > 0) || !std::isfinite(out.MetersPerUnit))
      throw std::runtime_error("Invalid stage metersPerUnit.");
    // Geometry and dimensional settings returned by this reader are always mm.
    // MetersPerUnit retains the original stage unit for callers/export.
    const double toMillimetres = out.MetersPerUnit * 1000.0;
    const auto metadata=stage->GetRootLayer()->GetCustomLayerData();
    const auto heightDomain=metadata.find("printheightdomain");
    if(heightDomain!=metadata.end()){
      if(!heightDomain->second.IsHolding<pxr::VtArray<double>>())
        throw std::runtime_error("customLayerData printheightdomain must be double[2].");
      const auto &domain=heightDomain->second.UncheckedGet<pxr::VtArray<double>>();
      if(domain.size()!=2)throw std::runtime_error("customLayerData printheightdomain must have two values.");
      for(double value:domain){const double mm=value*toMillimetres;
        if(!std::isfinite(mm)||mm<0)throw std::runtime_error("printheightdomain must contain finite nonnegative lengths.");out.PrintHeightDomain.push_back(mm);}
      if(out.PrintHeightDomain[0]>out.PrintHeightDomain[1])throw std::runtime_error("printheightdomain minimum exceeds maximum.");
    }
    const auto scaleSetting = [&](const char *group, const char *key) {
      auto &values = settings[group];
      if (!values.contains(key)) return;
      auto scaleNumber = [&](nlohmann::json &v) {
        if (!v.is_number()) throw std::runtime_error(std::string("Nonnumeric length setting: ") + key);
        const double scaled = v.get<double>() * toMillimetres;
        if (!std::isfinite(scaled)) throw std::runtime_error(std::string("Length setting overflow: ") + key);
        v = scaled;
      };
      auto &value = values[key];
      if (value.is_array()) for (auto &v : value) scaleNumber(v);
      else scaleNumber(value);
    };
    for (const char *key : {"BoundaryWidth", "BracingWidth", "BracingCurveWidths", "OverlapWidth", "FieldCellSize", "LayerSpacing"})
      scaleSetting("slicer", key);
    for (const char *key : {"ResampleLength", "SurfaceChordTolerance", "SeamOffset", "SeamSize"})
      scaleSetting("printSyn", key);
    out.SettingsJSON = settings.dump();
    const auto data = stage->GetPrimAtPath(pxr::SdfPath("/Slicer_data"));
    if (!data || !data.IsA<pxr::UsdGeomScope>())
      throw std::runtime_error("Missing Scope /Slicer_data.");
    const auto block = usdGeometry(stage, blockPath, "block");
    const auto bracing = usdGeometry(stage, bracingPath, "bracing");
    const auto start = usdGeometry(stage, startPath, "slotpoint");
    pxr::UsdGeomXformCache cache(pxr::UsdTimeCode::Default());
    pxr::UsdGeomMesh mesh(block);
    auto positions = usdWorldPoints(mesh, cache, toMillimetres);
    pxr::VtArray<int> counts, indices, tags, holes;
    if (!mesh.GetFaceVertexCountsAttr().Get(&counts) || !mesh.GetFaceVertexIndicesAttr().Get(&indices) || counts.empty())
      throw std::runtime_error("Block has no valid polygon topology.");
    mesh.GetHoleIndicesAttr().Get(&holes);
    if (!holes.empty()) throw std::runtime_error("Block contains USD hole faces; bake a closed polygon mesh first.");
    tags = usdDataIds(data, "block:FaceTag", counts.size());
    size_t total = 0;
    for (int count : counts) {
      if (count < 3) throw std::runtime_error("Block contains a face with fewer than three vertices.");
      total += static_cast<size_t>(count);
    }
    if (total != indices.size()) throw std::runtime_error("Block face counts do not match connectivity.");
    for (int id : indices) if (id < 0 || static_cast<size_t>(id) >= positions.size())
      throw std::runtime_error("Block has an out-of-range vertex index.");
    const auto boundaryWidths = usdWidthTable(settings, "BoundaryWidth");
    const auto bracingWidths = usdWidthTable(settings, "BracingWidth");
    const auto faceWidthIds = usdDataIds(data, "block:widthid", counts.size());
    for (size_t f = 0; f < tags.size(); ++f) {
      if (tags[f] != Bottom && tags[f] != Top && tags[f] != Side)
        throw std::runtime_error("FaceTag must use Bottom=1, Top=2, Side=3.");
      if (tags[f] == Side) {
        if (faceWidthIds[f] < 0 || static_cast<size_t>(faceWidthIds[f]) >= boundaryWidths.size())
          throw std::runtime_error("block face widthid is outside BoundaryWidth.");
      } else if (faceWidthIds[f] != -1)
        throw std::runtime_error("top/bottom face widthid must be -1.");
    }
    out.FaceWidthIds.assign(faceWidthIds.begin(), faceWidthIds.end());
    out.FaceTags.assign(tags.begin(), tags.end());
    zIntArray polygonCounts(counts.begin(), counts.end()), connects(indices.begin(), indices.end());
    pxr::TfToken orientation;
    mesh.GetOrientationAttr().Get(&orientation);
    const bool reverse = (orientation == pxr::UsdGeomTokens->leftHanded) !=
        (cache.GetLocalToWorldTransform(block).GetDeterminant() < 0);
    if (reverse) {
      size_t offset = 0;
      for (int count : polygonCounts) { std::reverse(connects.begin()+offset, connects.begin()+offset+count); offset += count; }
    }
    zFnMesh fn(out.Block);
    fn.create(positions, polygonCounts, connects);
    if (fn.numVertices() != static_cast<int>(positions.size()) || fn.numPolygons() != static_cast<int>(counts.size()))
      throw std::runtime_error("Mesh creation changed topology; face attributes cannot be aligned.");
    pxr::UsdGeomBasisCurves curves(bracing);
    pxr::TfToken curveType, wrap;
    curves.GetTypeAttr().Get(&curveType);
    curves.GetWrapAttr().Get(&wrap);
    if (curveType != pxr::UsdGeomTokens->linear ||
        (wrap != pxr::UsdGeomTokens->nonperiodic && wrap != pxr::UsdGeomTokens->periodic))
      throw std::runtime_error("bracing requires linear periodic or nonperiodic BasisCurves.");
    const auto curvePoints = usdWorldPoints(curves, cache, toMillimetres);
    pxr::VtArray<int> curveCounts;
    if (!curves.GetCurveVertexCountsAttr().Get(&curveCounts) || curveCounts.empty())
      throw std::runtime_error("bracing has no curves.");
    const bool closed = wrap == pxr::UsdGeomTokens->periodic;
    size_t pointCount = 0;
    for (int count : curveCounts) {
      if (count < (closed ? 3 : 2)) throw std::runtime_error("bracing curve has too few points.");
      pointCount += size_t(count);
    }
    if (pointCount != curvePoints.size()) throw std::runtime_error("bracing curve counts do not match points.");
    const auto widthIds = usdDataIds(data, "bracing:widthid", curveCounts.size(), true);
    out.BracingWidthIds.assign(widthIds.begin(), widthIds.end());
    out.BracingCurves.resize(curveCounts.size());
    size_t offset = 0;
    for (size_t c = 0; c < curveCounts.size(); ++c) {
      const int count = curveCounts[c], id = widthIds[c];
      if (id < 0 || size_t(id) >= bracingWidths.size())
        throw std::runtime_error("bracing curve widthid is outside BracingWidth.");
      zPointArray ribPoints(curvePoints.begin() + offset, curvePoints.begin() + offset + count);
      zIntArray ribEdges;
      for (int j = 0; j < (closed ? count : count - 1); ++j) {
        const int next = (j + 1) % count;
        if (zTsUtils::vecDistSq(ribPoints[j], ribPoints[next]) <= 0)
          throw std::runtime_error("bracing contains a zero-length segment.");
        ribEdges.push_back(j); ribEdges.push_back(next);
      }
      zFnGraph(out.BracingCurves[c]).create(ribPoints, ribEdges);
      out.BracingPrintWidths.push_back(bracingWidths[id]);
      offset += count;
    }
    const auto anchors = usdWorldPoints(pxr::UsdGeomPoints(start), cache, toMillimetres);
    if (anchors.size() != 1) throw std::runtime_error("slotpoint must contain exactly one point.");
    out.StartPoint = anchors.front();
    return true;
  } catch (const std::exception &error) {
    out = SlicerUSDInput();
    if (outError) *outError = error.what();
    return false;
  }
#endif
}

} // namespace zTsIO
} // namespace zSpace
