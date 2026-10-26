#include <zspace/slicing.h>
#include <zspace/toolsets.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace zSpace;

namespace {
void require(bool condition, const std::string& message) {
  if (!condition) throw std::runtime_error(message);
}

zObjectMesh block() {
  zObjectMesh mesh;
  zPointArray points = {
    {0,0,0},{4,0,0},{4,4,0},{0,4,0},
    {0,0,1},{4,0,1},{4,4,1},{0,4,1},
    {0,0,2},{4,0,2},{4,4,2},{0,4,2}
  };
  zIntArray counts(10, 4);
  zIntArray vertices = {
    0,3,2,1, 8,9,10,11,
    0,1,5,4, 1,2,6,5, 2,3,7,6, 3,0,4,7,
    4,5,9,8, 5,6,10,9, 6,7,11,10, 7,4,8,11
  };
  zFnMesh(mesh).create(points, counts, vertices);
  return mesh;
}

void configure(zTsSlicing& session, zObjectMesh& mesh) {
  zIntArray tags(10, Side);
  tags[0] = Bottom;
  tags[1] = Top;
  session.setInputMesh(mesh, tags);
  session.setSlotAnchor(zPoint(0, 0, 0));
  auto& settings = session.getSettings();
  settings.LayerSpacing = 0.5f;
  settings.BoundaryWidth = {0.4f};
  settings.FieldCellSize = 0.1f;
  session.getPrintSynSetting().ResampleLength = 0.2f;
}

void checkPaths(const PrintPathArray& layers) {
  require(layers.size() == 4, "Expected four printable layers; exact top cap must be omitted");
  for (const auto& layer : layers) {
    require(!layer.empty(), "Empty print-path layer");
    for (const auto& path : layer) {
      auto& graph = const_cast<zObjGraph&>(path.graph);
      zFnGraph fn(graph);
      require(fn.numVertices() > 2 && fn.numEdges() > 1, "Degenerate path");
      require(path.widths.size() == static_cast<size_t>(fn.numVertices()), "Missing per-vertex print widths");
      require(path.heights.size() == static_cast<size_t>(fn.numVertices()), "Missing per-vertex print heights");
      for (float height : path.heights)
        require(std::isfinite(height) && height > 0, "Non-positive or non-finite print height");
      for (float width : path.widths)
        require(std::isfinite(width) && width > 0, "Non-positive or non-finite print width");
      const auto* points = fn.getRawVertexPositions();
      for (int i = 0; i < fn.numVertices(); ++i) {
        const auto& p = points[i];
        require(std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z), "Non-finite path position");
        require(p.x >= -0.01 && p.x <= 4.01 && p.y >= -0.01 && p.y <= 4.01 && p.z >= -0.01 && p.z <= 2.01,
                "Path escaped the input block");
      }
    }
  }
}

void topology() {
  auto mesh = block();
  zTsSlicing session;
  configure(session, mesh);
  require(session.computeSlicer(false), session.getStatus());
  require(session.getIsoMeshes().size() == 5, "Unexpected slice count");
  require(session.getPathFields().size() == 5, "Missing layer fields");
  require(session.computePrintPath(false), session.getStatus());
  checkPaths(session.get3DPrintPaths());
  require(session.computePostProcess(false), session.getStatus());
  checkPaths(session.get3DPrintPaths());
  auto& combined = const_cast<PrintPath&>(session.getResultantPrintPath());
  require(zFnGraph(combined.graph).numEdges() > 0, "PrintSyn did not produce a combined path");

  require(session.computeSlicer(false), session.getStatus());
  require(session.computePostProcess(false), session.getStatus());
  zTsSlicing moved(std::move(session));
  require(moved.computeSlicer(false), moved.getStatus());
  require(moved.computePostProcess(false), moved.getStatus());
  checkPaths(moved.get3DPrintPaths());
  configure(session, mesh);
  require(session.computeSlicer(false), "Moved-from session could not be reused");
}

void invalidInput() {
  auto mesh = block();
  zTsSlicing session;
  configure(session, mesh);
  session.clearSlotAnchor();
  require(!session.computeSlicer(false), "Missing anchor was accepted");
  require(session.getStatus().find("Slot") != std::string::npos, "Missing anchor diagnostic");
  session.setSlotAnchor(zPoint(0,0,0));
  session.setInputMesh(mesh, zIntArray(1, Bottom));
  require(!session.computeSlicer(false), "Invalid face tags were accepted");
  require(session.get3DPrintPaths().empty(), "Failure retained stale print paths");
  zTsIO::SlicerUSDInput usdInput;
  usdInput.Block = mesh;
  std::string ioError;
  require(!zTsIO::readSlicerInputUSD("zspace-missing-slicer-input.usda", usdInput, &ioError), "Missing USD input was accepted");
  require(!ioError.empty() && zFnMesh(usdInput.Block).numVertices() == 0, "USD failure did not clear data and report a diagnostic");
  require(!zTsSlicing::supportsSlicerMode(static_cast<zTsSlicerMode>(999)), "Unknown mode reported supported");
  bool rejected = false;
  try { session.setSlicerMode(static_cast<zTsSlicerMode>(999)); }
  catch (const std::invalid_argument&) { rejected = true; }
  require(rejected && session.getSlicerMode() == zTsSlicerMode::Topology, "Mode validation changed the session");
#ifndef ZSPACE_TOOLSETS_WITH_TETGEN
  rejected = false;
  try { session.setSlicerMode(zTsSlicerMode::Harmonic); }
  catch (const std::invalid_argument&) { rejected = true; }
  require(rejected, "Disabled harmonic mode was accepted");
#endif
}

void features() {
  auto mesh = block();
  zTsSlicing session;
  configure(session, mesh);
  zIntArray tags(10, Side), widthIds(10, 0);
  tags[0] = Bottom; tags[1] = Top;
  widthIds[2] = widthIds[6] = 1;
  session.setInputMesh(mesh, tags, widthIds);
  session.getSettings().BoundaryWidth = {0.4f, 0.6f};
  session.getSettings().DoubleBoundary = true;
  session.getSettings().FieldCellSize = 0.05f;
  zObjGraph bracing;
  zPointArray points = {{0,2,0},{4,2,0}};
  zIntArray edges = {0,1};
  zFnGraph(bracing).create(points, edges);
  session.setBracingGraph(bracing);
  require(session.getWidthGroupCount() == 2, "Width groups were lost");
  require(session.computeSlicer(false), session.getStatus());
  require(session.getTopologySlicer()->getBracingGraphs2D().size() == 5, "Bracing was not mapped to each layer");
  const auto& widths = session.getWidthFields();
  bool varied = false;
  for (const auto& layer : widths) {
    if (layer.empty()) continue;
    for (double width : layer) if (std::abs(width - layer.front()) > 1e-3) varied = true;
  }
  require(varied, "Per-group print widths were lost in field computation");
  require(session.computePrintPath(false), session.getStatus());
  checkPaths(session.get3DPrintPaths());
}

void parallelParity() {
  auto mesh = block();
  zTsSlicing serial, parallel;
  configure(serial, mesh);
  configure(parallel, mesh);
  require(serial.computeSlicer(false), serial.getStatus());
  require(parallel.computeSlicer(true), parallel.getStatus());
  require(serial.getPathFields() == parallel.getPathFields(), "Serial and parallel path fields differ");
  require(serial.getWidthFields() == parallel.getWidthFields(), "Serial and parallel width fields differ");
  require(serial.computePostProcess(false), serial.getStatus());
  require(parallel.computePostProcess(true), parallel.getStatus());
  const auto& a = serial.get3DPrintPaths();
  const auto& b = parallel.get3DPrintPaths();
  require(a.size() == b.size(), "Serial and parallel layer counts differ");
  for (size_t layer = 0; layer < a.size(); ++layer) {
    require(a[layer].size() == b[layer].size(), "Serial and parallel path counts differ");
    for (size_t i = 0; i < a[layer].size(); ++i) {
      zPointArray pa, pb;
      zFnGraph(const_cast<zObjGraph&>(a[layer][i].graph)).getVertexPositions(pa);
      zFnGraph(const_cast<zObjGraph&>(b[layer][i].graph)).getVertexPositions(pb);
      require(pa.size() == pb.size(), "Serial and parallel vertex counts differ");
      require(a[layer][i].widths == b[layer][i].widths, "Serial and parallel widths differ");
      require(a[layer][i].heights == b[layer][i].heights, "Serial and parallel heights differ");
      for (size_t j = 0; j < pa.size(); ++j)
        require(zTsUtils::vecDistSq(pa[j], pb[j]) < 1e-12, "Serial and parallel positions differ");
    }
  }
}

void tetRefinement() {
#ifdef ZSPACE_TOOLSETS_WITH_TETGEN
  auto mesh = block();
  zPointArray original;
  zIntArray faces, counts;
  zFnMesh(mesh).getVertexPositions(original);
  zFnMesh(mesh).getPolygonData(faces, counts);
  std::vector<std::array<int,3>> triangles;
  int cursor = 0;
  for (int count : counts) {
    for (int i = 1; i + 1 < count; ++i)
      triangles.push_back({faces[cursor], faces[cursor+i], faces[cursor+i+1]});
    cursor += count;
  }
  const double sourceEdge = zTsUtils::meanEdgeLength(original, triangles);
  size_t coarseCount = 0;
  require(SlicerSettings{}._HTetEdgeLengthScale == 1.0, "Unexpected default tet refinement");
  for (double scale : {1.0, 0.75, 0.5}) {
    zTsSlicing session;
    configure(session, mesh);
    session.setSlicerMode(zTsSlicerMode::Harmonic);
    session.getSettings()._HTetEdgeLengthScale = scale;
    require(session.computeSlicer(false), session.getStatus());
    const auto &tet = session.getHarmonicSlicer()->getTetVolume();
    const double targetEdge = sourceEdge * scale;
    const double limit = targetEdge*targetEdge*targetEdge / (6*std::sqrt(2.0));
    double maximum = 0;
    for (const auto &element : tet.Volumes) if (element.Valid) {
      require(std::isfinite(element.Volume) && element.Volume > 0, "Invalid refined tetrahedron");
      maximum = (std::max)(maximum, element.Volume);
    }
    if (scale == 1) coarseCount = tet.Volumes.size();
    if (scale == 0.5) require(tet.Volumes.size() > coarseCount, "Refinement did not increase tet density");
    std::cout << "Tet edge scale=" << scale << ", count=" << tet.Volumes.size()
        << ", maximum volume=" << maximum << ", target volume=" << limit << '\n';
    require(maximum <= limit*(1+1e-5), "TetGen exceeded the requested volume bound");
    for (const auto &vertex : tet.Vertices) {
      const auto &p = vertex.Position;
      require(p.x >= -1e-6 && p.x <= 4+1e-6 && p.y >= -1e-6 && p.y <= 4+1e-6 &&
          p.z >= -1e-6 && p.z <= 2+1e-6, "Tet refinement changed the source domain");
    }
    require(session.getHarmonicSlicer()->getIsoMeshes().size() == 5 && session.getSettings().LayerSpacing == 0.5f &&
        session.getSettings().FieldCellSize == 0.1f && session.getPrintSynSetting().ResampleLength == 0.2f,
        "Tet refinement changed unrelated slicing settings");
  }
  zPointArray after;
  zIntArray afterFaces, afterCounts;
  zFnMesh(mesh).getVertexPositions(after);
  zFnMesh(mesh).getPolygonData(afterFaces, afterCounts);
  require(after.size() == original.size() && afterFaces == faces && afterCounts == counts,
      "Input surface topology was modified");
  for (size_t i = 0; i < after.size(); ++i)
    require(zTsUtils::vecDistSq(after[i], original[i]) == 0, "Input surface vertices moved");

  zTsSlicing invalid;
  configure(invalid, mesh);
  invalid.setSlicerMode(zTsSlicerMode::Harmonic);
  require(invalid.computeSlicer(false), invalid.getStatus());
  // A failed UV trace does not invalidate a geometrically valid cap face.
  auto *slicer = invalid.getHarmonicSlicer();
  auto &tet = const_cast<HarmonicField &>(slicer->getTetVolume());
  int topFaces = 0, failedVertex = -1;
  for (size_t f = 0; f < tet.Faces.size(); ++f)
    if (tet.FaceTags[f] == HarmonicField::TetFaceTag::Top) {
      ++topFaces; failedVertex = tet.Faces[f].VertexIds[0];
    }
  require(topFaces > 0 && failedVertex >= 0, "No cap for UV-failure regression");
  tet.UVWValid[failedVertex] = false;
  tet.UVW[failedVertex].x = tet.UVW[failedVertex].y = std::numeric_limits<float>::quiet_NaN();
  require(slicer->computeContours(false), "Failed UV trace removed valid cap geometry");
  require(zFnMesh(slicer->getIsoMeshes().back()).numPolygons() == topFaces,
      "Cap triangles were discarded solely because UV tracing failed");
  require(slicer->computeField(false), "Failed UV trace prevented planar cap unfolding");
  require(invalid.computePostProcess(false), invalid.getStatus());
  for (double scale : {0.0, -1.0, 1.1, std::numeric_limits<double>::quiet_NaN(),
       std::numeric_limits<double>::infinity(), 1e-200}) {
    invalid.getSettings()._HTetEdgeLengthScale = scale;
    require(!invalid.computeSlicer(false), "Invalid tet edge-size scale accepted");
    require(invalid.getHarmonicSlicer()->getTetVolume().Volumes.empty() && invalid.get3DPrintPaths().empty(),
        "Invalid tet refinement retained stale results");
  }
#endif
}

void harmonic() {
#ifdef ZSPACE_TOOLSETS_WITH_TETGEN
  auto mesh = block();
  zTsSlicing session;
  configure(session, mesh);
  session.setSlicerMode(zTsSlicerMode::Harmonic);
  require(session.getHarmonicSlicer() != nullptr, "Harmonic operator unavailable");
  require(session.computeSlicer(false), session.getStatus());
  require(session.computePrintPath(false), session.getStatus());
  checkPaths(session.get3DPrintPaths());
  const auto checkUniform = [](const PrintPathArray &paths, double expected) {
    for (const auto &layer : paths) for (const auto &path : layer) {
      require(path.heights.size() == size_t(zFnGraph(const_cast<zObjGraph &>(path.graph)).numVertices()),
              "Height/point alignment lost");
      for (float h : path.heights)
        require(std::abs(h - expected) < 2e-4, "Incorrect analytic harmonic spacing");
    }
  };
  checkUniform(session.get3DPrintPaths(), 0.5);
  require(std::abs(session.getPrintHeightDomain().min - 0.5) < 2e-4 &&
          std::abs(session.getPrintHeightDomain().max - 0.5) < 2e-4,
          "Harmonic height domain is not the measured spacing");
  require(session.getPrintVolume() > 0, "Harmonic volume is zero despite positive heights");
  auto *slicer = session.getHarmonicSlicer();
  const auto &heights = slicer->getVertexLayerHeights();
  require(heights.size() == 5, "Missing per-layer harmonic height tables");
  for (size_t l = 0; l < heights.size(); ++l)
    require(heights[l].size() == size_t(zFnMesh(const_cast<zObjMesh &>(slicer->getIsoMeshes()[l])).numVertices()),
            "Height table does not match its independent iso mesh");
  require(session.computePostProcess(false), session.getStatus());
  checkPaths(session.get3DPrintPaths());
  checkUniform(session.get3DPrintPaths(), 0.5);
  checkUniform(PrintPathArray{{session.getResultantPrintPath()}}, 0.5);
  require(session.getPrintVolume() > 0, "PrintSyn discarded harmonic height/volume");
  const auto serialCombinedHeights = session.getResultantPrintPath().heights;

  // A distinct constant on each independent section makes an accidental shared
  // height table observable, including after resampling and seam alignment.
  std::vector<zFloatArray> perLayer = heights;
  for (size_t l = 0; l < perLayer.size(); ++l)
    perLayer[l].assign(perLayer[l].size(), 0.2f + 0.1f * float(l));
  require(session.computePrintPath(false), session.getStatus());
  auto paths2D = session.get2DPrintPaths();
  PrintPathArray projected;
  PrintPath combined;
  std::string status;
  require(zPrintSyn::computePipeline(slicer->getIsoMeshes2D(), slicer->getIsoMeshes(), nullptr,
      zPrintSyn::SeamMapping::WorldNearest, session.getPrintSynSetting(), paths2D, projected,
      combined, status, false, &perLayer), status);
  for (size_t l = 0; l < projected.size(); ++l)
    checkUniform(PrintPathArray{projected[l]}, 0.2 + 0.1 * l);
  require(!combined.heights.empty(), "Linked path lost per-layer heights");
  for (float h : combined.heights)
    require(h >= 0.2f - 1e-5f && h <= 0.6f + 1e-5f, "Linked height escaped the section ranges");
  const auto rejectTables = [&]() {
    require(!zPrintSyn::projectPrintPaths(paths2D, slicer->getIsoMeshes2D(), slicer->getIsoMeshes(),
        nullptr, projected, status, nullptr, false, &perLayer), "Invalid height tables accepted");
    require(projected.empty(), "Invalid projection retained stale paths");
  };
  perLayer[0].pop_back();
  rejectTables();
  perLayer = heights;
  perLayer[1][0] = std::numeric_limits<float>::quiet_NaN();
  rejectTables();
  perLayer[1][0] = -1.0f;
  rejectTables();
  perLayer.pop_back();
  rejectTables();

  zTsSlicing parallel;
  configure(parallel, mesh);
  parallel.setSlicerMode(zTsSlicerMode::Harmonic);
  require(parallel.computeSlicer(true), parallel.getStatus());
  require(parallel.computePrintPath(true), parallel.getStatus());
  checkUniform(parallel.get3DPrintPaths(), 0.5);
  require(heights == parallel.getHarmonicSlicer()->getVertexLayerHeights(), "Parallel harmonic height tables differ");
  require(parallel.computePostProcess(true), parallel.getStatus());
  checkUniform(parallel.get3DPrintPaths(), 0.5);
  require(serialCombinedHeights == parallel.getResultantPrintPath().heights,
          "Serial/parallel linked harmonic heights differ");

  // Rigid rotation and uniform scaling have independent analytic expectations:
  // local normal height doubles, despite a changed world-Z separation.
  auto transformed = block();
  zFnMesh transformedFn(transformed);
  auto *points = transformedFn.getRawVertexPositions();
  const double c = std::cos(0.5), s = std::sin(0.5);
  for (int i = 0; i < transformedFn.numVertices(); ++i) {
    const zPoint p = points[i];
    points[i] = zPoint(2 * (c*p.x + s*p.z), 2*p.y, 2 * (-s*p.x + c*p.z));
  }
  zTsSlicing rotated;
  configure(rotated, transformed);
  rotated.getSettings().LayerSpacing = 1.1f; // ceil(4 / 1.1) intervals give actual height 1.0.
  rotated.getSettings().BoundaryWidth = {0.8f};
  rotated.getSettings().FieldCellSize = 0.2f;
  rotated.setSlicerMode(zTsSlicerMode::Harmonic);
  require(rotated.computeSlicer(false), rotated.getStatus());
  size_t rotationLayer=0;
  for (const auto &layer : rotated.getHarmonicSlicer()->getVertexLayerHeights()) {
    require(!layer.empty(), "Missing rotated section heights");
    for (size_t i=0;i<layer.size();++i) {float h=layer[i];
      require(std::abs(h - 1.0) < 2e-4, "Incorrect rotated harmonic height on layer " + std::to_string(rotationLayer) +
          ", vertex " + std::to_string(i) + ": " + std::to_string(h));}
    ++rotationLayer;
  }

  auto tapered = block();
  zFnMesh taperedFn(tapered);
  points = taperedFn.getRawVertexPositions();
  for (int i = 0; i < taperedFn.numVertices(); ++i) points[i].x *= 1.0 + 0.15 * points[i].z;
  zTsSlicing varying;
  configure(varying, tapered);
  varying.setSlicerMode(zTsSlicerMode::Harmonic);
  require(varying.computeSlicer(false), varying.getStatus());
  require(varying.computePrintPath(false), varying.getStatus());
  auto range = varying.getPrintHeightDomain();
  require(std::isfinite(range.min) && range.min > 0 && range.max - range.min > 1e-4,
          "Nonuniform harmonic spacing was replaced by a constant height");
  require(varying.computePostProcess(false), varying.getStatus());
  require(varying.getPrintHeightDomain().min > 0 && varying.getPrintVolume() > 0,
          "Varying harmonic heights were lost during post-processing");

  // An unusable gradient must fail and clear the previous successful output.
  auto &tet = const_cast<HarmonicField &>(slicer->getTetVolume());
  for (size_t t = 0; t < tet.Volumes.size(); ++t) if (tet.Volumes[t].Valid) {
    tet.VolumeFieldGradients[t].setZero();
    break;
  }
  require(!slicer->computeContours(false), "Zero harmonic gradient accepted for height calculation");
  require(slicer->getVertexLayerHeights().empty() && slicer->getPrintPaths3D().empty() &&
          slicer->getPrintHeightDomain().max == 0 && slicer->getPrintVolume() == 0,
          "Failed height computation retained stale output/statistics");
  require(session.computeSlicer(false), session.getStatus());
  require(session.computePrintPath(false), session.getStatus());
  checkUniform(session.get3DPrintPaths(), 0.5);
  const auto settings = slicer->getSettings();
  slicer->setSettings(settings);
  require(slicer->getVertexLayerHeights().empty() && slicer->getPrintPaths3D().empty() &&
          slicer->getPrintHeightDomain().max == 0 && slicer->getPrintVolume() == 0,
          "Settings invalidation retained stale harmonic heights/statistics");
  require(!slicer->computePrintPath(false), "Paths computed without height tables/contours");

  // The full harmonic entry point validates processed targets, including
  // exact inclusive bounds, rather than merely clamping their attributes.
  for(const auto &bounds:std::vector<std::pair<double,double>>{{.45,.55},{.45,.5},{double(std::nextafter(.5f,0.f)),.55}}){
    zTsSlicing bounded;configure(bounded,mesh);bounded.setSlicerMode(zTsSlicerMode::Harmonic);
    auto &target=bounded.getSettings().HPrintHeight;target.Enabled=true;target.Minimum=bounds.first;target.Maximum=bounds.second;
    require(bounded.computeSlicer(false),bounded.getStatus());require(bounded.computePostProcess(false),bounded.getStatus());
    require(bounded.getHarmonicSlicer()->getPrintHeightOptimizationReport().Converged,"Successful harmonic output is not compliant");
    for(float h:bounded.getResultantPrintPath().heights)require(h>=bounds.first&&h<=bounds.second,"Final harmonic target escaped inclusive bounds");
  }
  zTsSlicing impossible;configure(impossible,mesh);impossible.setSlicerMode(zTsSlicerMode::Harmonic);
  auto &target=impossible.getSettings().HPrintHeight;target.Enabled=true;target.Minimum=.8;target.Maximum=1;target.MaxDisplacement=.001;target.AllowLayerCountChange=false;
  require(impossible.computeSlicer(false),impossible.getStatus());
  require(!impossible.computePostProcess(false),"Infeasible harmonic domain accepted");
  require(zFnGraph(const_cast<PrintPath&>(impossible.getResultantPrintPath()).graph).numVertices()==0&&!impossible.hasPrintStatistics(),"Infeasible output retained a fabrication path/statistics");
  zTsSlicing adaptive;configure(adaptive,mesh);adaptive.setSlicerMode(zTsSlicerMode::Harmonic);
  auto &adaptiveTarget=adaptive.getSettings().HPrintHeight;adaptiveTarget=target;adaptiveTarget.AllowLayerCountChange=true;
  require(adaptive.computeSlicer(false)&&adaptive.computePostProcess(false),adaptive.getStatus());
  require(adaptive.getHarmonicSlicer()->getIsoMeshes().size()<impossible.getHarmonicSlicer()->getIsoMeshes().size(),"Adaptive count did not remove an unnecessary layer");
  for(float h:adaptive.getResultantPrintPath().heights)require(h>=.8&&h<=1,"Count adaptation claimed a false domain");
  zTsSlicing scheduled;configure(scheduled,mesh);scheduled.setSlicerMode(zTsSlicerMode::Harmonic);
  auto &schedule=scheduled.getSettings().HPrintHeight;schedule=target;schedule.OptimizeSurfaces=false;
  require(scheduled.computeSlicer(false)&&scheduled.computePostProcess(false),scheduled.getStatus());
  require(scheduled.getHarmonicSlicer()->getPrintHeightOptimizationReport().AcceptedIterations==0,"Scheduling-only mode moved surfaces for the target domain");
  target.Maximum=target.Minimum;require(!impossible.computeSlicer(false),"Degenerate harmonic domain accepted");
  zTsSlicing unaffected;configure(unaffected,mesh);unaffected.getSettings().HPrintHeight.Enabled=true;unaffected.getSettings().HPrintHeight.Minimum=-1;
  require(unaffected.computeSlicer(false)&&unaffected.computePostProcess(false),"Harmonic-only target changed topology slicing");

  auto roof=block();zFnMesh roofFn(roof);roofFn.getRawVertexPositions()[10].z=2.7f;
  zTsSlicing coarse,capFixed;configure(coarse,roof);configure(capFixed,roof);coarse.setSlicerMode(zTsSlicerMode::Harmonic);capFixed.setSlicerMode(zTsSlicerMode::Harmonic);
  capFixed.getSettings().HPrintHeight.Enabled=true;capFixed.getSettings().HPrintHeight.Minimum=.2;capFixed.getSettings().HPrintHeight.Maximum=.7;
  require(coarse.computeBasicStructure(false)&&capFixed.computeBasicStructure(false),"Roof cap refinement failed");
  const auto &before=coarse.getHarmonicSlicer()->getTetVolume(),&after=capFixed.getHarmonicSlicer()->getTetVolume();
  const auto capOnly=[](const HarmonicField &t){int count=0;for(const auto &cell:t.Volumes)if(cell.Valid){const int tag=t.vertextag[cell.VertexIds[0]];bool same=tag==Bottom||tag==Top;for(int id:cell.VertexIds)same&=t.vertextag[id]==tag;count+=same;}return count;};
  const int repaired=capOnly(before);require(capOnly(after)==0,"Finite-volume fixed-cap pocket survived refinement");
  require(after.Vertices.size()==before.Vertices.size()+repaired&&after.Volumes.size()==before.Volumes.size()+4*repaired,"Cap refinement was not a local internal-face split");
  require(before.Surface.Triangles==after.Surface.Triangles&&before.Surface.Positions.size()==after.Surface.Positions.size(),"Cap repair changed exterior connectivity");
  for(size_t i=0;i<before.Surface.Positions.size();++i){const auto &a=before.Surface.Positions[i],&b=after.Surface.Positions[i];require(a.x==b.x&&a.y==b.y&&a.z==b.z,"Cap repair moved the input surface");}
  double volumeBefore=0,volumeAfter=0;for(const auto &t:before.Volumes)if(t.Valid)volumeBefore+=t.Volume;for(const auto &t:after.Volumes)if(t.Valid)volumeAfter+=t.Volume;
  require(std::abs(volumeAfter-volumeBefore)<1e-8*volumeBefore,"Cap repair changed enclosed volume");
  require(capFixed.getHarmonicSlicer()->computeHarmonicField(),"Cap-repaired field solve failed");
  for(size_t i=before.Vertices.size();i<after.Vertices.size();++i)require(after.vertextag[i]==Interior&&after.UVW[i].z>0&&after.UVW[i].z<1,"New cap-interface node did not gain a free field value");
#endif
}
}

int main(int argc, char** argv) {
  try {
    const std::string test = argc > 1 ? argv[1] : "topology";
    if (test == "topology") topology();
    else if (test == "invalid") invalidInput();
    else if (test == "features") features();
    else if (test == "parallel") parallelParity();
    else if (test == "harmonic") harmonic();
    else if (test == "tet_refinement") tetRefinement();
    else throw std::runtime_error("Unknown test case");
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
