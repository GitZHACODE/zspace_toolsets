#include <zspace/slicing.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <thread>

using namespace zSpace;
namespace {
void require(bool value, const std::string &message) {
  if (!value) throw std::runtime_error(message);
}
void timingChecks() {
  using namespace zTsUtils;
  setSlicingTimingsEnabled(false);
  clearSlicingTimings();
  { ScopedSlicingTiming disabled("disabled"); }
  require(getSlicingTimings().empty(), "Disabled instrumentation recorded a sample");
  setSlicingTimingsEnabled(true);
  {
    ScopedSlicingTiming stale("stale");
    clearSlicingTimings();
  }
  require(getSlicingTimings().empty(), "A cleared timing leaked into the next run");
  std::vector<std::thread> workers;
  for (int worker = 0; worker < 4; ++worker) workers.emplace_back([] {
    for (int i = 0; i < 100; ++i) { ScopedSlicingTiming timed("parallel", true); }
  });
  for (auto &worker : workers) worker.join();
  try { ScopedSlicingTiming failed("exception"); throw std::runtime_error("expected"); }
  catch (const std::runtime_error &) {}
  setSlicingTimingsEnabled(false);
  const auto rows = getSlicingTimings();
  const auto parallel = std::find_if(rows.begin(), rows.end(), [](const auto &row) { return row.name == "parallel"; });
  const auto failed = std::find_if(rows.begin(), rows.end(), [](const auto &row) { return row.name == "exception"; });
  require(parallel != rows.end() && parallel->calls == 400 && parallel->workerAggregate,
      "Parallel timing lost samples or its aggregate marker");
  require(parallel->totalMilliseconds >= parallel->maxMilliseconds &&
      parallel->maxMilliseconds >= 0 && std::isfinite(parallel->totalMilliseconds),
      "Invalid parallel timing statistics");
  require(failed != rows.end() && failed->calls == 1, "Exception exit lost the method timing");
  clearSlicingTimings();
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
  settings.ResampleLength = spacing;settings.InnerResampleLength=spacing; settings.FeatureAngleTol = angle;
  std::string status;
  std::vector<int> bad;
  require(zPrintSyn::resamplePrintPaths(layers, settings, status, &bad, false), status);
  require(bad.empty(), "Resampling opened a closed curve");
  aligned(layers[0][0]);
  return layers[0][0];
}
zObjMesh chart();
void junctionFeatures(){
  // Two geometric turns bracket a straight side-ownership junction. Only the
  // latter is a fabrication feature, independent of the angle setting.
  auto source=path({{.1f,.1f,0},{.6f,.1f,0},{.85f,.35f,0},{1.1f,.6f,0},{1.6f,.6f,0}});
  source.edgeFaceTags={6,6,3,3};source.trackIds.assign(5,1);
  source.featureFlags={0,PrintPath::AngleFeature,0,PrintPath::IntersectionFeature,0};
  require(zPrintSyn::cornerVertexIds(source,150)==zIntArray({1,3}),"Junction fixture lost its ordinary geometric turns");
  for(float angle:{0.f,90.f,150.f,180.f})require(zPrintSyn::featureVertexIds(source,angle)==zIntArray{2},"Junction selection depends on angle or legacy masks");
  require(zPrintSyn::refreshJunctionFeatures(source)&&source.featureFlags==zIntArray({0,0,4,0,0}),"Dedicated junction flags are incorrect");
  for(auto tags:{zIntArray{},zIntArray{6,0,0,3},zIntArray{3,3,3,3},zIntArray{6,1,2,3}}){
    auto excluded=source;excluded.edgeFaceTags=tags;require(zPrintSyn::featureVertexIds(excluded,150).empty(),"Slot/internal/missing ownership created a feature");}
  auto inner=source;inner.trackIds.assign(5,0);require(zPrintSyn::featureVertexIds(inner,150).empty(),"Inner track became a feature");
  // Vertex/edge storage order is deliberately unrelated to path adjacency.
  const auto points=positions(source);zPointArray permuted{points[3],points[0],points[4],points[2],points[1]};
  auto shuffled=path(permuted);zIntArray shuffledEdges{3,0,1,4,0,2,4,3};
  zFnGraph(shuffled.graph).create(permuted,shuffledEdges);shuffled.edgeFaceTags={3,6,3,6};shuffled.trackIds.assign(5,1);
  require(zPrintSyn::featureVertexIds(shuffled,0)==zIntArray{3},"Feature selection used vertex order instead of graph incidence");
  auto sampled=resample(source,.17f,150);const auto junction=findPoint(sampled,points[2]);
  require(sampled.featureFlags[junction]==PrintPath::OuterJunctionFeature,"Straight junction was dropped by sampling");
  for(int id:{1,3})require(sampled.featureFlags[findPoint(sampled,points[id])]==0,"Geometric support was promoted to a feature");
  auto originalBevel=path({{0,0,0},{1,0,0},{1.8f,0,0},{2,.2f,0},{2,1,0},{2,2,0}});
  originalBevel.edgeFaceTags.assign(5,3);originalBevel.trackIds.assign(6,1);
  auto bevel=resample(originalBevel,.17f,150);
  findPoint(bevel,{1.8f,0,0});findPoint(bevel,{2,.2f,0});
  for(const auto &p:positions(bevel))require(zTsUtils::vecDistSq(p,{2,0,0})>1e-12,"Production sampling inserted a reconstructed intersection");
  originalBevel.edgeFaceTags.assign(5,3);originalBevel.trackIds.assign(6,0);
  auto innerBevel=resample(originalBevel,.17f,150);findPoint(innerBevel,{2,0,0});
  require(zPrintSyn::featureVertexIds(innerBevel,150).empty(),"Inner reconstruction created an outer junction");
  PrintPathArray inheritedInner{{originalBevel}};PrintSynSetting independent;
  independent.ResampleLength=.17f;independent.InnerResampleLength=20;independent.FeatureAngleTol=150;std::string innerStatus;
  require(zPrintSyn::resamplePrintPaths(inheritedInner,independent,innerStatus,nullptr,false),innerStatus);
  findPoint(inheritedInner[0][0],{2,0,0});
  require(positions(inheritedInner[0][0]).size()==3&&std::all_of(inheritedInner[0][0].edgeFaceTags.begin(),inheritedInner[0][0].edgeFaceTags.end(),[](int tag){return tag==0;}),
    "Pure-inner inherited side labels bypassed independent20mm reconstruction/sampling");
  PrintPathArray mapped;PrintPath linked;std::string status;
  require(zPrintSyn::projectPrintPaths({{sampled},{sampled}},{chart(),chart()},{chart(),chart()},nullptr,mapped,status,nullptr,false),status);
  for(const auto &layer:mapped)require(layer[0].featureFlags==sampled.featureFlags,"Projection changed junction identity");
  require(zPrintSyn::linkPrintPaths(mapped,linked,status),status);
  require(zPrintSyn::featureVertexIds(linked,0).size()==2&&
    std::count(linked.featureFlags.begin(),linked.featureFlags.end(),PrintPath::OuterJunctionFeature)==2,
    "Link travel edge created or removed a junction");
  // Both windings exercise label-preserving reversal and seam splitting.
  for(bool clockwise:{false,true}){
    auto ring=path({{.2f,.2f,0},{1,.2f,0},{1.8f,.2f,0},{1.8f,1.8f,0},{.2f,1.8f,0}},true);
    ring.edgeFaceTags={6,3,3,3,6};ring.trackIds.assign(5,1);
    if(clockwise){auto p=positions(ring);std::reverse(p.begin(),p.end());ring=path(p,true);ring.edgeFaceTags={3,3,3,6,6};ring.trackIds.assign(5,1);}
    ring=resample(ring,.17f,150);PrintPathArray layers{{ring}};
    require(zPrintSyn::alignStartPoints(layers,{chart()},.11,0,status,nullptr,false),status);
    const auto &out=layers[0][0];require(zPrintSyn::featureVertexIds(out,0).size()==2,"Reversal/seam changed outer junction count");
    for(size_t v=0;v<out.featureFlags.size();++v)require(out.featureFlags[v]==0||out.featureFlags[v]==PrintPath::OuterJunctionFeature,"Seam exported legacy angle flags");
  }
}
void surfaceContourOutput(){
  zPointArray p;for(int i=0;i<100;++i){const double angle=i*6.283185307179586/100;p.emplace_back(1+.8*std::cos(angle),1+.8*std::sin(angle),0);}
  auto raw=path(p,true);raw.trackIds.assign(p.size(),0);raw.featureFlags.assign(p.size(),0);raw.featureFlags[50]=PrintPath::AngleFeature|PrintPath::IntersectionFeature;
  zObjMeshArray charts{chart()},sections{chart()};PrintSynSetting settings;settings.ResampleLength=.5f;settings.InnerResampleLength=.5f;settings.FeatureAngleTol=0;settings.SeamOffset=0;settings.SeamSize=0;
  PrintPathArray first{{raw}},second=first,a,b,hidden;PrintPath ca,cb;std::string sa,sb;
  require(zPrintSyn::computePipeline(charts,sections,nullptr,zPrintSyn::SeamMapping::WorldNearest,settings,first,a,ca,sa,false),sa);
  zTsUtils::clearSlicingTimings(); zTsUtils::setSlicingTimingsEnabled(true);
  require(zPrintSyn::computePipeline(charts,sections,nullptr,zPrintSyn::SeamMapping::WorldNearest,settings,second,b,cb,sb,false,nullptr,nullptr,&hidden),sb);
  zTsUtils::setSlicingTimingsEnabled(false);
  const auto timingRows=zTsUtils::getSlicingTimings();
  for(const char *name:{"PrintSyn.computePipeline","PrintSyn.resamplePrintPaths",
      "PrintSyn.reconstructCorners","PrintSyn.alignStartPoints","PrintSyn.projectPrintPaths",
      "PrintSyn.linkPrintPaths"})
    require(std::any_of(timingRows.begin(),timingRows.end(),[&](const auto &row){return row.name==name&&row.calls>0;}),
        std::string("Missing postprocessing timing: ")+name);
  require(std::none_of(timingRows.begin(),timingRows.end(),[](const auto &row){
      return row.name=="PrintSyn.resamplePrintPaths3D" || row.name=="PrintSyn.resampleMapped";}),
      "Production pipeline runs a second resampling pass");
  zTsUtils::clearSlicingTimings();
  const auto compare=[](const PrintPath &a,const PrintPath &b){const auto p=positions(a),q=positions(b);require(p.size()==q.size(),"Hidden output changed target count");
    for(size_t i=0;i<p.size();++i)require(zTsUtils::vecDistSq(p[i],q[i])==0,"Hidden output changed target geometry");
    require(a.featureFlags==b.featureFlags&&a.heights==b.heights&&a.widths==b.widths&&a.trackIds==b.trackIds,"Hidden output changed attributes/features");
    require(a.normals.size()==b.normals.size(),"Hidden output changed normal count");for(size_t i=0;i<a.normals.size();++i)require(zTsUtils::vecDistSq(a.normals[i],b.normals[i])==0,"Hidden output changed normals");
    zIntArray ea,eb;zFnGraph(const_cast<zObjGraph&>(a.graph)).getEdgeData(ea);zFnGraph(const_cast<zObjGraph&>(b.graph)).getEdgeData(eb);require(ea==eb,"Hidden output changed connectivity");};
  compare(ca,cb);compare(first[0][0],second[0][0]);compare(a[0][0],b[0][0]);
  // The observer owns the final processed 2D targets; mapping must not
  // redistribute them. Angular preprocessing intentionally replaces raw points.
  PrintPathArray expected{{raw}};std::vector<int> bad;
  require(zPrintSyn::resamplePrintPaths(expected,settings,sa,&bad,false)&&bad.empty(),sa);
  require(zPrintSyn::alignStartPoints(expected,charts,sections,settings.SeamOffset,settings.SeamSize,sa,&bad,false)&&bad.empty(),sa);
  require(hidden.size()==1&&hidden[0].size()==1,"Hidden output lost the processed curve");
  compare(hidden[0][0],expected[0][0]);
  compare(hidden[0][0],second[0][0]);
  // Stretching the surface changes 3D lengths. Final targets must remain the
  // mapped 2D samples, rather than being redistributed by those new lengths.
  zObjMesh stretched=chart();zFnMesh stretchedFn(stretched);
  auto *world=stretchedFn.getRawVertexPositions();
  for(int i=0;i<stretchedFn.numVertices();++i){world[i].x*=3;world[i].y*=2;}
  PrintPathArray flat{{raw}},mapped,observed;PrintPath combined;
  require(zPrintSyn::computePipeline(charts,zObjMeshArray{stretched},nullptr,
      zPrintSyn::SeamMapping::WorldNearest,settings,flat,mapped,combined,sa,false,
      nullptr,nullptr,&observed),sa);
  require(observed.size()==1&&observed[0].size()==1,"Stretched mapping lost sampled targets");
  compare(observed[0][0],flat[0][0]);
  const auto chartPoints=positions(flat[0][0]),worldPoints=positions(mapped[0][0]);
  require(chartPoints.size()==worldPoints.size() &&
      flat[0][0].featureFlags==mapped[0][0].featureFlags,
      "3D mapping changed sampled target count or feature identities");
  for(size_t i=0;i<chartPoints.size();++i)
    require(zTsUtils::vecDistSq(worldPoints[i],zPoint(3*chartPoints[i].x,2*chartPoints[i].y,0))<1e-10,
        "Final target moved away from its sampled 2D surface mapping");
  raw.featureFlags.pop_back();first={{raw}};second=first;
  require(!zPrintSyn::computePipeline(charts,sections,nullptr,zPrintSyn::SeamMapping::WorldNearest,settings,first,a,ca,sa,false),"Old pipeline accepted invalid feature data");
  require(!zPrintSyn::computePipeline(charts,sections,nullptr,zPrintSyn::SeamMapping::WorldNearest,settings,second,b,cb,sb,false,nullptr,nullptr,&hidden),"Observed pipeline accepted invalid feature data");
  require(hidden.empty()&&sa==sb&&positions(ca).empty()&&positions(cb).empty(),"Failed pipeline retained hidden/combined results");
}
void resampling() {
  // Strict threshold semantics: equality is intentionally excluded.
  const auto corner = path({{0,0,0},{1.005,0,0},{1.005,1,0}});
  require(zPrintSyn::cornerVertexIds(corner, 89).empty(), "Below threshold selected a corner");
  require(zPrintSyn::cornerVertexIds(corner, 90).empty(), "Threshold equality changed");
  require(zPrintSyn::cornerVertexIds(corner, 91) == zIntArray{1}, "Above threshold missed a corner");
  auto sampled = resample(corner, 0.2f, 91);
  size_t id = findPoint(sampled, {1.005,0,0});
  require(sampled.featureFlags[id]==0&&zPrintSyn::featureVertexIds(sampled,91).empty(), "Angle support became a junction feature");
  require(std::abs(sampled.widths[id] - corner.widths[1]) < 1e-7 &&
      sampled.heights[id] == corner.heights[1] && sampled.trackIds[id] == corner.trackIds[1],
      "Nearby regular sample replaced a corner's attributes");

  // A near corner pair becomes one target; a tiny whole path cannot retain
  // valid connectivity under a 0.01 mm weld.
  sampled=resample(path({{0,0,0},{1,0,0},{1,.005f,0},{2,1,0}}),.2f,150);
  const auto welded=positions(sampled);zIntArray weldedEdges;zFnGraph(sampled.graph).getEdgeData(weldedEdges);
  for(size_t e=0;e<weldedEdges.size();e+=2)require(zTsUtils::vecDistSq(welded[weldedEdges[e]],welded[weldedEdges[e+1]])>=.0001,
      "Nearby targets did not merge at 0.01 mm");
  PrintPathArray tiny{{path({{0,0,0},{.0004f,0,0},{.0004f,.0003f,0},{0,.0003f,0}},true)}};
  PrintSynSetting tinySettings;tinySettings.ResampleLength=.0001f;tinySettings.InnerResampleLength=.0001f;tinySettings.FeatureAngleTol=91;std::string tinyStatus;
  require(!zPrintSyn::resamplePrintPaths(tiny,tinySettings,tinyStatus,nullptr,false),"A sub-tolerance loop retained invalid targets");

  // Historical metadata does not turn a straight point into an anchor.
  auto coincident = path({{0,0,0},{1,0,0},{1,0,0},{2,0,0}});
  coincident.featureFlags = {0, PrintPath::AngleFeature, PrintPath::IntersectionFeature, 0};
  sampled = resample(coincident, 0.3f, 0);
  require(zPrintSyn::featureVertexIds(sampled,0).empty()&&
      std::none_of(sampled.featureFlags.begin(),sampled.featureFlags.end(),[](int f){return f!=0;}),
      "Historical feature masks survived geometry re-detection");

  // Reconstruction intersects the horizontal/vertical supports at (2,0).
  auto bevel = path({{0,0,0},{1,0,0},{1.8,0,0},{2,0.2,0},{2,1,0},{2,2,0}});
  require(zPrintSyn::reconstructCorners(bevel, 150), "Corner reconstruction failed");
  id = findPoint(bevel, {2,0,0});
  require(bevel.featureFlags[id] & PrintPath::AngleFeature, "Reconstructed corner not re-detected");
  require(zPrintSyn::featureVertexIds(bevel,80).empty(),"Feature API used an old intersection mask");
  // A lower threshold deliberately stops preserving this 90-degree corner.
  sampled = resample(bevel, 0.3f, 80);
  require(std::none_of(sampled.featureFlags.begin(),sampled.featureFlags.end(),[](int f){return f&PrintPath::IntersectionFeature;}),
      "Resampling kept historical intersection identities");
  // Valid endpoints still survive spacing larger than the whole path.
  sampled = resample(path({{0,0,0},{0.02f,0,0}}), 1, 90);
  require(positions(sampled).size() == 2, "Open path lost a separated endpoint");
  sampled=resample(path({{0,0,0},{1,0,0},{1,1,0},{0,1,0}},true),20,90);
  require(positions(sampled).size()>=3&&zTsUtils::isClosedGraph(sampled.graph),
    "Coarse nominal inner spacing collapsed a generic closed path");
  for(bool reverse:{false,true}){
    zPointArray points{{.005f,0,0},{1,0,0},{1,1,0},{0,1,0},{0,0,0}};
    if(reverse){std::reverse(points.begin(),points.end());
      std::rotate(points.begin(),points.end()-1,points.end());}
    auto input=path(points,true);
    size_t cornerId=findPoint(input,{0,0,0});input.heights[cornerId]=17;
    sampled=resample(input,.2f,91);
    require(sampled.heights[findPoint(sampled,{0,0,0})]==17,
      "An arbitrary seam start replaced the nearby geometric corner");
  }

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

  // Both seam mapping strategies preserve current geometric features outside
  // explicit gaps. The threshold must select these corners; metadata alone
  // no longer forces them to remain anchors.
  for (auto strategy : {zPrintSyn::SeamMapping::SharedTopology, zPrintSyn::SeamMapping::WorldNearest}) {
    PrintPathArray paths{{input},{input}};
    PrintSynSetting settings;
    settings.ResampleLength = 0.17f;settings.InnerResampleLength=.17f; settings.FeatureAngleTol = 91;
    settings.SeamOffset = 0.11f;
    require(zPrintSyn::computePipeline(charts, sections, nullptr, strategy, settings,
        paths, mapped, combined, status, false, &heights), status);
    for (int l = 0; l < 2; ++l) {
      aligned(paths[l][0]); aligned(mapped[l][0]);
      for (int i : zPrintSyn::cornerVertexIds(input,91)) {
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
  square.trackIds.assign(4, 1);square.edgeFaceTags={3,4,5,6};
  require(zPrintSyn::refreshJunctionFeatures(square),"Cannot initialize junction metadata");
  PrintPathArray gap{{square}};
  std::vector<int> bad;
  require(zPrintSyn::alignStartPoints(gap, zObjMeshArray{chart()}, 0, 1.8,
      status, &bad, false) && bad.empty(), status);
  aligned(gap[0][0]);
  require(!zTsUtils::isClosedGraph(gap[0][0].graph) &&
      std::count_if(gap[0][0].featureFlags.begin(), gap[0][0].featureFlags.end(),
          [](int flag) { return flag != 0; }) == 2, "Explicit seam cut wrong junction set");
  require(gap[0][0].featureFlags[findPoint(gap[0][0],{.2f,.2f,0})]==0,"Open seam endpoint became a junction");
  for (const zPoint &p : zPointArray{{1.8,0.2,0},{1.8,1.8,0}})
    require(gap[0][0].featureFlags[findPoint(gap[0][0], p)]==PrintPath::OuterJunctionFeature,
        "Seam removed an outer junction outside the cut arc");
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
  for(int f=2;f<10;++f)tags[f]=3+(f-2)%4;
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
    timingChecks(); junctionFeatures(); surfaceContourOutput(); resampling(); mapping(); slicerIntegration(zTsSlicerMode::Topology);
#ifdef ZSPACE_TOOLSETS_WITH_TETGEN
    slicerIntegration(zTsSlicerMode::Harmonic);
#endif
    return 0;
  }
  catch (const std::exception &e) { std::cerr << e.what() << '\n'; return 1; }
}
