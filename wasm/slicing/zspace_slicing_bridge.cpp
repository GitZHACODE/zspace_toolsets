#include <omp.h>
#include <emscripten/emscripten.h>
#include <emscripten/heap.h>
#include <zspace/slicing.h>
#include "toolsets_binary_transfer.h"
#include <algorithm>
#include <cmath>
#include <map>
#include <memory>
#include <limits>
#include <cstdlib>

// The calling slicing Worker is itself one OpenMP thread. The rest of the
// team must fit in the workers created at module initialization.
EM_JS(int, slicingThreadCapacity, (), {
  return Module["slicingThreadCount"];
});

namespace {
using namespace zSpace;
using Json = nlohmann::json;
struct Session { std::unique_ptr<zTsSlicing> engine; std::vector<uint8_t> result; };
std::map<int, Session> sessions;
int nextHandle = 1;
std::string lastError;
Session& session(int handle) {
  const auto found = sessions.find(handle);
  if (found == sessions.end()) throw std::invalid_argument("Invalid slicing session.");
  return found->second;
}
template<class F> int guard(F work) {
  try { lastError.clear(); work(); return 1; }
  catch(const std::exception& error) { lastError = error.what(); }
  catch(...) { lastError = "Slicing failed with an unknown native error."; }
  return 0;
}
double number(const Json& input, const char* key, double fallback, double min, double max) {
  const double n = input.value(key, fallback);
  if (!std::isfinite(n) || n < min || n > max) throw std::invalid_argument(std::string("Invalid slicing setting: ") + key);
  return n;
}
zIntArray integers(const Json& input, size_t max, int minimum = 0) {
  if(!input.is_array() || input.size() > max) throw std::invalid_argument("Invalid slicing index array size.");
  zIntArray result; result.reserve(input.size());
  for(const auto& value : input) {
    if(!value.is_number()) throw std::invalid_argument("Slicing indices must be integers.");
    const double n = value.get<double>();
    if(!std::isfinite(n) || n != std::floor(n) || n < minimum || n > INT32_MAX) throw std::invalid_argument("Slicing indices are outside the allowed integer range.");
    result.push_back(static_cast<int>(n));
  }
  return result;
}
zPointArray points(const Json& input) {
  if(!input.is_array() || input.empty() || input.size()%3 || input.size()>150000) throw std::invalid_argument("Invalid slicing coordinates (maximum 50,000 vertices).");
  zPointArray result;
  for(size_t i=0;i<input.size();i+=3) {
    double xyz[3];
    for(int k=0;k<3;++k) {
      xyz[k]=input[i+k].get<double>();
      if(!std::isfinite(xyz[k]) || std::abs(xyz[k])>1e6) throw std::invalid_argument("Invalid slicing coordinate.");
    }
    result.emplace_back(xyz[0],xyz[1],xyz[2]);
  }
  return result;
}
Json mesh(zObjMesh& object) {
  zFnMesh fn(object); zPointArray p; zIntArray counts, connects;
  fn.getVertexPositions(p); fn.getPolygonData(connects,counts);
  std::vector<double> flat; flat.reserve(p.size()*3);
  for(const auto& v:p) {flat.push_back(v.x);flat.push_back(v.y);flat.push_back(v.z);}
  return {{"vertexPositions",flat},{"faceCounts",counts},{"faceConnects",connects}};
}
Json path(const PrintPath& output) {
  zFnGraph fn(const_cast<zObjGraph&>(output.graph)); zPointArray p; zIntArray edges;
  fn.getVertexPositions(p); fn.getEdgeData(edges);
  std::vector<double> flat,normals;
  for(const auto& v:p) {flat.push_back(v.x);flat.push_back(v.y);flat.push_back(v.z);}
  for(const auto& v:output.normals) {normals.push_back(v.x);normals.push_back(v.y);normals.push_back(v.z);}
  return {{"vertexPositions",flat},{"edgeConnects",edges},{"widths",output.widths},
    {"heights",output.heights},{"normals",normals},{"trackIds",output.trackIds}};
}
void compute(Session& current, const Json& input) {
  const auto& geometry=input.at("mesh");
  auto positions=points(geometry.at("vertexPositions"));
  auto counts=integers(geometry.at("faceCounts"),50000),connects=integers(geometry.at("faceConnects"),200000);
  if(counts.empty()) throw std::invalid_argument("The slicing mesh has no faces.");
  size_t total=0;
  for(int count:counts) { if(count<3 || count>4) throw std::invalid_argument("Slicing input requires triangle or quad faces."); total+=count; }
  if(total!=connects.size()) throw std::invalid_argument("Slicing face counts do not match connectivity.");
  for(int index:connects) if(size_t(index)>=positions.size()) throw std::invalid_argument("Slicing vertex index is out of range.");
  auto tags=integers(input.at("faceTags"),counts.size());
  if(tags.size()!=counts.size()) throw std::invalid_argument("Provide one slicing tag per face.");
  for(int tag:tags) if(tag<Bottom || tag>Side) throw std::invalid_argument("Face tags must be Bottom (1), Top (2), or Side (3).");
  auto widthIds=integers(input.value("faceWidthIds",Json::array()),counts.size(),-1);
  if(!widthIds.empty() && widthIds.size()!=counts.size()) throw std::invalid_argument("Width groups must match the mesh faces.");
  const int mode=static_cast<int>(number(input,"mode",0,0,1));
  if(input.contains("mode") && input.at("mode").get<double>()!=mode) throw std::invalid_argument("Invalid slicing mode.");
  const auto options=input.value("settings",Json::object());
  auto candidate=std::make_unique<zTsSlicing>();
  candidate->setSlicerMode(static_cast<zTsSlicerMode>(mode));
  zObjMesh source; zFnMesh(source).create(positions,counts,connects);
  candidate->setInputMesh(source,tags,widthIds);
  const auto anchor=input.at("anchor");
  if(!anchor.is_array() || anchor.size()!=3) throw std::invalid_argument("Provide a three-coordinate slot anchor.");
  for(int k=0;k<3;++k) if(!anchor[k].is_number() || !std::isfinite(anchor[k].get<double>())) throw std::invalid_argument("Invalid slot anchor.");
  candidate->setSlotAnchor(zPoint(anchor[0].get<double>(),anchor[1].get<double>(),anchor[2].get<double>()));
  if(input.contains("bracing")) {
    const auto& b=input.at("bracing");
    if(!b.at("vertexPositions").empty()) {
      auto bp=points(b.at("vertexPositions")); auto be=integers(b.at("edgeConnects"),200000);
      if(be.size()%2) throw std::invalid_argument("Invalid bracing edge pairs.");
      for(int v:be) if(size_t(v)>=bp.size()) throw std::invalid_argument("Bracing vertex index is out of range.");
      zObjGraph bracing; zFnGraph(bracing).create(bp,be); candidate->setBracingGraph(bracing);
    }
  }
  auto& settings=candidate->getSettings();
  settings.LayerSpacing=number(options,"layerSpacing",0.5,0.001,1e6);
  settings.FieldCellSize=number(options,"fieldCellSize",0.1,0.001,1e6);
  settings.FieldScaleUp=number(options,"fieldScaleUp",1.05,1.01,2);
  settings.BoundaryWidth.clear();
  const auto widths=options.value("boundaryWidths",Json::array({0.4}));
  if(!widths.is_array() || widths.empty() || widths.size()>256) throw std::invalid_argument("Provide between one and 256 boundary widths.");
  for(const auto& width:widths) {
    const double n=width.get<double>(); if(!std::isfinite(n) || n<=0 || n>1e6) throw std::invalid_argument("Boundary widths must be positive and finite.");
    settings.BoundaryWidth.push_back(n);
  }
  settings.BracingWidth=number(options,"bracingWidth",0.4,0.001,1e6);
  if(options.contains("bracingCurveWidths")) {
    const auto& widths=options.at("bracingCurveWidths");
    if(!widths.is_array() || widths.size()>50000) throw std::invalid_argument("Invalid bracing curve widths.");
    for(const auto& width:widths) {
      const double n=width.get<double>(); if(!std::isfinite(n)||n<=0||n>1e6)throw std::invalid_argument("Bracing widths must be positive and finite.");
      settings.BracingCurveWidths.push_back(n);
    }
  }
  settings.OverlapWidth=number(options,"overlapWidth",0,0,1e6);
  settings.DoubleBoundary=options.value("doubleBoundary",false);
  settings.SlotOffset=number(options,"slotOffset",number(options,"boundarySlotOffset",1,-100,100),-100,100);
  // Older clients may provide a shared even/odd value. Different values cannot
  // express the new width-based offset rule; require the new settings explicitly.
  const double even = number(options,"bracingSlotEven",0.5,0,1);
  const double odd = number(options,"bracingSlotOdd",even,0,1);
  if (!options.contains("bracingSlot") && even != odd)
    throw std::invalid_argument("Use bracingSlot and slotOffset instead of separate even/odd slot positions.");
  settings.BracingSlot=number(options,"bracingSlot",even,0,1);
  settings._HTraceStepFraction=number(options,"traceStepFraction",0.2,0.001,1);
  settings._HTraceMaxSteps=number(options,"traceMaxSteps",100000,100,1000000);
  settings._HEqualiseIterations=number(options,"equaliseIterations",2,0,20);
  settings._HArapIterations=number(options,"arapIterations",10,0,100);
  settings._HUnrollStart=number(options,"unrollStart",0,0,1);
  auto& post=candidate->getPrintSynSetting();
  post.ResampleLength=number(options,"resampleLength",0.2,0.001,1e6);
  post.InnerResampleLength=number(options,"innerResampleLength",20.0,0.001,1e6);
  post.FeatureAngleTol=number(options,"featureAngle",90,0,180);
  post.SeamOffset=number(options,"seamOffset",0,-1e6,1e6);
  post.SeamSize=number(options,"seamSize",0,0,1e6);
  zPoint min,max; zFnMesh(source).getBounds(min,max);
  const double diagonal=zTsUtils::vecSub(max,min).length();
  if(settings.FieldCellSize<diagonal/1024 || settings.LayerSpacing<diagonal/1000)
    throw std::invalid_argument("Increase the field cell size or layer spacing for this model.");
  const bool useParallel=input.value("useParallel",false);
  // Use all workers in the hardware-sized pool plus the calling thread.
  // Disable nested teams so they cannot exhaust the pool.
  omp_set_dynamic(0);
  omp_set_num_threads(useParallel ? slicingThreadCapacity() : 1);
  omp_set_max_active_levels(1);
  int threadCount=1;
#pragma omp parallel if(useParallel) reduction(max:threadCount)
  { threadCount=omp_get_num_threads(); }
  if(!candidate->computeSlicer(useParallel)) throw std::runtime_error(candidate->getStatus());
  const bool processed=input.value("postProcess",true);
  const bool success=processed?candidate->computePostProcess(useParallel):candidate->computePrintPath(useParallel);
  if(!success) throw std::runtime_error(candidate->getStatus());
  zObjMeshArray* slices; zObjMeshArray* charts; zObjMeshScalarField* field;
  const std::vector<zScalarArray>* pathFields; const std::vector<zScalarArray>* widthFields;
  if(mode==0) {
    auto* s=candidate->getTopologySlicer(); slices=&s->getIsoMeshes(); charts=&s->getIsoMeshes2D();
    field=&s->getScalarField(); pathFields=&s->getPathFields(); widthFields=&s->getWidthFields();
  } else {
#if defined(ZSPACE_TOOLSETS_WITH_TETGEN)
    auto* s=candidate->getHarmonicSlicer(); slices=&s->getIsoMeshes();
    charts=&const_cast<zObjMeshArray&>(s->getIsoMeshes2D()); field=&const_cast<zObjMeshScalarField&>(s->getScalarField());
    pathFields=&s->getPathFields(); widthFields=&s->getWidthFields();
#else
    throw std::invalid_argument("Harmonic slicing is unavailable in this build.");
#endif
  }
  Json layers=Json::array(); const auto& paths=candidate->get3DPrintPaths();
  for(size_t i=0;i<slices->size();++i) {
    Json layerPaths=Json::array(); if(i>0&&i-1<paths.size()) for(const auto& p:paths[i-1]) layerPaths.push_back(path(p));
    layers.push_back({{"slice",mesh((*slices)[i])},{"chart",mesh((*charts)[i])},
      {"paths",layerPaths},{"pathField",(*pathFields)[i]},{"widthField",(*widthFields)[i]}});
  }
  Json result={{"useParallel",useParallel},{"threadCount",threadCount},{"mode",mode},{"printSectionOffset",1},{"status",candidate->getStatus()},{"layers",layers},
    {"field",mesh(*field)},{"combined",path(candidate->getResultantPrintPath())},{"postProcessed",processed}};
  auto bytes=zspaceTransfer::encode(result);
  current.engine=std::move(candidate); current.result=std::move(bytes);
}
}
extern "C" {
EMSCRIPTEN_KEEPALIVE int zspace_slicing_abi_version(){return 1;}
EMSCRIPTEN_KEEPALIVE int zspace_slicing_capabilities(){return zTsSlicing::supportsSlicerMode(zTsSlicerMode::Harmonic)?3:1;}
EMSCRIPTEN_KEEPALIVE int zspace_toolsets_transfer_version(){return 1;}
EMSCRIPTEN_KEEPALIVE void* zspace_alloc(size_t bytes){return std::malloc(bytes);}
EMSCRIPTEN_KEEPALIVE void zspace_free(void* pointer){std::free(pointer);}
EMSCRIPTEN_KEEPALIVE const char* zspace_toolsets_last_error(){return lastError.c_str();}
EMSCRIPTEN_KEEPALIVE int zspace_slicing_create(){int h=0;guard([&]{h=nextHandle++;sessions.emplace(h,Session{});});return h;}
EMSCRIPTEN_KEEPALIVE int zspace_slicing_destroy(int handle){return guard([&]{session(handle);sessions.erase(handle);});}
EMSCRIPTEN_KEEPALIVE int zspace_slicing_compute_from_buffers(int handle,const void* bytes,size_t size){
  return guard([&]{
    auto address=reinterpret_cast<uintptr_t>(bytes); auto heap=emscripten_get_heap_size();
    if(!address || address>heap || size>heap-address)throw std::invalid_argument("Slicing input buffer is outside WASM memory.");
    compute(session(handle),zspaceTransfer::decode(bytes,size));
  });
}
EMSCRIPTEN_KEEPALIVE const uint8_t* zspace_slicing_result_buffers(int handle){
  const uint8_t* pointer=nullptr;guard([&]{const auto& data=session(handle).result;if(data.empty())throw std::invalid_argument("Compute slicing first.");pointer=data.data();});return pointer;
}
EMSCRIPTEN_KEEPALIVE size_t zspace_toolsets_result_buffer_size(int handle){size_t size=0;guard([&]{size=session(handle).result.size();});return size;}
}
