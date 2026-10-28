#include <omp.h>
#include <pxr/base/work/threadLimits.h>
#include <emscripten/emscripten.h>
#include <emscripten/heap.h>
#include <zspace/slicing.h>
#include <array>
#include <vector>
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
using NumericArray = std::vector<double>;
using Arrays = std::array<NumericArray, 12>;
struct Record { int kind=0, layer=-1; std::array<NumericArray,7> arrays; };
struct Session {
  Arrays input;
  std::unique_ptr<zTsSlicing> engine;
  std::vector<Record> records;
  std::vector<std::string> warnings;
  std::string status;
  double seconds=0;
  int mode=0, threads=1, layerCount=0;
  bool parallel=false, processed=false, hasResult=false;
  PrintPath exportPath;
  PrintPathArray exportLayers;
};
std::map<int,Session> sessions;
int nextHandle=1;
std::string lastError;
Session& session(int h) { auto it=sessions.find(h); if(it==sessions.end()) throw std::invalid_argument("Invalid slicing session."); return it->second; }
template<class F> int guard(F work) { try { lastError.clear(); work(); return 1; } catch(const std::exception& e) {lastError=e.what();} catch(...) {lastError="Slicing failed.";} return 0; }
// All boundary numbers are float64 arrays; indices are validated before conversion.
const char* keys[]={"layerSpacing","fieldCellSize","fieldScaleUp","bracingWidth","overlapWidth","doubleBoundary","slotOffset","bracingSlot","traceStepFraction","traceMaxSteps","equaliseIterations","arapIterations","unrollStart","resampleLength","featureAngle","seamOffset","seamSize"};
double number(const NumericArray& options,const char* key,double fallback,double min,double max) {
  size_t i=0; for(;i<17;++i) if(std::string(keys[i])==key) break;
  if(i==17) throw std::invalid_argument("Unknown setting.");
  const double n=options.empty()?fallback:options.at(i);
  if(!std::isfinite(n)||n<min||n>max) throw std::invalid_argument(std::string("Invalid slicing setting: ")+key);
  if((i==5||i==9||i==10||i==11||i==12) && n!=std::floor(n)) throw std::invalid_argument("Expected integral setting.");
  return n;
}
zIntArray integers(const NumericArray& input,size_t max,int minimum=0) {
  if(input.size()>max) throw std::invalid_argument("Invalid slicing index array size.");
  zIntArray out;out.reserve(input.size());
  for(double n:input) { if(!std::isfinite(n)||n!=std::floor(n)||n<minimum||n>INT32_MAX) throw std::invalid_argument("Invalid integer index.");out.push_back(int(n)); }
  return out;
}
zPointArray points(const NumericArray& input) {
  if(input.empty()||input.size()%3||input.size()>150000) throw std::invalid_argument("Invalid slicing coordinates (maximum 50,000 vertices).");
  zPointArray out;out.reserve(input.size()/3);
  for(size_t i=0;i<input.size();i+=3) {
    for(int k=0;k<3;++k) if(!std::isfinite(input[i+k])||std::abs(input[i+k])>1e6) throw std::invalid_argument("Invalid coordinate.");
    out.emplace_back(input[i],input[i+1],input[i+2]);
  } return out;
}
template<class T> NumericArray flat(const T& points) {NumericArray a;a.reserve(points.size()*3);for(const auto& p:points){a.push_back(p.x);a.push_back(p.y);a.push_back(p.z);}return a;}
Record mesh(zObjMesh& object,int kind,int layer=-1) {
  zPointArray p;zIntArray c,i;zFnMesh fn(object);fn.getVertexPositions(p);fn.getPolygonData(i,c);
  Record r;r.kind=kind;r.layer=layer;r.arrays[0]=flat(p);r.arrays[1]=NumericArray(c.begin(),c.end());r.arrays[2]=NumericArray(i.begin(),i.end());return r;
}
Record path(const PrintPath& p,int kind,int layer=-1) {
  zPointArray v;zIntArray e;zFnGraph fn(const_cast<zObjGraph&>(p.graph));fn.getVertexPositions(v);fn.getEdgeData(e);
  Record r;r.kind=kind;r.layer=layer;r.arrays[0]=flat(v);r.arrays[1]=NumericArray(e.begin(),e.end());r.arrays[2]=flat(p.normals);
  r.arrays[3]=NumericArray(p.widths.begin(),p.widths.end());r.arrays[4]=NumericArray(p.heights.begin(),p.heights.end());r.arrays[5]=NumericArray(p.trackIds.begin(),p.trackIds.end());r.arrays[6]=NumericArray(p.faceTags.begin(),p.faceTags.end());return r;
}
void compute(Session& current,int mode,bool useParallel,bool processed) {
  const auto& input=current.input;
  auto positions=points(input[0]);
  auto counts=integers(input[1],50000),connects=integers(input[2],200000);
  if(counts.empty()) throw std::invalid_argument("The slicing mesh has no faces.");
  size_t total=0;
  for(int count:counts){if(count<3||count>4)throw std::invalid_argument("Slicing input requires triangle or quad faces.");total+=count;}
  if(total!=connects.size())throw std::invalid_argument("Connectivity length mismatch.");
  for(int i:connects)if(size_t(i)>=positions.size())throw std::invalid_argument("Vertex index out of range.");
  auto tags=integers(input[3],counts.size());
  if(tags.size()!=counts.size())throw std::invalid_argument("Provide one tag per face.");
  for(int t:tags)if(t<Bottom||t>6)throw std::invalid_argument("Invalid face tag.");
  auto widthIds=integers(input[4],counts.size(),-1);
  if(!widthIds.empty()&&widthIds.size()!=counts.size())throw std::invalid_argument("Width groups must match faces.");
  if(mode<0||mode>1)throw std::invalid_argument("Invalid slicing mode.");
  const auto& options=input[8];
  if(!options.empty()&&options.size()!=17)throw std::invalid_argument("Expected 17 settings.");
  auto candidate=std::make_unique<zTsSlicing>();
  candidate->setSlicerMode(static_cast<zTsSlicerMode>(mode));
  zObjMesh source; zFnMesh(source).create(positions,counts,connects);
  candidate->setInputMesh(source,tags,widthIds);
  const auto& anchor=input[5];
  if(anchor.size()!=3) throw std::invalid_argument("Provide a three-coordinate slot anchor.");
  for(int k=0;k<3;++k) if(!std::isfinite(anchor[k])) throw std::invalid_argument("Invalid slot anchor.");
  candidate->setSlotAnchor(zPoint(anchor[0],anchor[1],anchor[2]));
  if(!input[6].empty()) {
    auto bp=points(input[6]);auto be=integers(input[7],200000);
    if(be.size()%2)throw std::invalid_argument("Invalid bracing edges.");
    for(int v:be)if(size_t(v)>=bp.size())throw std::invalid_argument("Bracing index out of range.");
    zObjGraph bracing;zFnGraph(bracing).create(bp,be);candidate->setBracingGraph(bracing);
  } else if(!input[7].empty()) throw std::invalid_argument("Bracing edges without vertices.");
  auto& settings=candidate->getSettings();
  settings.LayerSpacing=number(options,"layerSpacing",0.5,0.001,1e6);
  settings.FieldCellSize=number(options,"fieldCellSize",0.1,0.001,1e6);
  settings.FieldScaleUp=number(options,"fieldScaleUp",1.05,1.01,2);
  settings.BoundaryWidth.clear();
  const auto widths=input[9].empty()?NumericArray{0.4}:input[9];
  if(widths.size()>256)throw std::invalid_argument("Too many boundary widths.");
  for(double n:widths){if(!std::isfinite(n)||n<=0||n>1e6)throw std::invalid_argument("Invalid boundary width.");settings.BoundaryWidth.push_back(n);}
  settings.BracingWidth=number(options,"bracingWidth",0.4,0.001,1e6);
  if(input[10].size()>50000)throw std::invalid_argument("Too many bracing widths.");
  for(double n:input[10]){if(!std::isfinite(n)||n<=0||n>1e6)throw std::invalid_argument("Invalid bracing width.");settings.BracingCurveWidths.push_back(n);}
  settings.OverlapWidth=number(options,"overlapWidth",0,0,1e6);
  settings.DoubleBoundary=number(options,"doubleBoundary",0,0,1)!=0;
  settings.SlotOffset=number(options,"slotOffset",1,-100,100);
  settings.BracingSlot=number(options,"bracingSlot",0.5,0,1);
  settings._HTraceStepFraction=number(options,"traceStepFraction",0.2,0.001,1);
  settings._HTraceMaxSteps=number(options,"traceMaxSteps",100000,100,1000000);
  settings._HEqualiseIterations=number(options,"equaliseIterations",2,0,20);
  settings._HArapIterations=number(options,"arapIterations",10,0,100);
  settings._HUnrollStart=number(options,"unrollStart",0,0,1);
  auto& post=candidate->getPrintSynSetting();
  post.ResampleLength=number(options,"resampleLength",0.2,0.001,1e6);
  post.FeatureAngleTol=number(options,"featureAngle",90,0,180);
  post.SeamOffset=number(options,"seamOffset",0,-1e6,1e6);
  post.SeamSize=number(options,"seamSize",0,0,1e6);
  zPoint min,max; zFnMesh(source).getBounds(min,max);
  const double diagonal=zTsUtils::vecSub(max,min).length();
  if(settings.FieldCellSize<diagonal/1024 || settings.LayerSpacing<diagonal/1000)
    throw std::invalid_argument("Increase the field cell size or layer spacing for this model.");
  // Use all workers in the hardware-sized pool plus the calling thread.
  // Disable nested teams so they cannot exhaust the pool.
  omp_set_dynamic(0);
  omp_set_num_threads(useParallel ? slicingThreadCapacity() : 1);
  omp_set_max_active_levels(1);
  int threadCount=1;
#pragma omp parallel if(useParallel) reduction(max:threadCount)
  { threadCount=omp_get_num_threads(); }
  const double computeStarted=emscripten_get_now();
  if(!candidate->computeSlicer(useParallel)) throw std::runtime_error(candidate->getStatus());
  const bool success=processed?candidate->computePostProcess(useParallel):candidate->computePrintPath(useParallel);
  if(!success) throw std::runtime_error(candidate->getStatus());
  const double computationSeconds=(emscripten_get_now()-computeStarted)/1000.0;
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
  std::vector<Record> records;const auto& paths=candidate->get3DPrintPaths();
  records.push_back(mesh(*field,0));records.push_back(path(candidate->getResultantPrintPath(),1));
  for(size_t i=0;i<slices->size();++i) {
    records.push_back(mesh((*slices)[i],2,int(i)));records.push_back(mesh((*charts)[i],3,int(i)));
    if(i<paths.size())for(const auto& p:paths[i])records.push_back(path(p,4,int(i)));
    Record r;r.kind=5;r.layer=int(i);r.arrays[0]=NumericArray((*pathFields)[i].begin(),(*pathFields)[i].end());r.arrays[1]=NumericArray((*widthFields)[i].begin(),(*widthFields)[i].end());records.push_back(std::move(r));
  }
  current.status=candidate->getStatus();current.records=std::move(records);current.engine=std::move(candidate);
  current.seconds=computationSeconds;current.mode=mode;current.parallel=useParallel;current.processed=processed;current.threads=threadCount;current.layerCount=int(slices->size());current.hasResult=true;
}
}
extern "C" {
EMSCRIPTEN_KEEPALIVE int zspace_slicing_abi_version(){return 2;}
EMSCRIPTEN_KEEPALIVE int zspace_slicing_capabilities(){return zTsSlicing::supportsSlicerMode(zTsSlicerMode::Harmonic)?3:1;}
EMSCRIPTEN_KEEPALIVE const char* zspace_toolsets_last_error(){return lastError.c_str();}
EMSCRIPTEN_KEEPALIVE int zspace_slicing_create(){int h=0;guard([&]{int id=nextHandle++;sessions.emplace(id,Session{});h=id;});return h;}
EMSCRIPTEN_KEEPALIVE int zspace_slicing_destroy(int h){return guard([&]{session(h);sessions.erase(h);});}
EMSCRIPTEN_KEEPALIVE int zspace_slicing_set_array(int h,int slot,const double* values,size_t count){
  return guard([&]{
    auto& a=session(h).input.at(slot);
    auto address=reinterpret_cast<uintptr_t>(values), heap=emscripten_get_heap_size();
    if(count>10000000 || (count && (!address || address%alignof(double) || address>heap || count>(heap-address)/sizeof(double))))
      throw std::invalid_argument("Invalid slicing array pointer or length.");
    if(count) a.assign(values,values+count); else a.clear();
  });
}
EMSCRIPTEN_KEEPALIVE int zspace_slicing_compute(int h,int mode,int parallel,int processed){
  return guard([&]{if((parallel!=0&&parallel!=1)||(processed!=0&&processed!=1))throw std::invalid_argument("Invalid compute flags.");compute(session(h),mode,parallel!=0,processed!=0);});
}
EMSCRIPTEN_KEEPALIVE int zspace_slicing_record_count(int h){
  int n=-1;guard([&]{auto& s=session(h);if(!s.hasResult)throw std::invalid_argument("Compute slicing first.");n=int(s.records.size());});return n;
}
EMSCRIPTEN_KEEPALIVE int zspace_slicing_record_kind(int h,int r){int n=-1;guard([&]{n=session(h).records.at(r).kind;});return n;}
EMSCRIPTEN_KEEPALIVE int zspace_slicing_record_layer(int h,int r){int n=-2;guard([&]{n=session(h).records.at(r).layer;});return n;}
EMSCRIPTEN_KEEPALIVE const double* zspace_slicing_array_ptr(int h,int r,int slot){const double* p=nullptr;guard([&]{p=session(h).records.at(r).arrays.at(slot).data();});return p;}
EMSCRIPTEN_KEEPALIVE int zspace_slicing_array_count(int h,int r,int slot){int n=-1;guard([&]{n=int(session(h).records.at(r).arrays.at(slot).size());});return n;}
EMSCRIPTEN_KEEPALIVE const double* zspace_slicing_input_ptr(int h,int slot){const double* p=nullptr;guard([&]{p=session(h).input.at(slot).data();});return p;}
EMSCRIPTEN_KEEPALIVE int zspace_slicing_input_count(int h,int slot){int n=-1;guard([&]{n=int(session(h).input.at(slot).size());});return n;}
EMSCRIPTEN_KEEPALIVE double zspace_slicing_result_number(int h,int field){
  double n=std::numeric_limits<double>::quiet_NaN();guard([&]{const auto& s=session(h);if(!s.hasResult)throw std::invalid_argument("Compute slicing first.");
  switch(field){case 0:n=s.seconds;break;case 1:n=s.mode;break;case 2:n=s.threads;break;case 3:n=s.layerCount;break;case 4:n=s.parallel;break;case 5:n=s.processed;break;default:throw std::invalid_argument("Invalid result field.");}});return n;
}
EMSCRIPTEN_KEEPALIVE const char* zspace_slicing_status(int h){const char* p=nullptr;guard([&]{p=session(h).status.c_str();});return p;}
EMSCRIPTEN_KEEPALIVE int zspace_slicing_warning_count(int h){int n=-1;guard([&]{n=int(session(h).warnings.size());});return n;}
EMSCRIPTEN_KEEPALIVE const char* zspace_slicing_warning(int h,int i){const char* p=nullptr;guard([&]{p=session(h).warnings.at(i).c_str();});return p;}
EMSCRIPTEN_KEEPALIVE int zspace_slicing_import_usd(int h,const char* file){
  return guard([&]{
    auto& target=session(h);
    pxr::WorkSetConcurrencyLimit(1);
    zTsIO::SlicerUSDInput input;std::string error;
    if(!file||!zTsIO::importUSD(file,input,&error))throw std::runtime_error(error);
    Arrays a;auto m=mesh(input.Block,0);a[0]=std::move(m.arrays[0]);a[1]=std::move(m.arrays[1]);a[2]=std::move(m.arrays[2]);
    a[3]=NumericArray(input.FaceTags.begin(),input.FaceTags.end());a[4]=NumericArray(input.FaceWidthIds.begin(),input.FaceWidthIds.end());
    a[5]={input.StartPoint.x,input.StartPoint.y,input.StartPoint.z};
    zPointArray bp;zIntArray be;zFnGraph ribs(input.Bracing);ribs.getVertexPositions(bp);ribs.getEdgeData(be);a[6]=flat(bp);a[7]=NumericArray(be.begin(),be.end());
    const auto& s=input.Settings;const auto& p=input.PrintSettings;
    a[8]={s.LayerSpacing,s.FieldCellSize,s.FieldScaleUp,s.BracingWidth,s.OverlapWidth,double(s.DoubleBoundary),s.SlotOffset,s.BracingSlot,
      s._HTraceStepFraction,double(s._HTraceMaxSteps),double(s._HEqualiseIterations),double(s._HArapIterations),double(s._HUnrollStart),
      p.ResampleLength,p.FeatureAngleTol,p.SeamOffset,p.SeamSize};
    a[9]=NumericArray(s.BoundaryWidth.begin(),s.BoundaryWidth.end());a[10]=NumericArray(input.BracingPrintWidths.begin(),input.BracingPrintWidths.end());
    target.input=std::move(a);target.warnings=std::move(input.UnknownSettings);
  });
}
EMSCRIPTEN_KEEPALIVE int zspace_slicing_export_begin(int h,int layers){
  return guard([&]{if(layers<0||layers>10000)throw std::invalid_argument("Invalid layer count.");auto& s=session(h);s.exportPath=PrintPath{};s.exportLayers.clear();s.exportLayers.resize(layers);s.input[5].clear();s.input[6].clear();});
}
EMSCRIPTEN_KEEPALIVE int zspace_slicing_export_add_path(int h,int layer){
  return guard([&]{
    auto& s=session(h);if(layer < -1 || layer>=int(s.exportLayers.size()))throw std::invalid_argument("Invalid export layer.");
    const auto& a=s.input;auto v=points(a[0]);auto e=integers(a[1],v.size()*4);
    if(e.size()%2)throw std::invalid_argument("Expected edge pairs.");
    for(int id:e)if(size_t(id)>=v.size())throw std::invalid_argument("Invalid path edge.");
    if(a[2].size()!=a[0].size()||a[3].size()!=v.size()||a[4].size()!=v.size())throw std::invalid_argument("Path attribute lengths must match vertices.");
    PrintPath p;p.normals=points(a[2]);
    if(!a[6].empty()) {
      p.faceTags=integers(a[6],v.size());
      if(p.faceTags.size()!=v.size() || std::any_of(p.faceTags.begin(),p.faceTags.end(),[](int t){return t>6;}))
        throw std::invalid_argument("Invalid per-point FaceTags.");
    }
    for(double x:a[3]){if(!std::isfinite(x)||x<0||x>1e6)throw std::invalid_argument("Invalid path width.");p.widths.push_back(x);}
    for(double x:a[4]){if(!std::isfinite(x)||x<0||x>1e6)throw std::invalid_argument("Invalid path height.");p.heights.push_back(x);}
    zFnGraph(p.graph).create(v,e);
    if(layer==-1)s.exportPath=std::move(p);else s.exportLayers.at(layer).push_back(std::move(p));
  });
}
EMSCRIPTEN_KEEPALIVE int zspace_slicing_export_usd(int h,const char* source,const char* output,int blockid){
  return guard([&]{auto& s=session(h);pxr::WorkSetConcurrencyLimit(1);std::string error;
    if(!source||!output||!zTsIO::exportUSD(s.exportPath,source,output,int(s.exportLayers.size()),&error,blockid,&s.exportLayers))throw std::runtime_error(error);
  });
}
}
