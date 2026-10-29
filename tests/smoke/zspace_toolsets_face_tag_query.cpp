#include "face_tag_reference.h"
#include <iostream>
#include <limits>
#include <map>
#include <thread>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace zSpace;
namespace {
void require(bool ok,const std::string &message){if(!ok)throw std::runtime_error(message);}
zTsField::FaceTagSource rectangle(bool doubleBoundary,bool varied) {
  zTsField::FaceTagSource source;
  zPointArray points;zIntArray edges;
  const zPoint corners[]={{0,0,0},{4,0,0},{4,4,0},{0,4,0}};
  for(int side=0;side<4;++side)for(int j=0;j<24;++j){const auto &a=corners[side],&b=corners[(side+1)%4];
    const float t=float(j)/24;points.emplace_back(a.x+t*(b.x-a.x),a.y+t*(b.y-a.y),0);
    source.edgeTags.push_back(3+side);source.edgeWidthIds.push_back(side);}
  for(int i=0;i<int(points.size());++i){edges.push_back(i);edges.push_back((i+1)%int(points.size()));}
  zFnGraph(source.boundary).create(points,edges);
  source.widths=varied?std::vector<float>{.4f,.8f,.4f,.6f}:std::vector<float>{.4f,.4f,.4f,.4f};
  source.doubleBoundary=doubleBoundary;source.overlap=.05f;
  return source;
}
zPointArray queryPoints(){
  zPointArray points;
  for(int y=0;y<=82;++y)for(int x=0;x<=82;++x)points.emplace_back(-1+6.f*x/82,-1+6.f*y/82,0);
  for(float delta:{-2e-7f,-1e-8f,0.f,1e-8f,2e-7f}){
    points.emplace_back(2+delta,2,0);points.emplace_back(.2f+delta,.2f,0);
    points.emplace_back(3.8f+delta,2,0);points.emplace_back(4+delta,4,0);}
  return points;
}
void comparePaths(const PrintPath &a,const PrintPath &b){
  zPointArray pa,pb;zIntArray ea,eb;
  zFnGraph(const_cast<zObjGraph&>(a.graph)).getVertexPositions(pa);zFnGraph(const_cast<zObjGraph&>(b.graph)).getVertexPositions(pb);
  zFnGraph(const_cast<zObjGraph&>(a.graph)).getEdgeData(ea);zFnGraph(const_cast<zObjGraph&>(b.graph)).getEdgeData(eb);
  require(pa.size()==pb.size()&&ea==eb,"Transition splitting changed graph topology");
  for(size_t i=0;i<pa.size();++i)require(zTsUtils::vecDistSq(pa[i],pb[i])==0,"Transition anchor moved");
  require(a.faceTags==b.faceTags&&a.edgeFaceTags==b.edgeFaceTags&&a.featureFlags==b.featureFlags,
    "Transition classification or feature identities changed");
  require(a.widths==b.widths&&a.heights==b.heights&&a.trackIds==b.trackIds&&a.normals.size()==b.normals.size(),
    "Transition point attributes changed");
  for(size_t i=0;i<a.normals.size();++i)require(zTsUtils::vecDistSq(a.normals[i],b.normals[i])==0,"Transition normals changed");
  require(a.faceTagSegments.size()==b.faceTagSegments.size(),"FaceTag piece count changed");
  for(size_t i=0;i<a.faceTagSegments.size();++i)
    require(a.faceTagSegments[i].tag==b.faceTagSegments[i].tag&&a.faceTagSegments[i].vertices==b.faceTagSegments[i].vertices,
      "FaceTag piece budgets or vertex identities changed");
}
void parity(const zTsField::FaceTagSource &source){
  const auto points=queryPoints();zIntArray expected,actual,wrapper;
  require(face_tag_reference::queryFaceTags(source,points,expected),"Reference classifier failed");
  const zTsField::FaceTagQuery query(source);require(query.valid(),"Prepared classifier failed");
  require(query.query(points,actual)&&zTsField::queryFaceTags(source,points,wrapper),"Bulk classification failed");
  require(actual==expected&&wrapper==expected,"Bulk FaceTags changed from the scalar baseline");
  for(size_t i=0;i<points.size();++i)require(query.label(points[i])==expected[i],"Scalar/prepared FaceTags differ");
  zPointArray pathPoints{{.15f,.15f,0},{3.85f,.15f,0},{3.85f,3.85f,0},{.15f,3.85f,0}};
  zIntArray edges{0,1,1,2,2,3,3,0};PrintPath original;zFnGraph(original.graph).create(pathPoints,edges);
  original.featureFlags={PrintPath::AngleFeature,0,PrintPath::IntersectionFeature,0};
  original.widths={.4f,.5f,.6f,.7f};original.heights={.2f,.3f,.4f,.5f};
  original.trackIds={0,1,0,1};original.normals.assign(4,zVector(0,0,1));
  for(bool split:{false,true}){
    auto before=original,after=original;
    require(face_tag_reference::labelFaceTagSegments(before,source,split),"Reference transition split failed");
    require(zTsField::labelFaceTagSegments(after,query,split),"Prepared transition split failed");comparePaths(before,after);}
  std::vector<zIntArray> parallel(4);std::vector<std::thread> workers;
  for(int i=0;i<4;++i)workers.emplace_back([&,i]{query.query(points,parallel[i]);});
  for(auto &worker:workers)worker.join();
  for(const auto &tags:parallel)require(tags==expected,"Concurrent prepared query changed labels");
  // Direct splitter rejects every malformed optional attribute transactionally.
  for(int attribute=0;attribute<4;++attribute){auto bad=original;
    if(attribute==0)bad.heights.pop_back();else if(attribute==1)bad.normals.pop_back();
    else if(attribute==2)bad.trackIds.pop_back();else bad.featureFlags.pop_back();
    const auto unchanged=bad;require(!zTsField::labelFaceTagSegments(bad,query),"Malformed splitter attributes accepted");
    comparePaths(bad,unchanged);
    PrintPathArray layers{{original},{bad}},before=layers;PrintSynSetting settings;settings.ResampleLength=.2f;settings.FeatureAngleTol=91;std::string status;
    require(!zPrintSyn::resampleFaceTaggedPaths(layers,{source,source},settings,status,nullptr,nullptr,true),"Malformed parallel sampler attributes accepted");
    for(size_t l=0;l<layers.size();++l)comparePaths(layers[l][0],before[l][0]);
  }
#ifdef _OPENMP
  // A layer-parallel caller must use the serial chunk path, with private scratch.
  std::vector<zIntArray> nested(2);
#pragma omp parallel for num_threads(2)
  for(int i=0;i<2;++i)query.query(points,nested[i]);
  for(const auto &tags:nested)require(tags==expected,"Nested caller changed bulk labels");
#endif
}
void parallelSampling(){
  const auto source=rectangle(false,false);PrintPath path;
  zPointArray points{{.15f,.15f,0},{3.85f,.15f,0},{3.85f,3.85f,0},{.15f,3.85f,0}};
  zIntArray edges{0,1,1,2,2,3,3,0};zFnGraph(path.graph).create(points,edges);
  path.widths.assign(4,.4f);path.heights.assign(4,.2f);path.normals.assign(4,zVector(0,0,1));
  path.trackIds.assign(4,1);path.featureFlags.assign(4,PrintPath::AngleFeature);
  PrintPathArray serial(6,std::vector<PrintPath>{path}),parallel=serial;
  const std::vector<zTsField::FaceTagSource> sources(6,source);
  PrintSynSetting settings;settings.ResampleLength=.2f;settings.FeatureAngleTol=91;
  std::string sa,sb;
  require(zPrintSyn::resampleFaceTaggedPaths(serial,sources,settings,sa),sa);
  require(zPrintSyn::resampleFaceTaggedPaths(parallel,sources,settings,sb,nullptr,nullptr,true),sb);
  require(sa==sb&&serial.size()==parallel.size(),"Parallel preparation changed sampling budgets");
  for(size_t l=0;l<serial.size();++l){require(serial[l].size()==parallel[l].size(),"Parallel preparation changed path count");
    for(size_t p=0;p<serial[l].size();++p)comparePaths(serial[l][p],parallel[l][p]);}
  PrintPathArray failed(6,std::vector<PrintPath>{path});failed[2][0].widths.clear();const auto original=failed;
  require(!zPrintSyn::resampleFaceTaggedPaths(failed,sources,settings,sb,nullptr,nullptr,true),
    "A failed parallel layer was not reported");
  for(size_t l=0;l<failed.size();++l)comparePaths(failed[l][0],original[l][0]);
}
PrintPath anchorPath(zPointArray points,bool closed=true){
  PrintPath path;zIntArray edges;
  for(int i=1;i<int(points.size());++i){edges.push_back(i-1);edges.push_back(i);}
  if(closed){edges.push_back(int(points.size())-1);edges.push_back(0);}
  zFnGraph(path.graph).create(points,edges);path.widths.assign(points.size(),.4f);
  path.heights.assign(points.size(),.2f);path.normals.assign(points.size(),zVector(0,0,1));path.trackIds.assign(points.size(),1);
  // Deliberately misleading metadata must not create protected targets.
  path.featureFlags.assign(points.size(),PrintPath::IntersectionFeature);
  return path;
}
void anchorSampling(){
  const auto source=rectangle(false,false);
  zPointArray chartPoints{{0,0,0},{4,0,0},{4,4,0},{0,4,0}};
  zIntArray triangleCounts{3,3},triangles{0,1,2,0,2,3};
  zObjMesh chart;zFnMesh(chart).create(chartPoints,triangleCounts,triangles);
  PrintSynSetting settings;settings.ResampleLength=.2f;settings.FeatureAngleTol=91;
  const auto totals=[](const PrintPath &path,bool labels){std::map<int,int> result;
    if(labels){for(int tag:path.faceTags)if(tag)++result[tag];}
    else for(const auto &s:path.faceTagSegments)result[s.tag]+=int(s.vertices.size());
    return result;};
  // Near-category anchors span the closed seam in one winding and the first
  // edge in the other. The original corner and its attributes must win.
  for(bool reverse:{false,true})for(float scale:{.1f,1.f,10.f}){
    const float offset=scale==.1f?.04f:.004f;
    zPointArray points{{.15f+offset,.15f,0},{3.85f,.15f,0},{3.85f,3.85f,0},{.15f,3.85f,0}};
    if(reverse)std::reverse(points.begin(),points.end());
    const zPoint corner(.15f+offset,.15f,0);auto path=anchorPath(points);
    for(size_t i=0;i<points.size();++i)if(zTsUtils::vecDistSq(points[i],corner)==0)path.heights[i]=1.25f;
    auto second=path;
    // Include dense ordinary input samples on a straight segment. They must
    // not form a transitive weld cluster with a distant protected corner.
    zPointArray dense=points;auto a=dense[0],b=dense[1];
    dense.insert(dense.begin()+1,{a+(b-a)*.001f,a+(b-a)*.002f,a+(b-a)*.003f});
    second=anchorPath(dense);
    for(size_t i=0;i<dense.size();++i)if(zTsUtils::vecDistSq(dense[i],corner)==0)second.heights[i]=1.25f;
    PrintPathArray serial{{path},{second}},parallel=serial;
    auto worldPoints=chartPoints;for(auto &p:worldPoints){p.x*=scale;p.y*=scale;}
    zObjMesh section;zFnMesh(section).create(worldPoints,triangleCounts,triangles);
    const zObjMeshArray charts{chart,chart},sections{section,section};std::string sa,sb;
    require(zPrintSyn::resampleFaceTaggedPaths(serial,{source,source},settings,sa,&charts,&sections,false),sa);
    require(zPrintSyn::resampleFaceTaggedPaths(parallel,{source,source},settings,sb,&charts,&sections,true),sb);
    require(sa==sb,"Anchor cleanup changed parallel budgets");
    for(size_t l=0;l<serial.size();++l){const auto &out=serial[l][0];comparePaths(out,parallel[l][0]);
      require(totals(out,false)==totals(serial[0][0],false)&&totals(out,true)==totals(serial[0][0],true),
        "Anchor cleanup changed cross-layer categorical counts");
      zPointArray p;zIntArray edges;zFnGraph fn(const_cast<zObjGraph&>(out.graph));fn.getVertexPositions(p);fn.getEdgeData(edges);
      bool exact=false;for(size_t i=0;i<p.size();++i)if(zTsUtils::vecDistSq(p[i],corner)==0){exact=true;
        require(out.heights[i]==1.25f,"Near transition replaced exact corner attributes");}
      require(exact&&zTsUtils::isClosedGraph(out.graph),"Near transition removed the corner or opened the seam");
      require(zPrintSyn::featureVertexIds(out,91)==zPrintSyn::cornerVertexIds(out,91),"Features were not re-detected");
      for(size_t e=0;e<edges.size();e+=2)require(zTsUtils::vecDistSq(p[edges[e]],p[edges[e+1]])>=.0001,
        "Chart targets violated 0.01 mm clearance");
      for(int f:out.featureFlags)require(!(f&PrintPath::IntersectionFeature),"Sticky feature mask survived cleanup");
    }
  }
  // Spatially near endpoints with a long intervening arc must remain separate.
  settings.FeatureAngleTol=0;
  auto detour=anchorPath({{1,.1f,0},{2,.1f,0},{2,.104f,0},{1,.104f,0}},false);
  PrintPathArray paths{{detour}};std::string status;
  require(zPrintSyn::resampleFaceTaggedPaths(paths,{source},settings,status),status);
  require(zFnGraph(paths[0][0].graph).numVertices()>8&&zPrintSyn::featureVertexIds(paths[0][0],0).empty(),
    "A long narrow detour collapsed or historical metadata became features");
  // A sub-tolerance categorical piece cannot remain a ghost run. If it is the
  // only occurrence of that tag on one layer, the equal-count request fails.
  settings.FeatureAngleTol=91;
  PrintPathArray missing{{anchorPath({{1,.15f,0},{.154f,.15f,0},{.154f,.158f,0}},false)},
                        {anchorPath({{1,.15f,0},{.154f,.15f,0},{.154f,.3f,0}},false)}};
  const auto original=missing;
  require(!zPrintSyn::resampleFaceTaggedPaths(missing,{source,source},settings,status)&&status.find("missing on layer")!=std::string::npos,
    "Removing the only short categorical run did not report its missing tag");
  for(size_t l=0;l<missing.size();++l)comparePaths(missing[l][0],original[l][0]);
  // Tag 3 still exists after its short second piece disappears. Shared-endpoint
  // ownership then makes equal segment AND unique point-label counts impossible.
  const auto splitSource=[](float end){zTsField::FaceTagSource s;
    zPointArray p{{0,0,0},{1,0,0},{2,0,0},{end,0,0},{4,0,0},{4,4,0},{0,4,0}};
    zIntArray edges;for(int i=0;i<int(p.size());++i){edges.push_back(i);edges.push_back((i+1)%int(p.size()));}
    zFnGraph(s.boundary).create(p,edges);s.edgeTags={3,4,3,5,5,5,5};s.edgeWidthIds.assign(7,0);s.widths={.4f};return s;};
  auto line=anchorPath({{.15f,.1f,0},{3.85f,.1f,0}},false);
  PrintPathArray changed{{line},{line}};const auto unchanged=changed;
  require(!zPrintSyn::resampleFaceTaggedPaths(changed,{splitSource(2.02f),splitSource(2.005f)},settings,status)&&
      status.find("equal segment and point-label counts are incompatible")!=std::string::npos,
      "Short-run removal silently broke unique point-label totals: "+status);
  for(size_t l=0;l<changed.size();++l)comparePaths(changed[l][0],unchanged[l][0]);
  settings.FeatureAngleTol=std::numeric_limits<float>::quiet_NaN();changed={{line}};
  require(!zPrintSyn::resampleFaceTaggedPaths(changed,{source},settings,status),"Nonfinite tagged feature threshold accepted");
}
void tieOrder(){
  zTsField::FaceTagSource source;source.widths={.4f};zPointArray points;zIntArray edges;
  for(int i=0;i<24;++i){const float y=i<3?float(1e-7+7e-9*i):float(100+i);
    points.emplace_back(-1,y,0);points.emplace_back(1,y,0);edges.push_back(2*i);edges.push_back(2*i+1);
    source.edgeTags.push_back(i<3?6-i:3);source.edgeWidthIds.push_back(0);}
  zFnGraph(source.boundary).create(points,edges);zPointArray probes{{0,0,0}};zIntArray expected,actual;
  require(face_tag_reference::queryFaceTags(source,probes,expected)&&expected[0]==4,"Tie-chain fixture lost its order dependence");
  require(zTsField::queryFaceTags(source,probes,actual)&&actual==expected,"Tolerance tie chain was reordered");
  std::reverse(edges.begin(),edges.end());std::reverse(source.edgeTags.begin(),source.edgeTags.end());
  zFnGraph(source.boundary).create(points,edges);
  require(face_tag_reference::queryFaceTags(source,probes,expected)&&expected[0]==6,"Reversed tie-chain fixture was ineffective");
  require(zTsField::queryFaceTags(source,probes,actual)&&actual==expected,"Reversed source order changed categorical ties");
}
void snapshots(){
  auto source=rectangle(false,false);const auto points=queryPoints();zIntArray initial,retained,updated;
  const zTsField::FaceTagQuery old(source);require(old.query(points,initial),"Snapshot query failed");
  for(auto &tag:source.edgeTags)tag=tag==3?6:tag==6?3:tag;
  require(old.query(points,retained)&&retained==initial,"Changing source mutated a prepared snapshot");
  require(zTsField::FaceTagQuery(source).query(points,updated)&&updated!=initial,"Rebuilding a snapshot retained stale ownership");
  source.edgeTags.pop_back();zIntArray invalid{42};
  require(!zTsField::FaceTagQuery(source).valid()&&!zTsField::queryFaceTags(source,points,invalid)&&invalid.empty(),
    "Invalid source retained stale result labels");
}
}
int main(){try{
  parity(rectangle(false,false));parity(rectangle(false,true));
  auto source=rectangle(true,true);source.slotWidth=.4f;
  zPointArray slot{{3.7f,1.8f,0},{4.3f,1.8f,0}};zIntArray pair{0,1};zFnGraph(source.slot).create(slot,pair);
  zPointArray rib{{.2f,.5f,0},{3.8f,3.5f,0}};zObjGraph graph;zFnGraph(graph).create(rib,pair);
  source.bracing.emplace_back(graph,1,.35f,0,.4f);parity(source);tieOrder();snapshots();parallelSampling();anchorSampling();
  std::cout<<"PASS prepared/batched FaceTag parity, transitions, ties, snapshots and concurrent queries\n";return 0;
}catch(const std::exception &error){std::cerr<<"FAIL "<<error.what()<<"\n";return 1;}}
