#include <zspace/slicing.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
using namespace zSpace;
namespace {
void require(bool ok,const std::string &s) { if(!ok) throw std::runtime_error(s); }
struct Input {
  zObjMesh mesh; zIntArray tags,widths; zObjGraph bracing;
  zPoint anchor; SlicerSettings settings; PrintSynSetting print;
};
Input synthetic() {
  Input in;
  zPointArray p{{0,0,0},{4,0,0},{4,4,0},{0,4,0},
    {0,0,1},{4.5,0,1},{4.5,4,1},{0,4,1},
    {0,0,2},{5,0,2},{5,4,2},{0,4,2}};
  zIntArray c(10,4),v{0,3,2,1,8,9,10,11,0,1,5,4,1,2,6,5,2,3,7,6,3,0,4,7,4,5,9,8,5,6,10,9,6,7,11,10,7,4,8,11};
  zFnMesh(in.mesh).create(p,c,v); in.tags={1,2,3,4,5,6,3,4,5,6};
  in.widths={-1,-1,0,0,0,0,0,0,0,0}; in.anchor={0,2,0};
  in.settings.LayerSpacing=.5f;in.settings.FieldCellSize=.04f;
  in.settings.BoundaryWidth={.4f};in.settings.DoubleBoundary=true;
  in.print.ResampleLength=.2f; in.print.FeatureAngleTol=91;
  return in;
}
std::vector<double> readArray(std::istream &s) {
  size_t n;require(bool(s>>n) && n<10000000,"Invalid fixture array length");
  std::vector<double> a(n);for(auto &v:a)require(bool(s>>v),"Truncated fixture");return a;
}
Input fixture(const std::string &filename) {
  std::ifstream stream(filename);require(bool(stream),"Cannot open USD numeric fixture");Input in;
  auto a=readArray(stream);zPointArray p;for(size_t i=0;i<a.size();i+=3)p.emplace_back(a[i],a[i+1],a[i+2]);
  a=readArray(stream);zIntArray c(a.begin(),a.end());a=readArray(stream);zIntArray v(a.begin(),a.end());
  zFnMesh(in.mesh).create(p,c,v);a=readArray(stream);in.tags.assign(a.begin(),a.end());
  a=readArray(stream);in.widths.assign(a.begin(),a.end());a=readArray(stream);in.anchor=zPoint(a.at(0),a.at(1),a.at(2));
  a=readArray(stream);p.clear();for(size_t i=0;i<a.size();i+=3)p.emplace_back(a[i],a[i+1],a[i+2]);
  a=readArray(stream);v.assign(a.begin(),a.end());zFnGraph(in.bracing).create(p,v);
  a=readArray(stream);in.settings.BoundaryWidth.assign(a.begin(),a.end());
  a=readArray(stream);in.settings.BracingCurveWidths.assign(a.begin(),a.end());a=readArray(stream);
  in.settings.LayerSpacing=float(a.at(0));in.settings.FieldCellSize=float(a.at(1));in.settings.OverlapWidth=float(a.at(2));
  in.settings.FieldScaleUp=float(a.at(3));in.settings.SlotOffset=float(a.at(4));in.settings.BracingSlot=float(a.at(5));
  in.settings.DoubleBoundary=bool(a.at(6));in.settings.BracingWidth=float(a.at(7));in.print.ResampleLength=float(a.at(8));
  in.print.FeatureAngleTol=float(a.at(9));in.print.SeamOffset=float(a.at(10));in.print.SeamSize=float(a.at(11));
  in.settings._HEqualiseIterations=int(a.at(12));in.settings._HArapIterations=int(a.at(13));return in;
}
void configure(zTsSlicing &s, Input &in,int mode) {
  s.setSlicerMode(static_cast<zTsSlicerMode>(mode));s.setInputMesh(in.mesh,in.tags,in.widths);
  s.setSlotAnchor(in.anchor);if(zFnGraph(in.bracing).numEdges())s.setBracingGraph(in.bracing);
  s.getSettings()=in.settings;s.getSettings().HPrintHeight.Enabled=false;s.getPrintSynSetting()=in.print;
}
std::map<int,int> check(const PrintPathArray &layers, bool requireSplit) {
  require(layers.size()>1,"No layers");std::map<int,int> expected,expectedPointLabels,previousPiece;bool varying=false,split=false;
  for(size_t l=0;l<layers.size();++l) {
    std::map<int,int> counts,pieces,first,pointLabels;
    for(const auto &path:layers[l]) {
      zPointArray p;zIntArray edges;zFnGraph fn(const_cast<zObjGraph &>(path.graph));fn.getVertexPositions(p);fn.getEdgeData(edges);
      require(path.featureFlags.size()==p.size(),"Feature masks lost alignment");
      require(path.faceTags.size()==p.size() && path.widths.size()==p.size(),"Point attributes are misaligned: vertices="+std::to_string(p.size())+" widths="+std::to_string(path.widths.size())+" tags="+std::to_string(path.faceTags.size()));
      require(path.edgeFaceTags.size()==edges.size()/2,"Edge provenance misaligned");
      for(int tag:path.faceTags)if(tag)++pointLabels[tag];
      for(const auto &seg:path.faceTagSegments) {
        require(seg.tag>=3 && seg.tag<=6 && seg.vertices.size()>=2,"Invalid region segment");
        counts[seg.tag]+=int(seg.vertices.size());++pieces[seg.tag];
        if(!first.count(seg.tag)) first[seg.tag]=int(seg.vertices.size());
        for(size_t j=1;j<seg.vertices.size();++j)
          require(zTsUtils::vecDistSq(p.at(seg.vertices[j-1]),p.at(seg.vertices[j]))>1e-18,"Collapsed sampling edge");
      }
    }
    require(counts.size()==4,"A side FaceTag is missing");
    if(l==0)expected=counts;else require(counts==expected,"FaceTag total changed across layers");
    if(l==0)expectedPointLabels=pointLabels;else require(pointLabels==expectedPointLabels,"Per-point FaceTag count changed across layers");
    for(const auto &entry:pieces)if(entry.second==2) {
      split=true;if(previousPiece.count(entry.first) && previousPiece[entry.first]!=first[entry.first]) varying=true;
    }
    previousPiece=first;
    std::cout<<"layer="<<l;for(auto entry:counts)std::cout<<" tag"<<entry.first<<"="<<entry.second<<"/"<<pieces[entry.first]<<"pieces";
    std::cout<<"\n";
  }
  require(!requireSplit || (split && varying),"Staggered slot did not redistribute its fixed budget");
  return expected;
}

void projectedAnchors() {
  zPointArray chartPoints{{0,0,0},{4,0,0},{4,4,0},{0,4,0}};
  zIntArray triangles{0,1,2,0,2,3},counts{3,3};
  zObjMesh chart,section;zFnMesh(chart).create(chartPoints,counts,triangles);
  auto worldPoints=chartPoints;for(auto &p:worldPoints){p.x+=1000;p.y+=1000;}
  zFnMesh(section).create(worldPoints,counts,triangles);
  zObjMeshArray charts{chart,chart},sections{section,section};
  zTsField::FaceTagSource source;source.edgeTags={3,4,5,6};source.edgeWidthIds={0,0,0,0};source.widths={.4f};
  zFnGraph(source.boundary).create(chartPoints,zIntArray{0,1,1,2,2,3,3,0});
  PrintPathArray paths(2);
  for(int l=0;l<2;++l) {
    zPointArray points{{.1f,.1f,0}};
    if(l==0)points.emplace_back(.1000015f,.1f,0); // distinct in the chart, coincident at the 3D offset
    points.emplace_back(3.9f,.1f,0);points.emplace_back(3.9f,3.9f,0);points.emplace_back(.1f,3.9f,0);
    zIntArray edges;for(int i=0;i<int(points.size());++i){edges.push_back(i);edges.push_back((i+1)%int(points.size()));}
    PrintPath path;zFnGraph(path.graph).create(points,edges);path.widths.assign(points.size(),.4f);
    path.featureFlags.assign(points.size(),0);path.featureFlags[0]=PrintPath::AngleFeature;if(!l)path.featureFlags[1]=PrintPath::IntersectionFeature;
    paths[l].push_back(std::move(path));
  }
  PrintSynSetting settings;settings.ResampleLength=.2f;settings.FeatureAngleTol=91;
  std::string status;require(zPrintSyn::resampleFaceTaggedPaths(paths,{source,source},settings,status,&charts,&sections),status);
  std::map<int,int> expected;
  for(int l=0;l<2;++l){std::map<int,int> actual;
    for(const auto &segment:paths[l][0].faceTagSegments)actual[segment.tag]+=int(segment.vertices.size());
    if(!l)expected=actual;else require(actual==expected,"Projected-anchor repair changed the cross-layer budget");
  }
  for(const auto &layer:paths){const auto &p=layer[0];
    zIntArray flagged;for(int i=0;i<int(p.featureFlags.size());++i)if(p.featureFlags[i]){
      require(p.featureFlags[i]==PrintPath::OuterJunctionFeature,"Projected repair retained an angle feature");flagged.push_back(i);}
    require(zPrintSyn::featureVertexIds(p,91)==flagged,"Projected repair junction metadata is stale");
    require(std::none_of(p.featureFlags.begin(),p.featureFlags.end(),[](int mask){return mask&PrintPath::IntersectionFeature;}),
      "Projected repair retained stale intersection identities");}
  PrintPathArray mapped(2);zFloatArray heights(4,.2f);
  for(int l=0;l<2;++l)require(zTsField::mapPrintPaths(paths[l],charts[l],sections[l],&heights,mapped[l],l),"Cannot map repaired anchors");
  zObjMesh mesh;require(zTsUtils::buildPrintMesh(mapped,mesh,&status),status);
}

void run(Input in,int mode,bool parallel,bool external) {
  zTsSlicing s;configure(s,in,mode);
  require(s.computeSlicer(parallel),s.getStatus());
  require(s.computePostProcess(parallel),s.getStatus());
  zObjMesh printMesh;std::string meshError;
  const bool meshOK=zTsUtils::buildPrintMesh(s.get3DPrintPaths(),printMesh,&meshError);
  if(!meshOK) {
    const auto &world=s.get3DPrintPaths(),&flat=s.get2DPrintPaths();int reported=0;
    for(size_t l=0;l<world.size() && reported<4;++l)for(size_t k=0;k<world[l].size() && reported<4;++k) {
      zPointArray xyz,uv;zIntArray edges;zFnGraph w(const_cast<zObjGraph&>(world[l][k].graph));w.getVertexPositions(xyz);w.getEdgeData(edges);
      zFnGraph(const_cast<zObjGraph&>(flat[l][k].graph)).getVertexPositions(uv);
      for(size_t e=0;e<edges.size();e+=2)if(zTsUtils::vecDistSq(xyz[edges[e]],xyz[edges[e+1]])<=1e-18) {
        const int a=edges[e],b=edges[e+1];
        std::cout<<"Coincident projection layer="<<l<<" vertices="<<a<<","<<b<<" tag="<<world[l][k].edgeFaceTags[e/2]
          <<" chartDistance="<<std::sqrt(zTsUtils::vecDistSq(uv[a],uv[b]))<<" world="<<xyz[a].x<<","<<xyz[a].y<<","<<xyz[a].z<<"\n";++reported;
      }
    }
  }
  require(meshOK,meshError);
  auto counts=check(s.get2DPrintPaths(),!external);
  require(check(s.get3DPrintPaths(),false)==counts,"Projection changed FaceTag budgets");
  const auto &combined=s.getResultantPrintPath();
  require(combined.faceTags.size()==size_t(zFnGraph(const_cast<zObjGraph &>(combined.graph)).numVertices()),"Linking lost FaceTags");
  require(combined.edgeFaceTags.size()==size_t(zFnGraph(const_cast<zObjGraph &>(combined.graph)).numEdges()),"Linking lost edge provenance");
  if(mode==0 && !external) {
    for(auto &tag:in.tags)if(isSideFaceTag(tag))tag=Side;
    zTsSlicing legacy;configure(legacy,in,mode);require(legacy.computeSlicer(parallel),legacy.getStatus());
    require(legacy.getPathFields()==s.getPathFields() && legacy.getWidthFields()==s.getWidthFields(),"FaceTag ownership altered geometry/width fields");
  }
  // Every sampled segment must preserve the same indices after 2D->3D mapping.
  const auto &a=s.get2DPrintPaths(),&b=s.get3DPrintPaths();
  for(size_t l=0;l<a.size();++l)for(size_t p=0;p<a[l].size();++p) {
    require(a[l][p].faceTags==b[l][p].faceTags && a[l][p].edgeFaceTags==b[l][p].edgeFaceTags,"Projection relabelled samples");
  }
  std::cout<<"PASS mode="<<mode<<" parallel="<<parallel<<" layers="<<a.size()<<" target2D="<<in.print.ResampleLength<<"\n";
}
}
int main(int argc,char **argv) {
  try {
    projectedAnchors();
    const int mode=argc>1?std::stoi(argv[1]):0;
    auto in=argc>2?fixture(argv[2]):synthetic();
    if(argc>3 && std::string(argv[3])=="slicer-only") {
      zTsSlicing session;configure(session,in,mode);require(session.computeSlicer(false),session.getStatus());
      std::cout<<"PASS slicing-only baseline\n";return 0;
    }
    run(in,mode,false,argc>2);run(in,mode,true,argc>2);
    return 0;
  }catch(const std::exception &e){std::cerr<<"FAIL: "<<e.what()<<"\n";return 1;}
}
