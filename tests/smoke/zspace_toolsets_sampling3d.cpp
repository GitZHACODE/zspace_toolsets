#include <zspace/slicing.h>
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace zSpace;
namespace {
void require(bool ok,const char *message){if(!ok)throw std::runtime_error(message);}
double distance(const zPoint &a,const zPoint &b){return (zTsUtils::toEigenVec(a)-zTsUtils::toEigenVec(b)).norm();}
void check(bool closed){
  zObjMesh chart,section;zPointArray flat={{0,0,0},{10,0,0},{10,4,0},{0,4,0}},surface={{0,0,0},{20,0,0},{20,4,4},{0,4,0}};
  zIntArray counts={3,3},faces={0,1,2,0,2,3};zFnMesh(chart).create(flat,counts,faces);zFnMesh(section).create(surface,counts,faces);
  PrintPath path;zPointArray p=closed?zPointArray{{1,1,0},{1.01f,1,0},{9,1,0},{9,3,0},{1,3,0}}:zPointArray{{1,1,0},{2,1,0},{2.01f,1,0},{6,3,0},{9,3,0}};
  zIntArray edges={0,1,1,2,2,3,3,4};if(closed){edges.push_back(4);edges.push_back(0);}zFnGraph(path.graph).create(p,edges);
  path.widths.assign(p.size(),.4f);path.trackIds.assign(p.size(),7);
  path.edgeFaceTags=closed?zIntArray{3,4,5,6,4}:zIntArray{3,4,5,5};
  require(zPrintSyn::refreshJunctionFeatures(path),"Cannot initialize standalone junction anchors");
  PrintPathArray two{{path}},three;zObjMeshArray charts{chart},sections{section};std::string error;
  std::vector<zFloatArray> h{{1,1,3,3}};std::vector<zVectorArray> n{zVectorArray(4,zVector(0,0,1))};
  require(zPrintSyn::projectPrintPaths(two,charts,sections,nullptr,three,error,nullptr,false,&h,&n),error.c_str());
  const auto sourceFeatures=three[0][0].featureFlags;
  require(sourceFeatures==path.featureFlags,"Projection changed source junction masks");
  zPointArray original;zFnGraph(three[0][0].graph).getVertexPositions(original);
  PrintSynSetting settings;settings.ResampleLength=3;
  require(zPrintSyn::resamplePrintPaths3D(two,three,charts,sections,settings,error,false,nullptr,&h,&n),error.c_str());
  zPointArray out,q;zFnGraph(three[0][0].graph).getVertexPositions(out);zFnGraph(two[0][0].graph).getVertexPositions(q);const auto &result=three[0][0];
  require(out.size()==q.size()&&result.widths.size()==out.size()&&result.heights.size()==out.size()&&result.normals.size()==out.size()&&result.trackIds.size()==out.size()&&result.featureFlags.size()==out.size(),"3D sampling attributes lost alignment");
  int ordinary=0;
  for(size_t j=0;j<out.size();++j){require(std::abs(out[j].z-(std::min)(double(out[j].y),.2*out[j].x))<3e-6,"New target left its folded section");
    require(std::abs(out[j].x-2*q[j].x)<3e-6&&std::abs(out[j].y-q[j].y)<3e-6,"Updated 2D/3D targets do not correspond");
    require(std::abs(result.heights[j]-(1+.5*out[j].y))<3e-6,"Existing mesh-height interpolation changed");
    require(result.trackIds[j]==7,"Track identity changed");
    if(result.featureFlags[j]==0){bool endpoint=!closed&&(distance(out[j],original.front())==0||distance(out[j],original.back())==0);if(endpoint)continue;++ordinary;
      for(size_t f=0;f<original.size();++f)if(sourceFeatures[f]||(!closed&&(f==0||f+1==original.size())))require(distance(out[j],original[f])>=1.5-2e-5,"Ordinary 3D sample crowds a mandatory target");}
  }
  for(size_t f=0;f<original.size();++f)if(sourceFeatures[f]||(!closed&&(f==0||f+1==original.size()))){bool found=false;
    for(size_t j=0;j<out.size();++j)if(distance(out[j],original[f])==0&&result.featureFlags[j]==sourceFeatures[f])found=true;
    require(found,"A mandatory short-span/closure feature moved or disappeared");}
  if(!closed)require(ordinary>0,"3D fixture produced no ordinary samples");
  std::vector<zTsUtils::GraphPolyline> curves;require(zTsUtils::graphPolylines(result.graph,curves)&&curves.size()==1&&curves[0].Closed==closed,"Second sampling changed curve topology");
}
}
int main(){try{check(false);check(true);return 0;}catch(const std::exception &error){std::cerr<<error.what()<<'\n';return 1;}}
