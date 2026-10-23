#include <zspace/slicing.h>
#include <zspace/zToolsets/fabrication/slicing/utils/IsoRelaxation.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <stdexcept>

using namespace zSpace;
using namespace zSpace::zTsRelax;
namespace {
void require(bool value,const std::string &message) { if(!value)throw std::runtime_error(message); }
double norm(const zPoint &a,const zPoint &b) {return (zTsUtils::toEigenVec(a)-zTsUtils::toEigenVec(b)).norm();}
Reference box() {
  Reference r;
  r.Positions={{0,0,0},{4,0,0},{4,4,0},{0,4,0},{0,0,2},{4,0,2},{4,4,2},{0,4,2}};
  r.Triangles={{0,2,1},{0,3,2},{4,5,6},{4,6,7},
               {0,1,5},{0,5,4},{1,2,6},{1,6,5},{2,3,7},{2,7,6},{3,0,4},{3,4,7}};
  r.WidthIds={-1,-1,-1,-1,0,0,0,0,0,0,1,1};return r;
}
std::vector<Layer> stack(bool bent) {
  std::vector<Layer> layers(5);
  for(int k=0;k<5;++k){auto &m=layers[k];
    for(int y=0;y<=4;++y)for(int x=0;x<=4;++x){double z=.5*k;
      if(bent&&k>0&&k<4)z+=(k==2?.30:.15)*(1-std::abs(x-2)/2.0)+.04*std::sin(y*3.141592653589793/4);
      m.Positions.emplace_back(x,y,z);m.UV.emplace_back(x/4.0,y/4.0,z/2);m.Heights.push_back(.5f);}
    for(int y=0;y<4;++y)for(int x=0;x<4;++x){int a=y*5+x,b=a+1,c=a+5,d=c+1;m.Triangles.push_back({a,b,d});m.Triangles.push_back({a,d,c});}
    if(k==0)for(auto &t:m.Triangles)std::swap(t[1],t[2]); // outward bottom cap
    for(int j=0;j<4;++j){m.BoundaryWidths[zTsUtils::makeEdgeKey(j,j+1)]=0;
      m.BoundaryWidths[zTsUtils::makeEdgeKey(20+j,21+j)]=0;
      m.BoundaryWidths[zTsUtils::makeEdgeKey(5*j,5*(j+1))]=1;
      m.BoundaryWidths[zTsUtils::makeEdgeKey(5*j+4,5*(j+1)+4)]=0;}
  }return layers;
}
bool sample(const zPoint &p,zPoint *uv) {
  if(p.x<-1e-6||p.x>4+1e-6||p.y<-1e-6||p.y>4+1e-6||p.z<-1e-6||p.z>2+1e-6)return false;
  if(uv)*uv=zPoint(p.x/4,p.y/4,p.z/2);return true;
}
// Independent geometric measure: jumps between adjacent unit face normals.
double bending(const Layer &m) {
  std::map<unsigned long long,std::vector<Eigen::Vector3d>> normals;
  for(const auto &t:m.Triangles){const auto a=zTsUtils::toEigenVec(m.Positions[t[0]]),b=zTsUtils::toEigenVec(m.Positions[t[1]]),c=zTsUtils::toEigenVec(m.Positions[t[2]]);
    const Eigen::Vector3d n=(b-a).cross(c-a).normalized();for(int j=0;j<3;++j)normals[zTsUtils::makeEdgeKey(t[j],t[(j+1)%3])].push_back(n);}
  double result=0;for(const auto &entry:normals)if(entry.second.size()==2)result+=(entry.second[0]-entry.second[1]).squaredNorm();return result;
}
void relaxation() {
  const Reference r=box();auto layers=stack(true);const auto original=layers;
  Settings settings;settings.MaxDeviation=.05;Report report;std::string error;
  require(relax(r,layers,settings,sample,report,error),error);
  require(report.AcceptedIterations>0&&report.MovedBoundaryVertices>0,"Boundary relaxation did not move any vertices");
  require(report.MaxDisplacement<=.05001&&report.MaxBoundaryError<1e-5,"Exceeded deviation or input-surface constraint");
  require(bending(layers[2])<bending(original[2]),"Intermediate face-normal variation was not reduced");
  for(size_t k=0;k<layers.size();++k){const auto &m=layers[k];
    require(m.Triangles==original[k].Triangles&&m.BoundaryWidths==original[k].BoundaryWidths,"Topology or side width groups changed");
    for(size_t i=0;i<m.Positions.size();++i){const auto &p=m.Positions[i];
      if(k==0||k+1==layers.size())require(norm(p,original[k].Positions[i])==0,"A fixed cap moved");
      require(std::isfinite(m.Heights[i])&&m.Heights[i]>0,"Invalid refreshed height");
      require(norm(m.UV[i],zPoint(p.x/4,p.y/4,p.z/2))<2e-6,"Transported UV was not refreshed");
    }
    for(const auto &entry:m.BoundaryWidths){const int a=int(entry.first>>32),b=int(entry.first&0xffffffffu);
      for(int id:{a,b}){const auto &p=m.Positions[id];const double side=(std::min)({std::abs(double(p.x)),std::abs(double(p.x)-4),std::abs(double(p.y)),std::abs(double(p.y)-4)});
        require(side<1e-6,"A boundary target left the reference box");}}
  }
  // Retain mandatory feature coordinates/identity through the actual PrintSyn
  // pipeline on the relaxed section, including smooth-normal interpolation.
  auto &middle=layers[2];zObjMesh chart,surface;zIntArray counts(middle.Triangles.size(),3),connects;
  for(const auto &t:middle.Triangles)connects.insert(connects.end(),t.begin(),t.end());
  auto flat=middle.Positions;for(auto &p:flat)p.z=0;
  zFnMesh(chart).create(flat,counts,connects);zFnMesh(surface).create(middle.Positions,counts,connects);
  const int ids[4]={6,8,18,16};zPointArray targets;for(int id:ids)targets.push_back(flat[id]);
  PrintPath path;zIntArray edges={0,1,1,2,2,3,3,0};zFnGraph(path.graph).create(targets,edges);
  path.widths.assign(4,.1f);path.featureFlags={PrintPath::AngleFeature|PrintPath::IntersectionFeature,PrintPath::AngleFeature,PrintPath::AngleFeature,PrintPath::AngleFeature};
  PrintPathArray paths2D{{path}},projected;PrintPath combined;PrintSynSetting ps;ps.ResampleLength=.3f;ps.FeatureAngleTol=100;ps.SeamOffset=0;ps.SeamSize=0;
  zObjMeshArray charts{chart},surfaces{surface};std::vector<zFloatArray> heights{middle.Heights};std::vector<zVectorArray> n{middle.Normals};
  require(zPrintSyn::computePipeline(charts,surfaces,nullptr,zPrintSyn::SeamMapping::WorldNearest,ps,paths2D,projected,combined,error,false,&heights,&n),error);
  zPointArray points;zFnGraph(combined.graph).getVertexPositions(points);
  for(int f=0;f<4;++f){bool found=false;for(size_t i=0;i<points.size();++i)if(norm(points[i],middle.Positions[ids[f]])<3e-6){
      found=true;require((combined.featureFlags[i]&path.featureFlags[f])==path.featureFlags[f],"Lost angle/intersection identity");
      require(norm(combined.normals[i],middle.Normals[ids[f]])<3e-6,"Projection replaced recovered normal with a face normal");}
    require(found,"Resampling or 3D projection lost a mandatory feature target");}
  const auto rejectNormals=[&](){require(!zPrintSyn::projectPrintPaths(paths2D,charts,surfaces,nullptr,projected,error,nullptr,false,&heights,&n),
      "Invalid smooth-normal table accepted");require(projected.empty(),"Invalid normals retained stale projected paths");};
  n[0].pop_back();rejectNormals();n[0]=middle.Normals;n[0][0]=zVector(0,0,0);rejectNormals();
  n[0][0]=zVector(std::numeric_limits<float>::quiet_NaN(),0,1);rejectNormals();
  n.clear();rejectNormals();
  std::cout<<"relaxation: normal variation "<<bending(original[2])<<" -> "<<bending(layers[2])
           <<", boundary moved="<<report.MovedBoundaryVertices<<", max movement="<<report.MaxDisplacement
           <<", boundary error="<<report.MaxBoundaryError<<'\n';
}
void constrainedInput() {
  auto layers=stack(true);const auto original=layers;Settings settings;Report report;std::string error;
  const auto missingUV=[](const zPoint &p,zPoint *uv){if(uv&&p.x>2.3)return false;return sample(p,uv);};
  require(relax(box(),layers,settings,missingUV,report,error),error);
  require(report.AcceptedIterations>0&&report.UVConstrainedVertices>0,"A local UV failure stopped all fairing");
  for(size_t k=1;k+1<layers.size();++k)for(size_t i=0;i<layers[k].Positions.size();++i)
    if(original[k].Positions[i].x>=3)require(norm(layers[k].Positions[i],original[k].Positions[i])==0,"A vertex without transported UV moved");
  layers=stack(false);layers[2].Positions=layers[1].Positions;
  require(!relax(box(),layers,settings,sample,report,error),"Touching layers were accepted");
  layers=stack(true);layers[2].Positions[0].x=.2;
  require(!relax(box(),layers,settings,sample,report,error),"Off-surface boundary was accepted");
}
void flatAndDisabled() {
  auto layers=stack(false);Settings settings;Report report;std::string error;
  require(relax(box(),layers,settings,sample,report,error),error);
  for(size_t k=0;k<layers.size();++k)for(size_t i=0;i<layers[k].Positions.size();++i){
    require(std::abs(layers[k].Positions[i].z-.5*k)<2e-6,"Flat layer moved out of plane");
    require(std::abs(layers[k].Heights[i]-.5)<2e-5,"Analytic normal layer spacing changed");}
  layers=stack(true);const auto original=layers;settings.Iterations=0;
  require(relax(box(),layers,settings,sample,report,error),error);
  for(size_t k=0;k<layers.size();++k){require(layers[k].Heights==original[k].Heights,"Disabled relaxation changed heights");
    for(size_t i=0;i<layers[k].Positions.size();++i)require(norm(layers[k].Positions[i],original[k].Positions[i])==0,"Disabled relaxation changed positions");}
  settings.Iterations=-1;require(!relax(box(),layers,settings,sample,report,error),"Negative iterations accepted");
  settings.Iterations=8;settings.Strength=0;require(!relax(box(),layers,settings,sample,report,error),"Zero strength accepted");
  settings.Strength=.25;settings.MaxDeviation=-1;require(!relax(box(),layers,settings,sample,report,error),"Negative deviation accepted");
}
}
int main() {try{relaxation();flatAndDisabled();constrainedInput();return 0;}catch(const std::exception &e){std::cerr<<e.what()<<'\n';return 1;}}
