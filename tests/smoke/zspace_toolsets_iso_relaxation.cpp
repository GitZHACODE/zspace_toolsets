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
  // Retain geometric corner coordinates through the actual PrintSyn pipeline
  // on the relaxed section, including smooth-normal interpolation. Historical
  // intersection metadata is deliberately cleared by angle re-detection.
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
      found=true;require(combined.featureFlags[i]==0,"An angle-only sampling support became a junction feature");
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
void heightDomain() {
  auto reference=box();for(auto &p:reference.Positions)p.z*=20;
  auto layers=stack(false);const double z[5]={0,5,10,27,40};
  for(int k=0;k<5;++k)for(auto &p:layers[k].Positions)p.z=float(z[k]);
  const auto original=layers;
  const auto inside=[](const zPoint &p,zPoint *uv){if(p.x<0||p.x>4||p.y<0||p.y>4||p.z<0||p.z>40)return false;
    if(uv)*uv=zPoint(p.x/4,p.y/4,p.z/40);return true;};
  PrintHeightOptimizationSettings settings;settings.Enabled=true;settings.Iterations=16;
  HeightReport report;std::string error;
  require(optimizeHeightDomain(reference,layers,settings,inside,false,report,error),error);
  std::cout<<"height domain "<<report.InitialMinimum<<","<<report.InitialMaximum<<" -> "<<report.Minimum<<","<<report.Maximum<<" iterations="<<report.AcceptedIterations<<'\n';
  require(report.Converged,"Feasible planar height domain did not converge");
  require(report.DegreesOfFreedom>0&&report.DegreesOfFreedom<=75,
          "Height correction did not use a reduced scalar surface system");
  require(report.InnerIterations>0&&report.LineSearchTrials>=report.AcceptedIterations,
          "Height solver iteration diagnostics are incomplete");
  require(report.MinimumDirectedGap>0&&report.MaximumDirectedGap>=report.MinimumDirectedGap,
          "Directed geometric gap diagnostics are invalid");
  for(size_t k=0;k<layers.size();++k){require(layers[k].Triangles==original[k].Triangles,"Height correction changed topology");
    if(k==0||k+1==layers.size())for(size_t i=0;i<layers[k].Positions.size();++i)require(norm(layers[k].Positions[i],original[k].Positions[i])==0,"Height optimizer moved a cap");
    if(k>0)for(float h:layers[k].Heights)require(h>=6&&h<=12,"Successful TCP previous-layer height outside supplied domain");
    require(bending(layers[k])<1e-5,"Constant-height correction introduced surface ripples");
  }
  auto repeat=original;HeightReport parallelReport;
  require(optimizeHeightDomain(reference,repeat,settings,inside,true,parallelReport,error),error);
  for(size_t k=0;k<layers.size();++k){require(layers[k].Heights==repeat[k].Heights,"Parallel height oracle changed output");
    for(size_t i=0;i<layers[k].Positions.size();++i)require(norm(layers[k].Positions[i],repeat[k].Positions[i])==0,"Parallel height optimizer changed geometry");}
  std::vector<double> levels;require(redistributeIsoValues(original,settings,levels,false,error),error);
  require(levels.front()==0&&levels.back()==1,"Iso schedule changed cap values");
  for(size_t k=1;k<levels.size();++k)require(levels[k]>levels[k-1],"Height iso schedule is not ordered");
  settings.Maximum=settings.Minimum;require(!optimizeHeightDomain(reference,repeat,settings,inside,false,report,error),"Degenerate height target accepted");
  settings.Maximum=12;settings.Minimum=-1;require(!optimizeHeightDomain(reference,repeat,settings,inside,false,report,error),"Negative height target accepted");
  settings.Minimum=std::numeric_limits<double>::quiet_NaN();require(!optimizeHeightDomain(reference,repeat,settings,inside,false,report,error),"Nonfinite target accepted");
  settings.Minimum=6;settings.MaxDisplacement=.01;repeat=original;
  require(optimizeHeightDomain(reference,repeat,settings,inside,false,report,error),error);
  require(!report.Converged&&report.WorstViolation>1,"Infeasible displacement budget falsely reported compliance");
  require(report.MaxDisplacement<=.01001,"Exceeded height correction movement budget");
  PrintHeightOptimizationSettings exact;exact.Enabled=true;
  for(const auto &bounds:std::vector<std::pair<double,double>>{{.5,.6},{.4,.5}}){exact.Minimum=bounds.first;exact.Maximum=bounds.second;
    auto flat=stack(false);const auto before=flat;
    require(optimizeHeightDomain(box(),flat,exact,sample,false,report,error),error);
    require(report.Converged&&report.AcceptedIterations==0&&report.Minimum==.5&&report.Maximum==.5,"Exact inclusive planar bound was changed/rejected");
    std::vector<double> schedule;require(redistributeIsoValues(before,exact,schedule,false,error),error);
    for(size_t k=0;k<schedule.size();++k)require(std::abs(schedule[k]-.25*k)<1e-12,"Feasible exact-bound schedule was perturbed");
  }
  exact.Minimum=.52;exact.Maximum=.6;exact.Tolerance=.05;auto tolerant=stack(false);
  require(optimizeHeightDomain(box(),tolerant,exact,sample,false,report,error),error);
  require(report.Converged&&report.AcceptedIterations==0&&report.Minimum==.5&&report.Maximum==.5,"Fabrication tolerance changed/clamped measured heights");
  exact.Tolerance=-.01;require(!optimizeHeightDomain(box(),tolerant,exact,sample,false,report,error),"Negative tolerance accepted");
  exact.Tolerance=std::numeric_limits<double>::quiet_NaN();require(!optimizeHeightDomain(box(),tolerant,exact,sample,false,report,error),"Nonfinite tolerance accepted");
  exact.Enabled=false;exact.Minimum=std::numeric_limits<double>::quiet_NaN();auto disabled=original;
  require(optimizeHeightDomain(reference,disabled,exact,inside,false,report,error),error);
  for(size_t k=0;k<disabled.size();++k){require(disabled[k].Heights==original[k].Heights,"Disabled height optimizer changed heights");
    for(size_t i=0;i<disabled[k].Positions.size();++i)require(norm(disabled[k].Positions[i],original[k].Positions[i])==0,"Disabled height optimizer moved geometry");}
}
void qualityRepair() {
  auto layers=stack(false);const auto original=layers;
  for(int k=1;k<4;++k){layers[k].Positions[12].x=1.001f;layers[k].UV[12].x=1.001f/4;}
  Settings settings;settings.Iterations=2;settings.MinEdgeLength=.02;settings.MaxDeviation=.02;
  Report report;std::string error;require(relax(box(),layers,settings,sample,report,error),error);
  for(int k=0;k<5;++k){const auto &m=layers[k];
    if(k==0||k==4){require(m.Triangles==original[k].Triangles,"Quality repair changed a fixed cap");for(size_t i=0;i<m.Positions.size();++i)require(norm(m.Positions[i],original[k].Positions[i])==0,"Quality repair moved a cap");}
    require(bending(m)<1e-5,"Quality repair bent a planar section");
    for(float h:m.Heights)require(std::abs(h-.5)<1e-5,"Quality repair changed analytic spacing");
    for(const auto &t:m.Triangles)for(int j=0;j<3;++j)require(norm(m.Positions[t[j]],m.Positions[t[(j+1)%3]])>.01,"Tiny interior edge survived quality repair");
    require(m.BoundaryWidths.size()==original[k].BoundaryWidths.size(),"Quality repair lost semantic boundary edges");
    for(const auto &e:m.BoundaryWidths)for(int id:{int(e.first>>32),int(e.first&0xffffffffu)}){const auto &p=m.Positions[id];
      require((std::min)({std::abs(p.x),std::abs(p.x-4),std::abs(p.y),std::abs(p.y-4)})<1e-6,"Quality boundary left input surface");}
  }
  settings.MinEdgeLength=-1;require(!relax(box(),layers,settings,sample,report,error),"Negative repair length accepted");
}
void heightGeometryCases(){
  for(bool narrow:{false,true}){
    auto reference=box();auto layers=stack(false);const double level[5]={0,5,17,29,40};
    for(auto &p:reference.Positions){p.z*=20;if(narrow)p.x*=.025f;else p.z*=1+.002f*p.x;}
    for(int k=0;k<5;++k)for(auto &p:layers[k].Positions){p.z=float(level[k]);if(narrow)p.x*=.025f;else p.z*=1+.002f*p.x;}
    const auto original=layers;const auto inside=[&](const zPoint &p,zPoint *uv){const double width=narrow?.1:4,top=narrow?40:40*(1+.002*p.x);
      if(p.x<-1e-6||p.x>width+1e-6||p.y<-1e-6||p.y>4+1e-6||p.z<-1e-6||p.z>top+1e-5)return false;
      if(uv)*uv=zPoint(p.x/width,p.y/4,p.z/top);return true;};
    PrintHeightOptimizationSettings settings;settings.Enabled=true;settings.Tolerance=.05;settings.Iterations=32;
    // Printed supports lie inside the trim. A steeper 0.8 top slope on this
    // 4 mm wide, 40 mm tall prism has no previous-layer ray support; the old
    // plane-extrapolated fixture was not a physically feasible TCP case.
    const int materialIds[]={6,7,8,11,12,13,16,17,18};
    std::vector<zPointArray> targets(layers.size()-1);
    for(size_t k=1;k<layers.size();++k)for(int id:materialIds)targets[k-1].push_back(layers[k].Positions[id]);
    HeightReport report;std::string error;require(optimizeHeightDomain(reference,layers,settings,inside,false,report,error,&targets),error);
    std::cout<<(narrow?"narrow":"tilted")<<" height: domain="<<report.Minimum<<","<<report.Maximum
             <<" directed="<<report.MinimumDirectedGap<<","<<report.MaximumDirectedGap
             <<" bending="<<report.InitialBendingEnergy<<" -> "<<report.FinalBendingEnergy
             <<" dofs="<<report.DegreesOfFreedom<<" ms="<<report.Milliseconds<<'\n';
    require(report.Converged,narrow?"Narrow stack failed height convergence":"Tilted-cap stack failed height convergence");
    for(size_t k=0;k<layers.size();++k){if(k==0||k+1==layers.size())for(size_t i=0;i<layers[k].Positions.size();++i)
        require(norm(layers[k].Positions[i],original[k].Positions[i])==0,"Geometry-case cap moved");
    }
    std::vector<zVectorArray> targetNormals(targets.size());
    for(size_t k=1;k<layers.size();++k){targets[k-1].clear();for(int id:materialIds){
      targets[k-1].push_back(layers[k].Positions[id]);targetNormals[k-1].push_back(layers[k].Normals[id]);}}
    std::vector<zFloatArray> actual;require(measureTargetHeights(layers,targets,targetNormals,actual,error,false),error);
    for(const auto &layer:actual)for(float height:layer)require(height>=5.95&&height<=12.05,"Geometry case contains an unmeasured actual TCP-gap violation");
  }
  auto bad=stack(false);bad[2].Triangles[0][1]=bad[2].Triangles[0][0];PrintHeightOptimizationSettings settings;settings.Enabled=true;
  HeightReport report;std::string error;require(!optimizeHeightDomain(box(),bad,settings,sample,false,report,error),"Degenerate height-support triangle accepted");
}
void heightMissingTransport() {
  auto reference=box();for(auto &p:reference.Positions)p.z*=20;
  auto layers=stack(false);const double level[5]={0,5,10,27,40};
  for(int k=0;k<5;++k)for(auto &p:layers[k].Positions)p.z=float(level[k]);
  const auto original=layers;
  const auto inside=[](const zPoint &p,zPoint *uv){
    if(p.x<0||p.x>4||p.y<0||p.y>4||p.z<0||p.z>40)return false;
    if(uv){if(p.x>2.3)return false;*uv=zPoint(p.x/4,p.y/4,p.z/40);}return true;};
  PrintHeightOptimizationSettings settings;settings.Enabled=true;settings.Iterations=4;
  HeightReport report;std::string error;
  require(optimizeHeightDomain(reference,layers,settings,inside,false,report,error),error);
  for(size_t k=1;k+1<layers.size();++k)for(size_t i=0;i<layers[k].Positions.size();++i)
    if(original[k].Positions[i].x>=3)
      require(norm(layers[k].Positions[i],original[k].Positions[i])==0,
              "Height correction moved a point without transported UV");
  require(!report.Converged,"Unavailable transported coordinates were reported as a feasible height correction");
}
void inactiveHeightLayers() {
  auto layers=stack(false);
  for(auto &p:layers[3].Positions)p.z+=float(.02*std::sin(p.x*3.141592653589793/4)*std::sin(p.y*3.141592653589793/4));
  const auto original=layers;
  std::vector<zPointArray> targets(layers.size()-1);targets[0]=layers[1].Positions;
  PrintHeightOptimizationSettings settings;settings.Enabled=true;settings.Minimum=.6;settings.Maximum=.75;settings.Tolerance=.001;
  HeightReport report;std::string error;
  require(optimizeHeightDomain(box(),layers,settings,sample,false,report,error,&targets),error);
  require(report.Converged,"Local cap-to-first-layer correction did not converge");
  for(size_t k=2;k<layers.size();++k)for(size_t i=0;i<layers[k].Positions.size();++i)
    require(norm(layers[k].Positions[i],original[k].Positions[i])==0,
            "A height batch moved an unrelated layer without active supports");
  auto parallel=original;HeightReport parallelReport;
  require(optimizeHeightDomain(box(),parallel,settings,sample,true,parallelReport,error,&targets),error);
  require(parallelReport.Converged,"Parallel local height correction did not converge");
  for(size_t k=0;k<layers.size();++k)for(size_t i=0;i<layers[k].Positions.size();++i)
    require(norm(layers[k].Positions[i],parallel[k].Positions[i])==0,
            "Local height correction changed under parallel execution");
}
void previousLayerMeasurement() {
  for(bool rotated:{false,true}){
    auto reference=box();auto layers=stack(false);const double levels[5]={0,.2,.9,1.2,2};
    for(int k=0;k<5;++k)for(auto &p:layers[k].Positions)p.z=float(levels[k]);
    const auto rotate=[](const zPoint &p){return zPoint(-p.z,p.y,p.x);};
    if(rotated){for(auto &p:reference.Positions)p=rotate(p);
      for(auto &layer:layers)for(auto &p:layer.Positions)p=rotate(p);}
    const auto inside=[&](const zPoint &p,zPoint *uv){return sample(rotated?zPoint(p.z,p.y,-p.x):p,uv);};
    const auto original=layers;PrintHeightOptimizationSettings settings;settings.Enabled=true;settings.Minimum=.1;settings.Maximum=1;
    HeightReport report;std::string error;require(optimizeHeightDomain(reference,layers,settings,inside,false,report,error),error);
    require(report.Converged&&report.AcceptedIterations==0,"Feasible previous-layer gaps were changed");
    for(size_t k=1;k<layers.size();++k)for(size_t i=0;i<layers[k].Heights.size();++i){
      require(std::abs(layers[k].Heights[i]-(levels[k]-levels[k-1]))<3e-6,
              "TCP height averaged neighboring gaps or assumed world +Z");
      require(norm(layers[k].Positions[i],original[k].Positions[i])==0,
              "Height measurement moved an existing TCP surface");}
    require(std::abs(report.Minimum-.2)<3e-6&&std::abs(report.Maximum-.8)<3e-6,
            "TCP height report included substrate proxy or omitted the top TCP");
  }
}
void exactTCPRays() {
  auto layers=stack(false);for(auto &layer:layers)for(auto &p:layer.Positions)p.x*=10;
  std::vector<zPointArray> targets(layers.size()-1);std::vector<zVectorArray> directions(targets.size());
  targets[0]={{20,2,.5}};directions[0]={{float(std::sqrt(.99)),0,.1f}};
  std::vector<zFloatArray> measured,parallel;std::string error;
  require(measureTargetHeights(layers,targets,directions,measured,error,false),error);
  require(measured[0].size()==1&&std::abs(measured[0][0]-5)<2e-5,
          "TCP ray used nodal spacing or stopped at four nearest distances");
  require(measureTargetHeights(layers,targets,directions,parallel,error,true),error);
  require(measured==parallel,"Parallel exact TCP heights changed");
  directions[0][0]=zVector(0,0,.25f);
  require(measureTargetHeights(layers,targets,directions,measured,error,false)&&std::abs(measured[0][0]-.5)<1e-6,
          "A scaled valid TCP axis changed its physical height");
  directions[0][0]=zVector(1,0,0);
  require(!measureTargetHeights(layers,targets,directions,measured,error,false)&&measured.empty(),
          "Unsupported horizontal TCP ray was extrapolated or retained stale heights");
  directions[0][0]=zVector(0,0,-1);
  require(!measureTargetHeights(layers,targets,directions,measured,error,false),"Reversed TCP ray acquired previous support");
  directions[0][0]=zVector(0,0,1);targets[0][0].z=.6f;
  require(!measureTargetHeights(layers,targets,directions,measured,error,false),"Off-section TCP target accepted");
}
}
int main() {try{relaxation();flatAndDisabled();constrainedInput();heightDomain();qualityRepair();heightGeometryCases();heightMissingTransport();inactiveHeightLayers();previousLayerMeasurement();exactTCPRays();return 0;}catch(const std::exception &e){std::cerr<<e.what()<<'\n';return 1;}}
