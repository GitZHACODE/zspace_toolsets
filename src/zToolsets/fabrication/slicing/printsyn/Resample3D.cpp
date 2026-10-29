#include <zspace/zToolsets/fabrication/slicing/printsyn/zPrintSyn.h>
#include <zspace/zToolsets/fabrication/slicing/utils/Parallel.h>
#include "../utils/SurfaceQuery.h"
#include "../utils/Interpolation.h"
#include <algorithm>
#include <cmath>
#include <sstream>

namespace zSpace {
namespace {
using V=Eigen::Vector3d;
using zTsUtils::toEigenVec;
zPoint point(const V &v){return zPoint(v.x(),v.y(),v.z());}

bool resampleMapped(PrintPath &path2D,PrintPath &path3D,const zPointArray &chart,
                    const zTsDetail::SurfaceQuery &surface,double spacing,
                    const zFloatArray *heights,const zVectorArray *normals,std::string &error){
  zTsUtils::ScopedSlicingTiming methodTiming("PrintSyn.resampleMapped", true);
  zPointArray p,q;zFnGraph(path3D.graph).getVertexPositions(p);zFnGraph(path2D.graph).getVertexPositions(q);
  const size_t n=p.size();if(n==0)return true;
  const auto valid=[&](size_t count){return count==0||count==n;};
  if(q.size()!=n||!valid(path2D.normals.size())||!valid(path3D.widths.size())||!valid(path3D.heights.size())||!valid(path3D.normals.size())||!valid(path3D.trackIds.size())||!valid(path3D.featureFlags.size())){error="Mapped attribute count mismatch.";return false;}
  std::vector<zTsUtils::GraphPolyline> curves;if(!zTsUtils::graphPolylines(path3D.graph,curves)||curves.size()!=1||curves[0].Vertices.size()!=n){error="Expected one connected unbranched mapped curve.";return false;}
  const auto &sequence=curves[0].Vertices;const bool closed=curves[0].Closed;
  const int segments=closed?int(n):int(n)-1;std::vector<double> arc(segments+1,0);
  for(int j=0;j<segments;++j)arc[j+1]=arc[j]+(toEigenVec(p[sequence[(j+1)%n]])-toEigenVec(p[sequence[j]])).norm();
  if(!std::isfinite(arc.back())||arc.back()<=0){error="Invalid mapped arc length.";return false;}
  std::vector<int> anchors{0};for(int j=1;j<int(n);++j)if((!path3D.featureFlags.empty()&&path3D.featureFlags[sequence[j]])||(!closed&&j==int(n)-1))anchors.push_back(j);
  if(closed)anchors.push_back(int(n));
  PrintPath result3,result2;zPointArray positions,charts;zIntArray edges;
  const bool hasH=!path3D.heights.empty(),hasN=!path3D.normals.empty(),hasT=!path3D.trackIds.empty();
  const auto value=[](const zFloatArray &values,int a,int b,double t,float fallback){return values.empty()?fallback:float(values[a]+t*(values[b]-values[a]));};
  const auto append=[&](int a,int b,double t,bool mandatory,int flags,const zTsDetail::SurfaceQuery::Hit *located=nullptr){
    V target=(1-t)*toEigenVec(p[a])+t*toEigenVec(p[b]);V flat=(1-t)*toEigenVec(q[a])+t*toEigenVec(q[b]);
    float height=value(path3D.heights,a,b,t,0);V normal=hasN?((1-t)*toEigenVec(path3D.normals[a])+t*toEigenVec(path3D.normals[b])).eval():V(0,0,1);
    if(!mandatory){zTsDetail::SurfaceQuery::Hit hit;if(located)hit=*located;else if(!surface.nearest(point(target),hit)){error="Cannot place new sample on its section.";return false;}
      target=hit.Position;const auto &tri=surface.triangle(hit.Triangle);flat.setZero();
      for(int k=0;k<3;++k)flat+=hit.Weights[k]*toEigenVec(chart[tri.Ids[k]]);
      if(heights&&!heights->empty())height=float(zTsDetail::interpolateScalar(
        (*heights)[tri.Ids[0]],(*heights)[tri.Ids[1]],(*heights)[tri.Ids[2]],
        hit.Weights[0],hit.Weights[1],hit.Weights[2]));
      if(normals&&!normals->empty()){normal.setZero();for(int k=0;k<3;++k)normal+=hit.Weights[k]*toEigenVec((*normals)[tri.Ids[k]]);}else normal=tri.Normal;
    }
    if(normal.norm()>0)normal.normalize();
    positions.push_back(point(target));charts.push_back(point(flat));result3.widths.push_back(value(path3D.widths,a,b,t,.04f));
    if(!path2D.normals.empty()){V chartNormal=(1-t)*toEigenVec(path2D.normals[a])+t*toEigenVec(path2D.normals[b]);if(chartNormal.norm()>0)chartNormal.normalize();result2.normals.push_back(point(chartNormal));}
    if(hasH)result3.heights.push_back(height);if(hasN)result3.normals.push_back(point(normal));if(hasT)result3.trackIds.push_back(value(path3D.trackIds,a,b,t,0));result3.featureFlags.push_back(flags);return true;
  };
  const auto flags=[&](int id){return path3D.featureFlags.empty()?0:path3D.featureFlags[id];};
  if(!append(sequence[0],sequence[0],0,true,flags(sequence[0])))return false;
  for(size_t s=0;s+1<anchors.size();++s){const int begin=anchors[s],end=anchors[s+1];const double length=arc[end]-arc[begin];
    int intervals=(std::max)(1,int(std::ceil(length/spacing)));if(closed&&anchors.size()==2)intervals=(std::max)(3,intervals);
    const int startId=sequence[begin%n],endId=sequence[end%n];
    for(int j=1;j<intervals;++j){const double d=arc[begin]+length*j/intervals;
      // A tiny span can have mandatory endpoints closer than the requested
      // spacing. Keep both; no ordinary target may crowd either endpoint.
      if(d-arc[begin]<.5*spacing-1e-8||arc[end]-d<.5*spacing-1e-8)continue;
      int seg=int(std::upper_bound(arc.begin()+begin,arc.begin()+end+1,d)-arc.begin())-1;seg=(std::max)(begin,(std::min)(seg,end-1));
      const int a=sequence[seg%n],b=sequence[(seg+1)%n];const double den=arc[seg+1]-arc[seg],t=den>0?(d-arc[seg])/den:0;
      const V target=(1-t)*toEigenVec(p[a])+t*toEigenVec(p[b]);zTsDetail::SurfaceQuery::Hit hit;
      if(!surface.nearest(point(target),hit)){error="Cannot project new 3D sample.";return false;}
      if((hit.Position-toEigenVec(p[startId])).norm()<.5*spacing-1e-5||(hit.Position-toEigenVec(p[endId])).norm()<.5*spacing-1e-5)continue;
      if(!append(a,b,t,false,0,&hit))return false;
    }
    if(!(closed&&end==int(n))&&!append(endId,endId,0,true,flags(endId)))return false;
  }
  if(positions.size()<(closed?3u:2u)){error="Feature clearance leaves too few targets for the curve topology.";return false;}
  for(int j=0;j<int(positions.size())-1;++j){edges.push_back(j);edges.push_back(j+1);}if(closed){edges.push_back(int(positions.size())-1);edges.push_back(0);}
  zFnGraph(result3.graph).create(positions,edges);zFnGraph(result2.graph).create(charts,edges);
  result2.widths=result3.widths;result2.heights=result3.heights;result2.trackIds=result3.trackIds;result2.featureFlags=result3.featureFlags;
  path3D=std::move(result3);path2D=std::move(result2);return true;
}
}

bool zPrintSyn::resamplePrintPaths3D(PrintPathArray &paths2D,PrintPathArray &paths3D,
    const zObjMeshArray &charts2D,const zObjMeshArray &sections3D,const PrintSynSetting &settings,
    std::string &status,bool useParallel,const zFloatArray *vertexHeights,
    const std::vector<zFloatArray> *perLayerHeights,const std::vector<zVectorArray> *perLayerNormals){
  zTsUtils::ScopedSlicingTiming methodTiming("PrintSyn.resamplePrintPaths3D");
  const int n=int(paths3D.size());if(n==0||paths2D.size()!=size_t(n)||charts2D.size()<size_t(n)||sections3D.size()<size_t(n)||!std::isfinite(settings.ResampleLength)||settings.ResampleLength<=0){status="Invalid second 3D sampling inputs.";return false;}
  if((perLayerHeights&&perLayerHeights->size()<size_t(n))||(perLayerNormals&&perLayerNormals->size()<size_t(n))){status="Missing mapped layer attribute tables.";return false;}
  std::vector<std::string> errors(n);zTsDetail::ParallelErrors failures(n);
  #pragma omp parallel for schedule(dynamic) if(useParallel)
  for(int k=0;k<n;++k){try{
      zPointArray chart;zFnMesh(const_cast<zObjMesh &>(charts2D[k])).getVertexPositions(chart);
      if(chart.size()!=size_t(zFnMesh(const_cast<zObjMesh &>(sections3D[k])).numVertices())){errors[k]="Chart/section vertex count mismatch.";continue;}
      zTsDetail::SurfaceQuery surface(const_cast<zObjMesh &>(sections3D[k]));
      const zFloatArray *h=perLayerHeights?&(*perLayerHeights)[k]:vertexHeights;const zVectorArray *normals=perLayerNormals?&(*perLayerNormals)[k]:nullptr;
      if(paths2D[k].size()!=paths3D[k].size()||(h&&!h->empty()&&h->size()!=chart.size())||(normals&&!normals->empty()&&normals->size()!=chart.size())){errors[k]="Mapped layer attribute table mismatch.";continue;}
      for(size_t p=0;p<paths3D[k].size();++p)if(!resampleMapped(paths2D[k][p],paths3D[k][p],chart,surface,settings.ResampleLength,h,normals,errors[k])){errors[k]="path "+std::to_string(p)+": "+errors[k];break;}
    }catch(...){failures.capture(k);}}
  failures.rethrow();for(int k=0;k<n;++k)if(!errors[k].empty()){status="Second 3D sampling failed on layer "+std::to_string(k)+", "+errors[k];return false;}
  status="Resampled mapped 3D paths between mandatory features.";return true;
}
}
