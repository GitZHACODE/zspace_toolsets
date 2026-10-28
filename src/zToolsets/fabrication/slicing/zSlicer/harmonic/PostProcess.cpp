#ifndef EIGEN_DONT_PARALLELIZE
#define EIGEN_DONT_PARALLELIZE
#endif

#include <zspace/zToolsets/fabrication/slicing/zSlicer/HarmonicSlicer.h>

#ifndef ZSPACE_INLINE
#define ZSPACE_INLINE ZSPACE_TOOLSETS_INLINE
#endif

#include <sstream>
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <zspace/zToolsets/fabrication/slicing/utils/Parallel.h>
#include <string>
#include <vector>
#include "../../utils/SurfaceQuery.h"
#include "../../utils/Interpolation.h"

namespace zSpace {

using namespace zTsUtils;

namespace {
using V=Eigen::Vector3d;
struct SurfaceFootprint {
  std::vector<zPointArray> Targets;
  double Minimum=std::numeric_limits<double>::infinity(),Maximum=0;
  int Violations=0;
};
// Within one chart triangle both the mapped position and the existing nodal
// height interpolation are affine. A convex hull of the original curve's
// material supports therefore bounds EVERY point of that footprint, without
// adding any fabrication targets or retaining every dense SDF vertex.
bool surfaceFootprint(const PrintPathArray &curves,const zObjMeshArray &charts,
    const zObjMeshArray &sections,const std::vector<zFloatArray> &heights,
    const PrintHeightOptimizationSettings &domain,bool parallel,SurfaceFootprint &out,std::string &error){
  out=SurfaceFootprint();out.Targets.resize(curves.size());std::vector<double> low(curves.size(),out.Minimum),high(curves.size(),0);
  std::vector<int> bad(curves.size(),0);std::vector<std::string> errors(curves.size());zTsDetail::ParallelErrors failures(curves.size());
  #pragma omp parallel for schedule(dynamic) if(parallel)
  for(int k=0;k<int(curves.size());++k){try{
    zTsDetail::SurfaceQuery query(const_cast<zObjMesh&>(charts[k]));zPointArray world,chart;
    zFnMesh(const_cast<zObjMesh&>(sections[k])).getVertexPositions(world);zFnMesh(const_cast<zObjMesh&>(charts[k])).getVertexPositions(chart);
    double scale=1;for(const auto &p:chart)scale=(std::max)(scale,toEigenVec(p).cwiseAbs().maxCoeff());
    const double tolerance=(std::max)(1e-7,8*std::numeric_limits<float>::epsilon()*scale);
    std::map<int,std::vector<V>> supports;
    for(const auto &curve:curves[k]){zPointArray p;zIntArray edges;zFnGraph fn(const_cast<zObjGraph&>(curve.graph));fn.getVertexPositions(p);fn.getEdgeData(edges);
      for(size_t e=0;e+1<edges.size();e+=2){const auto &a=p[edges[e]],&b=p[edges[e+1]];std::vector<double> breaks;
        if(!query.chartSegment(a,b,tolerance,breaks)){errors[k]="Continuous height footprint leaves its chart or has overlapping supports.";break;}
        for(double t:breaks){const V q=(1-t)*toEigenVec(a)+t*toEigenVec(b);zTsDetail::SurfaceQuery::Hit hit;
          if(!query.uniqueChartHit(zPoint(q.x(),q.y(),q.z()),hit,tolerance)){errors[k]="Cannot map continuous height footprint.";break;}
          V w=hit.Weights.cwiseMax(0);w/=w.sum();const auto &ids=query.triangle(hit.Triangle).Ids;
          const double h=zTsDetail::interpolateScalar(heights[k][ids[0]],heights[k][ids[1]],heights[k][ids[2]],w[0],w[1],w[2]);
          if(!std::isfinite(h)){errors[k]="Nonfinite continuous-curve height.";break;}
          low[k]=(std::min)(low[k],h);high[k]=(std::max)(high[k],h);
          if(h<domain.Minimum-domain.Tolerance||h>domain.Maximum+domain.Tolerance)++bad[k];supports[hit.Triangle].push_back(w);}
        if(!errors[k].empty())break;}
      if(!errors[k].empty())break;}
    if(!errors[k].empty())continue;
    const auto cross=[](const V &a,const V &b,const V &c){return (b[1]-a[1])*(c[2]-a[2])-(b[2]-a[2])*(c[1]-a[1]);};
    for(auto &entry:supports){auto &points=entry.second;
      // A rebuilt chart moves the material footprint slightly. Cover a small
      // TWO-dimensional neighbourhood of the curve, clipped to its face.
      // A single sample towards the centroid does not cover that region.
      const auto original=points;points.clear();points.reserve(3*original.size());
      for(const auto &w:original)for(int j=0;j<3;++j)points.push_back(.88*w+.12*V::Unit(j));
      std::sort(points.begin(),points.end(),[](const V &a,const V &b){return a[1]<b[1]||(a[1]==b[1]&&a[2]<b[2]);});
      points.erase(std::unique(points.begin(),points.end(),[](const V &a,const V &b){return (a-b).norm()<1e-10;}),points.end());
      std::vector<V> hull;for(const auto &p:points){while(hull.size()>1&&cross(hull[hull.size()-2],hull.back(),p)<=1e-12)hull.pop_back();hull.push_back(p);}
      const size_t lower=hull.size();if(points.size()>1){for(int i=int(points.size())-2;i>=0;--i){while(hull.size()>lower&&cross(hull[hull.size()-2],hull.back(),points[i])<=1e-12)hull.pop_back();hull.push_back(points[i]);}hull.pop_back();}
      const auto &ids=query.triangle(entry.first).Ids;for(const auto &w:hull){V p=V::Zero();for(int j=0;j<3;++j)p+=w[j]*toEigenVec(world[ids[j]]);out.Targets[k].emplace_back(p.x(),p.y(),p.z());}}
  }catch(...){failures.capture(k);}}
  failures.rethrow();for(size_t k=0;k<curves.size();++k){if(!errors[k].empty()){error=errors[k];return false;}
    out.Minimum=(std::min)(out.Minimum,low[k]);out.Maximum=(std::max)(out.Maximum,high[k]);out.Violations+=bad[k];}
  if(!std::isfinite(out.Minimum)){error="No continuous print-curve supports.";return false;}return true;
}
}

//--------------------------
//---- PUBLIC METHODS: zSlicer
//--------------------------

ZSPACE_INLINE bool HarmonicSlicer::computeSlicer(bool useParallel) {
  const bool success = computeSlicerImpl(useParallel);
  if (success) zTsUtils::printSuccess("[HarmonicSlicer::computeSlicer] Success");
  else zTsUtils::printError("[HarmonicSlicer::computeSlicer] Failed");
  return success;
}

ZSPACE_INLINE bool HarmonicSlicer::computeSlicerImpl(bool useParallel) {
  if (!computeBasicStructure(useParallel)) return false;
  if (!computeHarmonicField()) return false;
  if (!computeVolumetricParameterization(useParallel)) return false;
  if (!computeContours(useParallel)) return false;
  return computeField(useParallel);
}

// Every layer has its own path/width/track field on the shared o_ScalarField
// grid, in its own o_IsoMeshes2D chart. The layers share no topology, so the seam is
// carried between them in world coordinates (zPrintSyn's charts/sections
// overload of alignStartPoints()), and every layer's 2D paths go to world
// through its own o_IsoMeshes2D / o_IsoMeshes pair, interpolating that layer's
// vertex heights after resampling and seam alignment.
//
// Stops at the first zPrintSyn step that fails or reports a problem on any
// layer: o_PrintPaths2D keeps the paths as that step left them, projected to
// o_PrintPaths3D for inspection (no o_ResultantPrintPath), _Status names the
// step and its layers, and false is returned.
ZSPACE_INLINE bool HarmonicSlicer::computePostProcess(const PrintSynSetting &settings, bool useParallel) {
  const bool success = computePostProcessImpl(settings, useParallel);
  if (success) updatePrintStatistics(true);
  else resetPrintStatistics();
  if (success) zTsUtils::printSuccess("[HarmonicSlicer::computePostProcess] Success");
  else zTsUtils::printError("[HarmonicSlicer::computePostProcess] Failed");
  return success;
}


bool HarmonicSlicer::rebuildSurfaceCharts(bool parallel){
  const int count=int(o_IsoMeshes.size());std::vector<char> valid(count,0);zTsDetail::ParallelErrors failures(count);
  #pragma omp parallel for schedule(dynamic) if(parallel)
  for(int k=0;k<count;++k){try{zPointArray world,flat;
    if(!unrollIsoMesh(o_IsoMeshes[k],o_IsoMeshUV[k],o_IsoMeshes2D[k],world,flat))continue;
    zIntArray connects,counts;zFnMesh(o_IsoMeshes[k]).getPolygonData(connects,counts);std::vector<std::string> debug;
    createIsoBoundaries(world,flat,connects,counts,_HarmonicIsoEdgeWidthIds[k],o_IsoMeshBoundary[k],o_IsoMeshBoundary2D[k],o_IsoMeshBoundaryWidthId[k],debug);decodeBoundaryGroups(k);valid[k]=1;
  }catch(...){failures.capture(k);}}
  failures.rethrow();if(std::find(valid.begin(),valid.end(),char(0))!=valid.end()){
    _Status="Surface feedback could not unfold every section.";return false;}
  return computeField(parallel);
}

ZSPACE_INLINE bool HarmonicSlicer::computePostProcessImpl(const PrintSynSetting &settings, bool useParallel) {
  const auto &domain=i_Settings.HPrintHeight;PrintPathArray surfaceContours;SurfaceFootprint footprint;
  const auto pipeline=[&](){const bool enabled=domain.Enabled&&domain.OptimizeSurfaces;
    if(!computePrintPath(useParallel))return false;PrintPathArray rawContours;if(enabled)rawContours=o_PrintPaths2D;
    if(!zPrintSyn::computePipeline(o_IsoMeshes2D,o_IsoMeshes,nullptr,
      zPrintSyn::SeamMapping::WorldNearest,settings,o_PrintPaths2D,o_PrintPaths3D,
      o_ResultantPrintPath,_Status,useParallel,&o_VertexLayerHeights,&o_VertexLayerNormals,enabled?&surfaceContours:nullptr,&o_FaceTagSources))return false;
    if(enabled){
      // Validate the entire original SDF curve, including the future seam
      // interval. Sampling/feature reconstruction cannot hide an excursion.
      for(size_t k=0;k<surfaceContours.size();++k)surfaceContours[k].insert(surfaceContours[k].end(),rawContours[k].begin(),rawContours[k].end());}
    if(enabled&&!surfaceFootprint(surfaceContours,o_IsoMeshes2D,o_IsoMeshes,o_VertexLayerHeights,domain,useParallel,footprint,_Status)){
      o_ResultantPrintPath=PrintPath();return false;}return true;};
  bool ok=pipeline();
  if(!ok||!domain.Enabled||!domain.OptimizeSurfaces)return ok;
  auto total=_PrintHeightOptimizationReport;
  const int initialSections=int(o_IsoMeshes.size());std::vector<int> candidates{initialSections};
  struct SurfaceState {
    zObjMeshArray Iso,Charts;zObjGraphArray Boundary,Boundary2D;
    std::vector<zPointArray> UV;std::vector<zFloatArray> Heights;std::vector<zVectorArray> Normals;
    std::vector<EdgeWidthIds> EdgeWidths;std::vector<zIntArray> Widths,RegionTags;std::vector<double> Levels;
    PrintPathArray Paths2D,Paths3D;zTsRelax::HeightReport Report;zTsRelax::Report Relaxation;
    bool Valid=false;
  } best;
  const auto retain=[&](){if(best.Valid&&best.Report.WorstViolation<=total.WorstViolation)return;
    best.Iso=o_IsoMeshes;best.Charts=o_IsoMeshes2D;best.Boundary=o_IsoMeshBoundary;best.Boundary2D=o_IsoMeshBoundary2D;
    best.UV=o_IsoMeshUV;best.Heights=o_VertexLayerHeights;best.Normals=o_VertexLayerNormals;best.EdgeWidths=_HarmonicIsoEdgeWidthIds;
    best.Widths=o_IsoMeshBoundaryWidthId;best.RegionTags=o_IsoMeshBoundaryFaceTags;best.Levels=_IsoValues;best.Paths2D=o_PrintPaths2D;best.Paths3D=o_PrintPaths3D;
    best.Report=total;best.Relaxation=_IsoRelaxationReport;best.Valid=true;};
  std::string lastCandidateError;
  for(size_t candidate=0;candidate<candidates.size();++candidate){
    if(candidate){_RequestedSectionCount=candidates[candidate];if(!computeContours(useParallel)||!computeField(useParallel)||!pipeline()){
        lastCandidateError=_Status;zTsUtils::printDebug("[surfaceCount] sections="+std::to_string(candidates[candidate])+" rejected: "+_Status);continue;}
      total.AcceptedIterations+=_PrintHeightOptimizationReport.AcceptedIterations;total.Milliseconds+=_PrintHeightOptimizationReport.Milliseconds;
      total.Factorizations+=_PrintHeightOptimizationReport.Factorizations;total.OuterIterations+=_PrintHeightOptimizationReport.OuterIterations;}
    // One surface-only feedback cycle is allowed for the actual fabrication
    // footprint. No individual target moves; regenerate charts/fields afterward.
    for(int feedback=0;feedback<=1;++feedback){double lo=std::numeric_limits<double>::infinity(),hi=0;int bad=0,points=0;
      for(const auto &layer:o_PrintPaths3D)for(const auto &path:layer)for(float h:path.heights){
        lo=(std::min)(lo,double(h));hi=(std::max)(hi,double(h));++points;
        if(!std::isfinite(h)||h<domain.Minimum-domain.Tolerance||h>domain.Maximum+domain.Tolerance)++bad;}
      lo=(std::min)(lo,footprint.Minimum);hi=(std::max)(hi,footprint.Maximum);bad+=footprint.Violations;
      total.Converged=points>0&&bad==0;total.Minimum=lo;total.Maximum=hi;total.ViolatingVertices=bad;
      total.WorstViolation=(std::max)({0.0,domain.Minimum-lo,hi-domain.Maximum});total.InitialSections=initialSections;total.FinalSections=int(o_IsoMeshes.size());
      if(total.Converged){total.LastRejectedReason.clear();_PrintHeightOptimizationReport=total;return true;}
      if(points==0){_Status="Surface optimization produced no final print targets.";o_ResultantPrintPath=PrintPath();return false;}
      retain();
      if(candidate==0&&feedback==0&&domain.AllowLayerCountChange){
        // Estimate from the measured distribution, without allowing one tight
        // corner or its feature density to determine the whole stack count.
        std::vector<double> layerMeans;
        for(const auto &layer:o_PrintPaths3D){double integral=0,length=0;
          for(const auto &path:layer){zPointArray positions;zIntArray edges;zFnGraph fn(const_cast<zObjGraph&>(path.graph));fn.getVertexPositions(positions);fn.getEdgeData(edges);
            for(size_t e=0;e+1<edges.size();e+=2){const int a=edges[e],b=edges[e+1];const double span=(toEigenVec(positions[b])-toEigenVec(positions[a])).norm();
              integral+=.5*span*(path.heights[a]+path.heights[b]);length+=span;}}
          if(length>0)layerMeans.push_back(integral/length);}
        double typical=.5*(lo+hi);if(!layerMeans.empty()){auto middle=layerMeans.begin()+layerMeans.size()/2;std::nth_element(layerMeans.begin(),middle,layerMeans.end());typical=*middle;}
        const int intervals=initialSections-1;const double estimate=intervals*2*typical/(domain.Minimum+domain.Maximum);
        const int low=(std::max)(2,int(std::floor(estimate))),high=(std::max)(2,int(std::ceil(estimate)));
        for(int count:{low+1,high+1})if(std::find(candidates.begin(),candidates.end(),count)==candidates.end())candidates.push_back(count);
      }
      if(feedback==1)break;
      auto targets=footprint.Targets;
      for(size_t k=0;k<targets.size();++k)for(const auto &path:o_PrintPaths3D[k]){zPointArray points;
        zFnGraph(const_cast<zObjGraph&>(path.graph)).getVertexPositions(points);targets[k].insert(targets[k].end(),points.begin(),points.end());}
      std::string error;if(!optimizeSurfaceHeightDomain(&targets,error,useParallel)){lastCandidateError="Surface feedback failed: "+error;break;}
      total.AcceptedIterations+=_PrintHeightOptimizationReport.AcceptedIterations;total.Milliseconds+=_PrintHeightOptimizationReport.Milliseconds;
      total.Factorizations+=_PrintHeightOptimizationReport.Factorizations;total.OuterIterations+=_PrintHeightOptimizationReport.OuterIterations;
      total.MaxDisplacement+=_PrintHeightOptimizationReport.MaxDisplacement;total.MaxBoundaryError=(std::max)(total.MaxBoundaryError,_PrintHeightOptimizationReport.MaxBoundaryError);
      total.FinalBendingEnergy=_PrintHeightOptimizationReport.FinalBendingEnergy;
      if(!rebuildSurfaceCharts(useParallel)||!pipeline()){lastCandidateError=_Status;break;}
    }
  }
  if(best.Valid){o_IsoMeshes=std::move(best.Iso);o_IsoMeshes2D=std::move(best.Charts);o_IsoMeshBoundary=std::move(best.Boundary);
    o_IsoMeshBoundary2D=std::move(best.Boundary2D);o_IsoMeshUV=std::move(best.UV);o_VertexLayerHeights=std::move(best.Heights);
    o_VertexLayerNormals=std::move(best.Normals);_HarmonicIsoEdgeWidthIds=std::move(best.EdgeWidths);o_IsoMeshBoundaryWidthId=std::move(best.Widths);o_IsoMeshBoundaryFaceTags=std::move(best.RegionTags);_IsoValues=std::move(best.Levels);
    _IsoRelaxationReport=best.Relaxation;computeField(useParallel);o_PrintPaths2D=std::move(best.Paths2D);o_PrintPaths3D=std::move(best.Paths3D);
    best.Report.AcceptedIterations=total.AcceptedIterations;best.Report.OuterIterations=total.OuterIterations;
    best.Report.Factorizations=total.Factorizations;best.Report.Milliseconds=total.Milliseconds;total=best.Report;}
  std::ostringstream message;message<<"Surface height domain not achieved after "<<candidates.size()<<" layer-count candidates: "
    <<total.ViolatingVertices<<" path/continuous-curve supports outside ["<<domain.Minimum<<","<<domain.Maximum<<"] +/-"<<domain.Tolerance
    <<", measured ["<<total.Minimum<<","<<total.Maximum<<"].";
  if(!lastCandidateError.empty())message<<" Last rejected candidate: "<<lastCandidateError;
  total.LastRejectedReason=message.str();_PrintHeightOptimizationReport=total;_Status=message.str();o_ResultantPrintPath=PrintPath();return false;
}

} // namespace zSpace
