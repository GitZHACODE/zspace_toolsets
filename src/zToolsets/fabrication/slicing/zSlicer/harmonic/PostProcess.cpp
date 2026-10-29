#ifndef EIGEN_DONT_PARALLELIZE
#define EIGEN_DONT_PARALLELIZE
#endif

#include <zspace/zToolsets/fabrication/slicing/zSlicer/HarmonicSlicer.h>

#ifndef ZSPACE_INLINE
#define ZSPACE_INLINE ZSPACE_TOOLSETS_INLINE
#endif

#include <sstream>
#include <array>
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
template<class Array> Array tcpSections(const Array &all) {
  // Keep section zero as an immutable support, not a robot target.
  return all.empty()?Array{}:Array(all.begin()+1,all.end());
}
bool measureTCPPaths(const std::vector<zTsRelax::Layer> &layers,PrintPathArray &paths,
    bool parallel,std::string &error) {
  zTsUtils::ScopedSlicingTiming timing("HarmonicPostProcess.measureTCPPaths");
  std::vector<zPointArray> targets(paths.size());std::vector<zVectorArray> normals(paths.size());
  for(size_t k=0;k<paths.size();++k)for(const auto &path:paths[k]){zPointArray p;
    zFnGraph(const_cast<zObjGraph&>(path.graph)).getVertexPositions(p);
    targets[k].insert(targets[k].end(),p.begin(),p.end());
    normals[k].insert(normals[k].end(),path.normals.begin(),path.normals.end());}
  std::vector<zFloatArray> heights;
  if(!zTsRelax::measureTargetHeights(layers,targets,normals,heights,error,parallel))return false;
  for(size_t k=0;k<paths.size();++k){size_t cursor=0;for(auto &path:paths[k]){
    const size_t n=size_t(zFnGraph(path.graph).numVertices());
    path.heights.assign(heights[k].begin()+cursor,heights[k].begin()+cursor+n);cursor+=n;}}
  return true;
}
// Compact hulls guide geometry correction. Physical TCP heights are nonlinear
// under normal interpolation: validate original and sampled curves with actual
// previous-section rays at chart crossings and at most 1 mm chart intervals.
// These diagnostic samples never change fabrication targets or their counts.
bool surfaceFootprint(const PrintPathArray &curves,const zObjMeshArray &charts,
    const zObjMeshArray &sections,const std::vector<zVectorArray> &normals,
    const std::vector<zTsRelax::Layer> &layers,
    const PrintHeightOptimizationSettings &domain,bool parallel,SurfaceFootprint &out,std::string &error){
  zTsUtils::ScopedSlicingTiming footprintTiming("HarmonicPostProcess.surfaceFootprint");
  out=SurfaceFootprint();out.Targets.resize(curves.size());
  std::vector<zPointArray> probes(curves.size());std::vector<zVectorArray> probeNormals(curves.size());
  std::vector<zIntArray> probeTriangles(curves.size());
  std::vector<std::string> errors(curves.size());zTsDetail::ParallelErrors failures(curves.size());
  #pragma omp parallel for schedule(dynamic) if(parallel)
  for(int k=0;k<int(curves.size());++k){try{
    zTsDetail::SurfaceQuery query(const_cast<zObjMesh&>(charts[k]));zPointArray world,chart;
    zFnMesh(const_cast<zObjMesh&>(sections[k])).getVertexPositions(world);zFnMesh(const_cast<zObjMesh&>(charts[k])).getVertexPositions(chart);
    if(size_t(k)>=normals.size()||(!normals[k].empty()&&normals[k].size()!=world.size())){errors[k]="TCP footprint normals do not match the section.";continue;}
    double scale=1;for(const auto &p:chart)scale=(std::max)(scale,toEigenVec(p).cwiseAbs().maxCoeff());
    const double tolerance=(std::max)(1e-7,8*std::numeric_limits<float>::epsilon()*scale);
    std::map<int,std::vector<V>> supports;
    for(const auto &curve:curves[k]){zPointArray p;zIntArray edges;zFnGraph fn(const_cast<zObjGraph&>(curve.graph));fn.getVertexPositions(p);fn.getEdgeData(edges);
      for(size_t e=0;e+1<edges.size();e+=2){const auto &a=p[edges[e]],&b=p[edges[e+1]];std::vector<double> breaks;
        if(!query.chartSegment(a,b,tolerance,breaks)){errors[k]="Continuous height footprint leaves its chart or has overlapping supports.";break;}
        const double step=(std::min)(1.0,.25*domain.Minimum);
        const double length=(toEigenVec(b)-toEigenVec(a)).norm();
        const double pieces=std::ceil(length/step);
        if(!std::isfinite(pieces)||pieces>2000000){errors[k]="TCP height validation sample budget exceeded.";break;}
        const int count=(std::max)(1,int(pieces));
        for(int j=0;j<=count;++j)breaks.push_back(double(j)/count);
        std::sort(breaks.begin(),breaks.end());
        breaks.erase(std::unique(breaks.begin(),breaks.end(),[](double a,double b){return std::abs(a-b)<1e-10;}),breaks.end());
        for(double t:breaks){const V q=(1-t)*toEigenVec(a)+t*toEigenVec(b);zTsDetail::SurfaceQuery::Hit hit;
          if(!query.uniqueChartHit(zPoint(q.x(),q.y(),q.z()),hit,tolerance)){errors[k]="Cannot map continuous height footprint.";break;}
          V w=hit.Weights.cwiseMax(0);w/=w.sum();const auto &ids=query.triangle(hit.Triangle).Ids;
          V p=V::Zero(),n=V::Zero();for(int j=0;j<3;++j){p+=w[j]*toEigenVec(world[ids[j]]);
            if(!normals[k].empty())n+=w[j]*toEigenVec(normals[k][ids[j]]);}
          if(normals[k].empty())n=(toEigenVec(world[ids[1]])-toEigenVec(world[ids[0]])).cross(toEigenVec(world[ids[2]])-toEigenVec(world[ids[0]]));
          if(!n.allFinite()||n.norm()<1e-12){errors[k]="Cannot recover a TCP footprint normal.";break;}n.normalize();
          probes[k].emplace_back(p.x(),p.y(),p.z());probeNormals[k].emplace_back(n.x(),n.y(),n.z());
          probeTriangles[k].push_back(hit.Triangle);supports[hit.Triangle].push_back(w);}
        if(!errors[k].empty())break;}
      if(!errors[k].empty())break;}
    if(!errors[k].empty())continue;
    const auto cross=[](const V &a,const V &b,const V &c){return (b[1]-a[1])*(c[2]-a[2])-(b[2]-a[2])*(c[1]-a[1]);};
    for(auto &entry:supports){auto &points=entry.second;
      // Hull vertices must be actual curve samples. Artificial padding towards
      // triangle corners can be unsupported even when the TCP curve is valid.
      std::sort(points.begin(),points.end(),[](const V &a,const V &b){return a[1]<b[1]||(a[1]==b[1]&&a[2]<b[2]);});
      points.erase(std::unique(points.begin(),points.end(),[](const V &a,const V &b){return (a-b).norm()<1e-10;}),points.end());
      std::vector<V> hull;for(const auto &p:points){while(hull.size()>1&&cross(hull[hull.size()-2],hull.back(),p)<=1e-12)hull.pop_back();hull.push_back(p);}
      const size_t lower=hull.size();if(points.size()>1){for(int i=int(points.size())-2;i>=0;--i){while(hull.size()>lower&&cross(hull[hull.size()-2],hull.back(),points[i])<=1e-12)hull.pop_back();hull.push_back(points[i]);}hull.pop_back();}
      const auto &ids=query.triangle(entry.first).Ids;for(const auto &w:hull){V p=V::Zero();for(int j=0;j<3;++j)p+=w[j]*toEigenVec(world[ids[j]]);out.Targets[k].emplace_back(p.x(),p.y(),p.z());}}
  }catch(...){failures.capture(k);}}
  failures.rethrow();for(const auto &message:errors)if(!message.empty()){error=message;return false;}
  std::vector<zFloatArray> measured;
  if(!zTsRelax::measureTargetHeights(layers,probes,probeNormals,measured,error,parallel))return false;
  for(size_t k=0;k<curves.size();++k){std::map<int,std::array<std::pair<double,size_t>,2>> worst;
    for(size_t i=0;i<measured[k].size();++i){const double h=measured[k][i];
      out.Minimum=(std::min)(out.Minimum,h);out.Maximum=(std::max)(out.Maximum,h);
      const double v=(std::max)({0.0,domain.Minimum-h,h-domain.Maximum});
      if(v>domain.Tolerance){++out.Violations;auto &entry=worst[probeTriangles[k][i]][h<domain.Minimum?0:1];if(v>entry.first)entry={v,i};}}
    // Actual ray height can peak inside a triangle. Include its worst observed
    // violation instead of assuming hull endpoint heights bound the interior.
    for(const auto &entry:worst)for(const auto &bound:entry.second)if(bound.first>0)out.Targets[k].push_back(probes[k][bound.second]);}
  if(!std::isfinite(out.Minimum)){error="No TCP curve validation samples.";return false;}return true;
}
}

//--------------------------
//---- PUBLIC METHODS: zSlicer
//--------------------------

ZSPACE_INLINE bool HarmonicSlicer::computeSlicer(bool useParallel) {
  zTsUtils::ScopedSlicingTiming methodTiming("HarmonicSlicer.computeSlicer");
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
  zTsUtils::ScopedSlicingTiming methodTiming("HarmonicSlicer.computePostProcess");
  const bool success = computePostProcessImpl(settings, useParallel);
  if (success) updatePrintStatistics(true);
  else resetPrintStatistics();
  if (success) zTsUtils::printSuccess("[HarmonicSlicer::computePostProcess] Success");
  else zTsUtils::printError("[HarmonicSlicer::computePostProcess] Failed");
  return success;
}


bool HarmonicSlicer::rebuildSurfaceCharts(bool parallel){
  zTsUtils::ScopedSlicingTiming methodTiming("HarmonicSlicer.rebuildSurfaceCharts");
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
    const auto charts=tcpSections(o_IsoMeshes2D),sections=tcpSections(o_IsoMeshes);
    const auto heights=tcpSections(o_VertexLayerHeights);const auto normals=tcpSections(o_VertexLayerNormals);
    const auto sources=tcpSections(o_FaceTagSources);
    if(!zPrintSyn::computePipeline(charts,sections,nullptr,
      zPrintSyn::SeamMapping::WorldNearest,settings,o_PrintPaths2D,o_PrintPaths3D,
      o_ResultantPrintPath,_Status,useParallel,&heights,&normals,enabled?&surfaceContours:nullptr,&sources))return false;
    const auto layerGeometry=i_Settings._HIsoRelaxIterations>0?relaxationLayers():std::vector<zTsRelax::Layer>{};
    if(i_Settings._HIsoRelaxIterations>0){
      if(!measureTCPPaths(layerGeometry,o_PrintPaths3D,useParallel,_Status)||
          !zPrintSyn::linkPrintPaths(o_PrintPaths3D,o_ResultantPrintPath,_Status)){
        o_ResultantPrintPath=PrintPath();return false;}}
    if(enabled){
      // Validate the entire original SDF curve, including the future seam
      // interval. Sampling/feature reconstruction cannot hide an excursion.
      for(size_t k=0;k<surfaceContours.size();++k)surfaceContours[k].insert(surfaceContours[k].end(),rawContours[k].begin(),rawContours[k].end());}
    if(enabled&&!surfaceFootprint(surfaceContours,charts,sections,normals,layerGeometry,domain,useParallel,footprint,_Status)){
      o_ResultantPrintPath=PrintPath();return false;}return true;};
  bool ok=pipeline();
  if(!ok||!domain.Enabled||!domain.OptimizeSurfaces)return ok;
  auto total=_PrintHeightOptimizationReport;
  // Iterations is one correction budget for the complete fabrication pipeline,
  // not a fresh budget for every feedback pass and layer-count candidate.
  const int iterationLimit=domain.Iterations;
  const auto remainingIterations=[&](){return (std::max)(0,iterationLimit-total.OuterIterations);};
  const auto accumulateWork=[&](const zTsRelax::HeightReport &r){
    if(total.OuterIterations==0&&total.Factorizations==0){total.InitialMinimum=r.InitialMinimum;
      total.InitialMaximum=r.InitialMaximum;total.InitialBendingEnergy=r.InitialBendingEnergy;}
    total.AcceptedIterations+=r.AcceptedIterations;total.Milliseconds+=r.Milliseconds;
    total.Factorizations+=r.Factorizations;total.OuterIterations+=r.OuterIterations;
    total.InnerIterations+=r.InnerIterations;total.LineSearchTrials+=r.LineSearchTrials;
    total.DegreesOfFreedom=(std::max)(total.DegreesOfFreedom,r.DegreesOfFreedom);
    total.InitialSliverVertices=r.InitialSliverVertices;total.FixedVertices=r.FixedVertices;
    total.UncontrollableSamples=r.UncontrollableSamples;
    total.UnavailableAuxiliaryGaps=r.UnavailableAuxiliaryGaps;
    total.MinimumDirectedGap=r.MinimumDirectedGap;total.MaximumDirectedGap=r.MaximumDirectedGap;
  };
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
  std::string lastCandidateError;int attemptedCandidates=1;
  for(size_t candidate=0;candidate<candidates.size();++candidate){
    if(candidate){
      if(remainingIterations()==0)break;
      ++attemptedCandidates;
      _RequestedSectionCount=candidates[candidate];const bool contoursReady=computeContours(useParallel);
      accumulateWork(_PrintHeightOptimizationReport);
      if(!contoursReady||!computeField(useParallel)||!pipeline()){
        lastCandidateError=_Status;zTsUtils::printDebug("[surfaceCount] sections="+std::to_string(candidates[candidate])+" rejected: "+_Status);continue;}}
    double candidateMovement=0;
    const double movementLimit=domain.MaxDisplacement>0?domain.MaxDisplacement:2*domain.Maximum;
    // At most two geometry corrections for the actual fabrication footprint.
    // The second can repair a shifted material footprint after chart rebuild;
    // both share the same outer-iteration budget. No target moves individually.
    for(int feedback=0;feedback<=2;++feedback){double lo=std::numeric_limits<double>::infinity(),hi=0;int bad=0,points=0;
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
      if(feedback==2)break;
      if(remainingIterations()==0){lastCandidateError="Surface correction iteration budget exhausted.";break;}
      const double remainingMovement=movementLimit-candidateMovement;
      if(remainingMovement<=1e-9){lastCandidateError="Surface correction displacement budget exhausted.";break;}
      auto targets=footprint.Targets;
      for(size_t k=0;k<targets.size();++k)for(const auto &path:o_PrintPaths3D[k]){zPointArray points;
        zFnGraph(const_cast<zObjGraph&>(path.graph)).getVertexPositions(points);targets[k].insert(targets[k].end(),points.begin(),points.end());}
      // A constrained solver can return valid geometry without moving it.
      // Rebuilding ARAP charts/SDFs/paths in that case cannot improve the same
      // footprint and used to repeat the most expensive downstream stages.
      std::vector<zPointArray> before(o_IsoMeshes.size());
      for(size_t k=1;k+1<o_IsoMeshes.size();++k)zFnMesh(o_IsoMeshes[k]).getVertexPositions(before[k]);
      std::string error;const bool corrected=optimizeSurfaceHeightDomain(&targets,error,useParallel,remainingIterations(),remainingMovement);
      accumulateWork(_PrintHeightOptimizationReport);
      if(!corrected){lastCandidateError="Surface feedback failed: "+error;break;}
      candidateMovement+=_PrintHeightOptimizationReport.MaxDisplacement;
      total.MaxDisplacement+=_PrintHeightOptimizationReport.MaxDisplacement;total.MaxBoundaryError=(std::max)(total.MaxBoundaryError,_PrintHeightOptimizationReport.MaxBoundaryError);
      total.FinalBendingEnergy=_PrintHeightOptimizationReport.FinalBendingEnergy;
      bool moved=false;
      for(size_t k=1;k+1<o_IsoMeshes.size()&&!moved;++k){zPointArray after;zFnMesh(o_IsoMeshes[k]).getVertexPositions(after);
        if(after.size()!=before[k].size()){moved=true;break;}
        for(size_t i=0;i<after.size();++i)if(zTsUtils::vecDistSq(after[i],before[k][i])!=0){moved=true;break;}}
      if(!moved){lastCandidateError="Surface feedback is constrained: no geometry changed.";
        zTsUtils::printDebug("[surfaceHeightOptimizer] unchanged surfaces; skipped chart/field/path rebuild.");break;}
      if(!rebuildSurfaceCharts(useParallel)||!pipeline()){lastCandidateError=_Status;break;}
    }
  }
  if(best.Valid){o_IsoMeshes=std::move(best.Iso);o_IsoMeshes2D=std::move(best.Charts);o_IsoMeshBoundary=std::move(best.Boundary);
    o_IsoMeshBoundary2D=std::move(best.Boundary2D);o_IsoMeshUV=std::move(best.UV);o_VertexLayerHeights=std::move(best.Heights);
    o_VertexLayerNormals=std::move(best.Normals);_HarmonicIsoEdgeWidthIds=std::move(best.EdgeWidths);o_IsoMeshBoundaryWidthId=std::move(best.Widths);o_IsoMeshBoundaryFaceTags=std::move(best.RegionTags);_IsoValues=std::move(best.Levels);
    _IsoRelaxationReport=best.Relaxation;computeField(useParallel);o_PrintPaths2D=std::move(best.Paths2D);o_PrintPaths3D=std::move(best.Paths3D);
    best.Report.AcceptedIterations=total.AcceptedIterations;best.Report.OuterIterations=total.OuterIterations;
    best.Report.Factorizations=total.Factorizations;best.Report.Milliseconds=total.Milliseconds;
    best.Report.InnerIterations=total.InnerIterations;best.Report.LineSearchTrials=total.LineSearchTrials;
    best.Report.DegreesOfFreedom=total.DegreesOfFreedom;total=best.Report;}
  std::ostringstream message;message<<"Surface height domain not achieved after "<<attemptedCandidates<<" layer-count candidates and "<<total.OuterIterations<<"/"<<iterationLimit<<" correction iterations: "
    <<total.ViolatingVertices<<" TCP/curve validation samples outside ["<<domain.Minimum<<","<<domain.Maximum<<"] +/-"<<domain.Tolerance
    <<", measured ["<<total.Minimum<<","<<total.Maximum<<"].";
  if(!lastCandidateError.empty())message<<" Last rejected candidate: "<<lastCandidateError;
  total.LastRejectedReason=message.str();_PrintHeightOptimizationReport=total;_Status=message.str();o_ResultantPrintPath=PrintPath();return false;
}

} // namespace zSpace
