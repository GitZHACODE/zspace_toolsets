#include <zspace/zToolsets/fabrication/slicing/printsyn/zPrintSyn.h>
#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>
#include <stdexcept>
#include <zspace/zToolsets/fabrication/slicing/utils/Parallel.h>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace zSpace {
namespace {
struct RegionRun {
  int tag=0;
  zIntArray ids;
  std::vector<double> arc;
  std::vector<double> worldArc;
  std::vector<int> anchors;
  double length() const { return arc.back(); }
};
struct RegionPath {
  PrintPath input;
  zPointArray points;
  zPointArray world;
  std::vector<RegionRun> runs;
  bool closed=false;
};


constexpr double AnchorClearance = 0.01;

void refreshFeatures(PrintPath &path,float threshold) {
  path.featureFlags.assign(zFnGraph(path.graph).numVertices(),0);
  for(int id:zPrintSyn::cornerVertexIds(path,threshold))path.featureFlags[id]=PrintPath::AngleFeature;
}

// Collapse one contiguous span to an existing exact representative. Keeping
// each surviving original edge's tag preserves ownership on both sides of a
// merged transition, including a span crossing the closed seam.
bool mergeAnchorSpan(PrintPath &path,zPointArray &world,const zIntArray &span,int keep,bool closed) {
  zPointArray points;zIntArray edges;
  zFnGraph(path.graph).getVertexPositions(points);zFnGraph(path.graph).getEdgeData(edges);
  const int n=int(points.size()),root=*std::min_element(span.begin(),span.end());
  std::vector<int> parent(n);std::iota(parent.begin(),parent.end(),0);
  for(int v:span)parent[v]=root;
  std::map<int,int> groups;std::vector<int> indices(n),representatives;
  for(int i=0;i<n;++i){auto entry=groups.emplace(parent[i],int(groups.size()));
    if(entry.second)representatives.push_back(parent[i]==root?keep:i);indices[i]=entry.first->second;
  }
  PrintPath out;zPointArray result,mapped;zIntArray connects;
  out.featureFlags.assign(representatives.size(),0);
  for(int representative:representatives) {
    result.push_back(points[representative]);out.widths.push_back(path.widths[representative]);
    if(!world.empty())mapped.push_back(world[representative]);
    if(!path.heights.empty())out.heights.push_back(path.heights[representative]);
    if(!path.normals.empty())out.normals.push_back(path.normals[representative]);
    if(!path.trackIds.empty())out.trackIds.push_back(path.trackIds[representative]);
  }
  out.faceTags.assign(result.size(),0);
  for(size_t e=0;e<edges.size()/2;++e){const int a=indices[edges[2*e]],b=indices[edges[2*e+1]];if(a==b)continue;
    connects.push_back(a);connects.push_back(b);const int tag=path.edgeFaceTags[e];out.edgeFaceTags.push_back(tag);
    for(int v:{a,b})if(tag && (!out.faceTags[v] || tag<out.faceTags[v]))out.faceTags[v]=tag;
  }
  if(result.size()<(closed?3u:2u) || connects.size()<(closed?6u:2u))return false;
  zFnGraph(out.graph).create(result,connects);path=std::move(out);world=std::move(mapped);return true;
}

// Normalize protected anchors before computing any fixed-count budget. Dense
// ordinary contour vertices remain intact. A short span may lose a categorical
// run; retaining both endpoints would contradict the requested 0.01 mm weld.
bool weldProjectedAnchors(PrintPath &path,float threshold,zPointArray &world,
    zObjMesh *chart,zObjMesh *section) {
  zTsUtils::ScopedSlicingTiming methodTiming("FaceTagSampling.weldProjectedAnchors", true);
  if(chart&&section){zObjGraph mapped=path.graph;
    if(!zTsUtils::remapGraph(mapped,*chart,*section))return false;
    zFnGraph(mapped).getVertexPositions(world);
  }
  for(;;) {
    zPointArray points;zFnGraph(path.graph).getVertexPositions(points);
    std::vector<zTsUtils::GraphPolyline> curves;
    if(!zTsUtils::graphPolylines(path.graph,curves)||curves.size()!=1)return false;
    const bool closed=curves[0].Closed;const int n=int(points.size()),e=int(path.edgeFaceTags.size());
    if(e!=n-(closed?0:1)||(!world.empty()&&world.size()!=points.size()))return false;
    refreshFeatures(path,threshold);
    zIntArray span;int keep=-1;
    // Retain the exact-world weld for chart vertices that map to one point,
    // even when an ordinary vertex separates the protected anchors.
    for(int a=0;a<e;++a){const int b=(a+1)%n;
      if(zTsUtils::vecDistSq(points[a],points[b])<=1e-18 ||
          (!world.empty()&&zTsUtils::vecDistSq(world[a],world[b])<=1e-18)){
        span={a,b};keep=path.featureFlags[b]&&!path.featureFlags[a]?b:a;break;
      }
    }
    if(span.empty()) {
      zIntArray anchors;
      for(int v=0;v<n;++v)if(path.featureFlags[v] || (!closed&&(v==0||v==n-1)) ||
          (closed?path.edgeFaceTags[(v+e-1)%e]!=path.edgeFaceTags[v]:
            v>0&&v<e&&path.edgeFaceTags[v-1]!=path.edgeFaceTags[v]))anchors.push_back(v);
      const int pairs=int(anchors.size())-(closed?0:1);
      for(int i=0;anchors.size()>1&&i<pairs;++i){const int a=anchors[i],b=anchors[(i+1)%anchors.size()];
        double chartLength=0,worldLength=0;zIntArray candidate{a};
        for(int v=a;v!=b;){const int next=(v+1)%n;
          chartLength+=std::sqrt(zTsUtils::vecDistSq(points[v],points[next]));
          if(!world.empty())worldLength+=std::sqrt(zTsUtils::vecDistSq(world[v],world[next]));
          candidate.push_back(next);v=next;
        }
        // The intervening arc check prevents welding opposite sides of a long
        // narrow detour whose endpoints happen to be spatially close.
        if(chartLength<AnchorClearance || (!world.empty()&&worldLength<AnchorClearance)){
          span=std::move(candidate);
          const bool endA=!closed&&(a==0||a==n-1),endB=!closed&&(b==0||b==n-1);
          keep=endB&&!endA?b:endA?a:path.featureFlags[b]&&!path.featureFlags[a]?b:a;
          break;
        }
      }
    }
    if(span.empty())return true;
    if(!mergeAnchorSpan(path,world,span,keep,closed))return false;
  }
}

// Integer intervals with a lower bound of one per protected span. Remaining
// rounding units go to largest remainders, with stable source-order ties.
std::vector<int> allocateIntervals(const std::vector<double> &lengths, int budget,
    const std::vector<int> &capacities) {
  zTsUtils::ScopedSlicingTiming methodTiming("FaceTagSampling.allocateIntervals", true);
  const int n=int(lengths.size());
  if(budget<n) throw std::runtime_error("FaceTag budget cannot preserve anchors");
  if(std::any_of(capacities.begin(),capacities.end(),[](int n){return n<1;}) ||
      std::accumulate(capacities.begin(),capacities.end(),int64_t(0))<budget)
    throw std::runtime_error("FaceTag fixed count cannot fit with 0.01 mm target clearance");
  std::vector<int> result(n,0); std::vector<double> quota(n);
  int remaining=budget; double length=std::accumulate(lengths.begin(),lengths.end(),0.0);
  bool changed=true;
  while(changed) {
    changed=false;
    // Saturate upper bounds first. Doing a lower bound first could strand
    // budget that should move to that span after a different span saturates.
    for(int i=0;i<n;++i) if(!result[i] && remaining*lengths[i]/length>capacities[i]) {
      result[i]=capacities[i];
      remaining-=result[i]; length-=lengths[i]; changed=true;
    }
    if(changed)continue;
    for(int i=0;i<n;++i) if(!result[i] && remaining*lengths[i]/length<1) {
      result[i]=1;--remaining;length-=lengths[i];changed=true;
    }
  }
  std::vector<int> order;
  for(int i=0;i<n;++i) if(!result[i]) {
    quota[i]=remaining*lengths[i]/length; order.push_back(i);
  }
  int used=0;
  for(int i:order) { result[i]=int(std::floor(quota[i])); used+=result[i]; }
  std::stable_sort(order.begin(),order.end(),[&](int a,int b){
    return quota[a]-std::floor(quota[a]) > quota[b]-std::floor(quota[b]);
  });
  for(int i=0;i<remaining-used;++i) ++result[order.at(i)];
  return result;
}
}

bool zPrintSyn::resampleFaceTaggedPaths(PrintPathArray &paths,
    const std::vector<zTsField::FaceTagSource> &sources,
    const PrintSynSetting &settings, std::string &status,
    const zObjMeshArray *charts2D, const zObjMeshArray *sections3D) {
  return resampleFaceTaggedPaths(paths,sources,settings,status,charts2D,sections3D,false);
}

bool zPrintSyn::resampleFaceTaggedPaths(PrintPathArray &paths,
    const std::vector<zTsField::FaceTagSource> &sources,
    const PrintSynSetting &settings, std::string &status,
    const zObjMeshArray *charts2D, const zObjMeshArray *sections3D, bool useParallel) {
  zTsUtils::ScopedSlicingTiming methodTiming("PrintSyn.resampleFaceTaggedPaths");
  try {
    if(paths.empty() || paths.size()>sources.size() ||
        !std::isfinite(settings.ResampleLength) || settings.ResampleLength<zTsField::MinPhysicalParameter ||
        !std::isfinite(settings.FeatureAngleTol) || settings.FeatureAngleTol<0 || settings.FeatureAngleTol>180)
      throw std::runtime_error("Invalid FaceTag layers or sampling distance");
    std::vector<std::vector<RegionPath>> prepared(paths.size());
    std::vector<std::map<int,double>> lengths(paths.size());
    std::vector<std::map<int,int>> pieces(paths.size()), spans(paths.size());
    // Resolve source geometry and any lazy Core reads before workers share it.
    std::vector<zTsField::FaceTagQuery> queries;queries.reserve(paths.size());
    for(size_t l=0;l<paths.size();++l) {
      if(sources[l].edgeTags.empty()) throw std::runtime_error("Missing FaceTag boundary on layer "+std::to_string(l));
      queries.emplace_back(sources[l]);
      if(charts2D&&sections3D){
        zPointArray points;zIntArray connects,counts;
        zFnMesh chart(const_cast<zObjMesh&>(charts2D->at(l))),section(const_cast<zObjMesh&>(sections3D->at(l)));
        chart.getVertexPositions(points);chart.getPolygonData(connects,counts);
        section.getVertexPositions(points);section.getPolygonData(connects,counts);
      }
    }
    const auto prepareLayer=[&](size_t l){
      const auto &query=queries[l];
      for(auto path:paths[l]) {
        zTsUtils::ScopedSlicingTiming prepareTiming("FaceTagSampling.preparePath", true);
        if(!path.featureFlags.empty()&&path.featureFlags.size()!=size_t(zFnGraph(path.graph).numVertices()))
          throw std::runtime_error("Invalid FaceTag feature metadata");
        refreshFeatures(path,settings.FeatureAngleTol);
        if(!zTsField::labelFaceTagSegments(path,query)) throw std::runtime_error("Cannot label FaceTag contour");
        RegionPath p; p.input=std::move(path);
        if(!weldProjectedAnchors(p.input,settings.FeatureAngleTol,p.world,
            charts2D&&sections3D?&const_cast<zObjMesh&>(charts2D->at(l)):nullptr,
            charts2D&&sections3D?&const_cast<zObjMesh&>(sections3D->at(l)):nullptr))
          throw std::runtime_error("Cannot merge FaceTag anchors within 0.01 mm without collapsing the path");
        zFnGraph(p.input.graph).getVertexPositions(p.points);
        std::vector<zTsUtils::GraphPolyline> curves;
        if(!zTsUtils::graphPolylines(p.input.graph,curves) || curves.size()!=1)
          throw std::runtime_error("Disconnected FaceTag path");
        p.closed=curves.front().Closed;
        const int n=int(p.points.size()), e=int(p.input.edgeFaceTags.size());
        if(e!=n-(p.closed?0:1)) throw std::runtime_error("FaceTag connectivity mismatch");
        auto corners=cornerVertexIds(p.input,settings.FeatureAngleTol);
        int start=0;
        if(p.closed) for(int i=0;i<e;++i) if(p.input.edgeFaceTags[i]!=p.input.edgeFaceTags[(i+e-1)%e]) { start=i; break; }
        for(int k=0;k<e;++k) {
          const int a=(start+k)%n,b=(a+1)%n,tag=p.input.edgeFaceTags[a];
          if(p.runs.empty() || p.runs.back().tag!=tag) {
            RegionRun r;r.tag=tag;r.ids.push_back(a);r.arc.push_back(0);r.worldArc.push_back(0);r.anchors.push_back(0);
            p.runs.push_back(std::move(r));
          }
          auto &r=p.runs.back();
          const double d=std::sqrt(zTsUtils::vecDistSq(p.points[a],p.points[b]));
          if(!std::isfinite(d) || d<=0) throw std::runtime_error("Degenerate FaceTag edge");
          r.ids.push_back(b);r.arc.push_back(r.arc.back()+d);
          r.worldArc.push_back(r.worldArc.back()+(p.world.empty()?d:
            std::sqrt(zTsUtils::vecDistSq(p.world[a],p.world[b]))));
          if(std::find(corners.begin(),corners.end(),b)!=corners.end()) r.anchors.push_back(int(r.ids.size())-1);
        }
        for(auto &r:p.runs) {
          const int end=int(r.ids.size())-1;
          if(r.anchors.back()!=end) r.anchors.push_back(end);
          if(r.tag) { lengths[l][r.tag]+=r.length(); ++pieces[l][r.tag]; spans[l][r.tag]+=int(r.anchors.size())-1; }
        }
        prepared[l].push_back(std::move(p));
      }
    };
    bool parallel=false;
#ifdef _OPENMP
    parallel=useParallel&&paths.size()>1&&!omp_in_parallel();
#endif
    if(parallel){
      zTsDetail::ParallelErrors failures(paths.size());
#pragma omp parallel for schedule(dynamic)
      for(int l=0;l<int(paths.size());++l){try{prepareLayer(size_t(l));}catch(...){failures.capture(size_t(l));}}
      failures.rethrow();
    }else for(size_t l=0;l<paths.size();++l)prepareLayer(l);
    // Only this calling thread writes the global budget map. The reduction's
    // order is unchanged, so floating-point means and rounded counts are exact.
    std::map<int,int> budgets;
    for(size_t l=0;l<paths.size();++l)for(const auto &entry:lengths[l])budgets.emplace(entry.first,0);
    if(budgets.empty()) throw std::runtime_error("No outer FaceTag segments were found");
    for(auto &entry:budgets) {
      zTsUtils::ScopedSlicingTiming budgetTiming("FaceTagSampling.computeTagBudget", true);
      const int tag=entry.first;double mean=0;int maxPieces=0,minimum=0;
      for(size_t l=0;l<paths.size();++l) {
        if(!lengths[l].count(tag)) throw std::runtime_error("FaceTag "+std::to_string(tag)+" missing on layer "+std::to_string(l));
        mean+=lengths[l][tag]/paths.size();
        maxPieces=(std::max)(maxPieces,pieces[l][tag]);
        minimum=(std::max)(minimum,spans[l][tag]+pieces[l][tag]);
      }
      const double intervals=std::round(mean/settings.ResampleLength);
      if(intervals>2000000) throw std::runtime_error("FaceTag sample budget exceeds 2000000 points");
      entry.second=(std::max)(minimum,int(intervals)+maxPieces);
    }
    PrintPathArray result(paths.size());
    for(size_t l=0;l<paths.size();++l) {
      zTsUtils::ScopedSlicingTiming samplingTiming("FaceTagSampling.sampleLayer", true);
      std::map<int,std::vector<double>> spanLengths;
      std::map<int,std::vector<int>> spanCapacities;
      for(const auto &p:prepared[l]) for(const auto &r:p.runs) if(r.tag)
        for(size_t a=1;a<r.anchors.size();++a) {
          spanLengths[r.tag].push_back(r.arc[r.anchors[a]]-r.arc[r.anchors[a-1]]);
          const double physical=(std::min)(spanLengths[r.tag].back(),
            r.worldArc[r.anchors[a]]-r.worldArc[r.anchors[a-1]]);
          spanCapacities[r.tag].push_back(int((std::min)(2000000.0,std::floor(physical/AnchorClearance))));
        }
      std::map<int,std::vector<int>> counts;std::map<int,size_t> cursors;
      for(const auto &b:budgets) counts[b.first]=allocateIntervals(spanLengths[b.first],b.second-pieces[l][b.first],spanCapacities[b.first]);
      for(const auto &p:prepared[l]) {
        PrintPath out;zPointArray points;zIntArray edges;
        const auto append=[&](const RegionRun &r,double arc) {
          auto it=std::upper_bound(r.arc.begin(),r.arc.end(),arc);
          const int k=(std::min)(int(r.ids.size())-2,(std::max)(0,int(it-r.arc.begin())-1));
          const int a=r.ids[k],b=r.ids[k+1];
          const double t=std::clamp((arc-r.arc[k])/(r.arc[k+1]-r.arc[k]),0.0,1.0);
          auto pa=p.points[a],pb=p.points[b];
          const auto point=t==0?pa:t==1?pb:pa+(pb-pa)*t;
          const auto value=[&](float a,float b){return t==0?a:t==1?b:float(a+(b-a)*t);};
          if(!points.empty() && zTsUtils::vecDistSq(points.back(),point)<1e-18) return int(points.size())-1;
          if(p.closed && arc==r.length() && &r==&p.runs.back() && !points.empty() &&
              zTsUtils::vecDistSq(points.front(),point)<1e-18) return 0;
          const int id=int(points.size());points.push_back(point);out.faceTags.push_back(0);
          out.featureFlags.push_back(0);
          out.widths.push_back(value(p.input.widths[a],p.input.widths[b]));
          if(!p.input.heights.empty()) out.heights.push_back(value(p.input.heights[a],p.input.heights[b]));
          if(!p.input.normals.empty()) { auto na=p.input.normals[a],nb=p.input.normals[b]; auto normal=t==0?na:t==1?nb:na+(nb-na)*t; if(normal.length()>0)normal.normalize();out.normals.push_back(normal); }
          if(!p.input.trackIds.empty()) out.trackIds.push_back((p.input.trackIds[a]+(p.input.trackIds[b]-p.input.trackIds[a])*t)>=.5?1.f:0.f);
          return id;
        };
        for(const auto &r:p.runs) {
          PrintPath::FaceTagSegment segment;segment.tag=r.tag;
          int previous=append(r,0);segment.vertices.push_back(previous);
          for(size_t a=1;a<r.anchors.size();++a) {
            const double lo=r.arc[r.anchors[a-1]],hi=r.arc[r.anchors[a]];
            const double physical=(std::min)(hi-lo,r.worldArc[r.anchors[a]]-r.worldArc[r.anchors[a-1]]);
            const int capacity=(std::max)(1,int((std::min)(2000000.0,std::floor(physical/AnchorClearance))));
            const int count=r.tag ? counts[r.tag].at(cursors[r.tag]++) :
                (std::min)(capacity,(std::max)(1,int(std::ceil((hi-lo)/settings.ResampleLength))));
            for(int j=1;j<=count;++j) {
              const double sampleArc=j==count?hi:lo+(hi-lo)*j/count;
              const int next=append(r,sampleArc);
              if(next==previous) {
                if(!r.tag) continue; // Float rounding can coalesce a final regular sample with its endpoint.
                throw std::runtime_error("Collapsed FaceTag sampling interval");
              }
              edges.push_back(previous);edges.push_back(next);out.edgeFaceTags.push_back(r.tag);
              previous=next;segment.vertices.push_back(next);
            }
          }
          if(r.tag) {
            for(int id:segment.vertices) if(!out.faceTags[id] || r.tag<out.faceTags[id]) out.faceTags[id]=r.tag;
            out.faceTagSegments.push_back(std::move(segment));
          }
        }
        zFnGraph(out.graph).create(points,edges);
        {
          zTsUtils::ScopedSlicingTiming clearanceTiming("FaceTagSampling.validateClearance",true);
          zPointArray world;
          if(charts2D&&sections3D){zObjGraph mapped=out.graph;
            if(!zTsUtils::remapGraph(mapped,const_cast<zObjMesh&>(charts2D->at(l)),const_cast<zObjMesh&>(sections3D->at(l))))
              throw std::runtime_error("Cannot validate sampled FaceTag target clearance");
            zFnGraph(mapped).getVertexPositions(world);
          }
          for(size_t e=0;e<edges.size();e+=2)if(
              zTsUtils::vecDistSq(points[edges[e]],points[edges[e+1]])<AnchorClearance*AnchorClearance ||
              (!world.empty()&&zTsUtils::vecDistSq(world[edges[e]],world[edges[e+1]])<AnchorClearance*AnchorClearance))
            throw std::runtime_error("FaceTag fixed count cannot preserve 0.01 mm target clearance on layer "+std::to_string(l));
        }
        refreshFeatures(out,settings.FeatureAngleTol);
        result[l].push_back(std::move(out));
      }
    }
    // Removing a short categorical run can change the number of shared tag-tag
    // endpoints. Under canonical minimum-incident-tag ownership it may then be
    // impossible to keep both segment budgets and point-label totals equal.
    // Report that topology change instead of silently relabelling ownership.
    std::map<int,int> expectedLabels;
    for(size_t l=0;l<result.size();++l){std::map<int,int> labels;
      for(const auto &path:result[l])for(int tag:path.faceTags)if(tag)++labels[tag];
      if(l==0)expectedLabels=labels;
      else if(labels!=expectedLabels)throw std::runtime_error(
        "FaceTag equal segment and point-label counts are incompatible after 0.01 mm anchor merging on layer "+std::to_string(l));
    }
    paths=std::move(result);status="FaceTag sampling (2D, fixed total across pieces):";
    for(const auto &b:budgets) status+=" "+std::to_string(b.first)+"="+std::to_string(b.second);
    return true;
  } catch(const std::exception &error) { status=std::string("failed: ")+error.what();return false; }
}
}
