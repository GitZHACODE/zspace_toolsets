#include <zspace/zToolsets/fabrication/slicing/printsyn/zPrintSyn.h>
#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>
#include <stdexcept>

namespace zSpace {
namespace {
struct RegionRun {
  int tag=0;
  zIntArray ids;
  std::vector<double> arc;
  std::vector<int> anchors;
  double length() const { return arc.back(); }
};
struct RegionPath {
  PrintPath input;
  zPointArray points;
  std::vector<RegionRun> runs;
  bool closed=false;
};


// Coalesce only consecutive anchors that become exactly coincident in the
// printable section. Count budgets are computed afterward, so no layer loses
// targets. This also retains the union of coincident feature identities.
bool weldProjectedAnchors(PrintPath &path,zObjMesh &chart,zObjMesh &section) {
  zObjGraph mapped=path.graph;
  if(!zTsUtils::remapGraph(mapped,chart,section))return false;
  zPointArray world,points;zIntArray edges;
  zFnGraph(mapped).getVertexPositions(world);
  zFnGraph(path.graph).getVertexPositions(points);zFnGraph(path.graph).getEdgeData(edges);
  const int n=int(points.size());if(world.size()!=points.size())return false;
  std::vector<int> parent(n);std::iota(parent.begin(),parent.end(),0);
  const auto find=[&](int v){while(parent[v]!=v){parent[v]=parent[parent[v]];v=parent[v];}return v;};
  bool changed=false;
  for(size_t e=0;e<edges.size();e+=2)if(zTsUtils::vecDistSq(world[edges[e]],world[edges[e+1]])<=1e-18) {
    const int a=find(edges[e]),b=find(edges[e+1]);if(a!=b){parent[(std::max)(a,b)]=(std::min)(a,b);changed=true;}
  }
  if(!changed)return true;
  std::map<int,int> groups;std::vector<int> indices(n),representatives;
  for(int i=0;i<n;++i){const int root=find(i);auto entry=groups.emplace(root,int(groups.size()));
    if(entry.second)representatives.push_back(i);indices[i]=entry.first->second;
  }
  PrintPath out;zPointArray result;zIntArray connects;
  for(int representative:representatives) {
    result.push_back(points[representative]);out.widths.push_back(path.widths[representative]);
    if(!path.heights.empty())out.heights.push_back(path.heights[representative]);
    if(!path.normals.empty())out.normals.push_back(path.normals[representative]);
    if(!path.trackIds.empty())out.trackIds.push_back(path.trackIds[representative]);
  }
  out.faceTags.assign(result.size(),0);
  for(size_t e=0;e<edges.size()/2;++e){const int a=indices[edges[2*e]],b=indices[edges[2*e+1]];if(a==b)continue;
    connects.push_back(a);connects.push_back(b);const int tag=path.edgeFaceTags[e];out.edgeFaceTags.push_back(tag);
    for(int v:{a,b})if(tag && (!out.faceTags[v] || tag<out.faceTags[v]))out.faceTags[v]=tag;
  }
  if(result.size()<2 || connects.empty())return false;
  zFnGraph(out.graph).create(result,connects);path=std::move(out);return true;
}

// Integer intervals with a lower bound of one per protected span. Remaining
// rounding units go to largest remainders, with stable source-order ties.
std::vector<int> allocateIntervals(const std::vector<double> &lengths, int budget) {
  const int n=int(lengths.size());
  if(budget<n) throw std::runtime_error("FaceTag budget cannot preserve anchors");
  std::vector<int> result(n,0); std::vector<double> quota(n);
  int remaining=budget; double length=std::accumulate(lengths.begin(),lengths.end(),0.0);
  bool changed=true;
  while(changed) {
    changed=false;
    for(int i=0;i<n;++i) if(!result[i] && remaining*lengths[i]/length<1) {
      result[i]=1; --remaining; length-=lengths[i]; changed=true;
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
  try {
    if(paths.empty() || paths.size()!=sources.size() ||
        !std::isfinite(settings.ResampleLength) || settings.ResampleLength<zTsField::MinPhysicalParameter)
      throw std::runtime_error("Invalid FaceTag layers or sampling distance");
    std::vector<std::vector<RegionPath>> prepared(paths.size());
    std::vector<std::map<int,double>> lengths(paths.size());
    std::vector<std::map<int,int>> pieces(paths.size()), spans(paths.size());
    std::map<int,int> budgets;
    for(size_t l=0;l<paths.size();++l) {
      if(sources[l].edgeTags.empty()) throw std::runtime_error("Missing FaceTag boundary on layer "+std::to_string(l));
      for(auto path:paths[l]) {
        if(!zTsField::labelFaceTagSegments(path,sources[l])) throw std::runtime_error("Cannot label FaceTag contour");
        if(charts2D && sections3D &&
            !weldProjectedAnchors(path,const_cast<zObjMesh&>(charts2D->at(l)),const_cast<zObjMesh&>(sections3D->at(l))))
          throw std::runtime_error("Cannot normalize projected FaceTag anchors");
        RegionPath p; p.input=std::move(path);
        zFnGraph(p.input.graph).getVertexPositions(p.points);
        std::vector<zTsUtils::GraphPolyline> curves;
        if(!zTsUtils::graphPolylines(p.input.graph,curves) || curves.size()!=1)
          throw std::runtime_error("Disconnected FaceTag path");
        p.closed=curves.front().Closed;
        const int n=int(p.points.size()), e=int(p.input.edgeFaceTags.size());
        if(e!=n-(p.closed?0:1)) throw std::runtime_error("FaceTag connectivity mismatch");
        const auto corners=cornerVertexIds(p.input,settings.FeatureAngleTol);
        int start=0;
        if(p.closed) for(int i=0;i<e;++i) if(p.input.edgeFaceTags[i]!=p.input.edgeFaceTags[(i+e-1)%e]) { start=i; break; }
        for(int k=0;k<e;++k) {
          const int a=(start+k)%n,b=(a+1)%n,tag=p.input.edgeFaceTags[a];
          if(p.runs.empty() || p.runs.back().tag!=tag) {
            RegionRun r;r.tag=tag;r.ids.push_back(a);r.arc.push_back(0);r.anchors.push_back(0);
            p.runs.push_back(std::move(r));
          }
          auto &r=p.runs.back();
          const double d=std::sqrt(zTsUtils::vecDistSq(p.points[a],p.points[b]));
          if(!std::isfinite(d) || d<=0) throw std::runtime_error("Degenerate FaceTag edge");
          r.ids.push_back(b);r.arc.push_back(r.arc.back()+d);
          if(std::find(corners.begin(),corners.end(),b)!=corners.end()) r.anchors.push_back(int(r.ids.size())-1);
        }
        for(auto &r:p.runs) {
          const int end=int(r.ids.size())-1;
          if(r.anchors.back()!=end) r.anchors.push_back(end);
          if(r.tag) { lengths[l][r.tag]+=r.length(); ++pieces[l][r.tag]; spans[l][r.tag]+=int(r.anchors.size())-1; budgets[r.tag]=0; }
        }
        prepared[l].push_back(std::move(p));
      }
    }
    if(budgets.empty()) throw std::runtime_error("No outer FaceTag segments were found");
    for(auto &entry:budgets) {
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
      std::map<int,std::vector<double>> spanLengths;
      for(const auto &p:prepared[l]) for(const auto &r:p.runs) if(r.tag)
        for(size_t a=1;a<r.anchors.size();++a)
          spanLengths[r.tag].push_back(r.arc[r.anchors[a]]-r.arc[r.anchors[a-1]]);
      std::map<int,std::vector<int>> counts;std::map<int,size_t> cursors;
      for(const auto &b:budgets) counts[b.first]=allocateIntervals(spanLengths[b.first],b.second-pieces[l][b.first]);
      for(const auto &p:prepared[l]) {
        PrintPath out;zPointArray points;zIntArray edges;
        const auto append=[&](const RegionRun &r,double arc) {
          auto it=std::upper_bound(r.arc.begin(),r.arc.end(),arc);
          const int k=(std::min)(int(r.ids.size())-2,(std::max)(0,int(it-r.arc.begin())-1));
          const int a=r.ids[k],b=r.ids[k+1];
          const double t=std::clamp((arc-r.arc[k])/(r.arc[k+1]-r.arc[k]),0.0,1.0);
          auto pa=p.points[a],pb=p.points[b];
          const auto point=pa+(pb-pa)*t;
          if(!points.empty() && zTsUtils::vecDistSq(points.back(),point)<1e-18) return int(points.size())-1;
          if(p.closed && arc==r.length() && &r==&p.runs.back() && !points.empty() &&
              zTsUtils::vecDistSq(points.front(),point)<1e-18) return 0;
          const int id=int(points.size());points.push_back(point);out.faceTags.push_back(0);
          out.widths.push_back(float(p.input.widths[a]+(p.input.widths[b]-p.input.widths[a])*t));
          if(!p.input.heights.empty()) out.heights.push_back(float(p.input.heights[a]+(p.input.heights[b]-p.input.heights[a])*t));
          if(!p.input.normals.empty()) { auto na=p.input.normals[a],nb=p.input.normals[b]; auto normal=na+(nb-na)*t; if(normal.length()>0)normal.normalize();out.normals.push_back(normal); }
          if(!p.input.trackIds.empty()) out.trackIds.push_back((p.input.trackIds[a]+(p.input.trackIds[b]-p.input.trackIds[a])*t)>=.5?1.f:0.f);
          return id;
        };
        for(const auto &r:p.runs) {
          PrintPath::FaceTagSegment segment;segment.tag=r.tag;
          int previous=append(r,0);segment.vertices.push_back(previous);
          for(size_t a=1;a<r.anchors.size();++a) {
            const double lo=r.arc[r.anchors[a-1]],hi=r.arc[r.anchors[a]];
            const int count=r.tag ? counts[r.tag].at(cursors[r.tag]++) :
                (std::max)(1,int(std::ceil((hi-lo)/settings.ResampleLength)));
            for(int j=1;j<=count;++j) {
              const double sampleArc=lo+(hi-lo)*j/count;
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
        zFnGraph(out.graph).create(points,edges);result[l].push_back(std::move(out));
      }
    }
    paths=std::move(result);status="FaceTag sampling (2D, fixed total across pieces):";
    for(const auto &b:budgets) status+=" "+std::to_string(b.first)+"="+std::to_string(b.second);
    return true;
  } catch(const std::exception &error) { status=std::string("failed: ")+error.what();return false; }
}
}
