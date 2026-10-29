#include <zspace/zToolsets/fabrication/slicing/utils/Utils.h>
#include <zspace/zToolsets/fabrication/slicing/printsyn/zPrintSyn.h>
#include <zspace/zToolsets/fabrication/slicing/utils/Parallel.h>
#include "FieldDistance.h"
#include <algorithm>
#include <cmath>
#include <limits>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace zSpace {
namespace zTsField {
using namespace zTsUtils;
namespace {
using zTsDetail::FastSegment;
using zTsDetail::computeFastSegments;
using zTsDetail::segmentDistanceSquared;

// Contiguous source ranges retain edge order, including the nontransitive
// 1e-8 categorical tie rule. Bounds cull work without reordering candidates.
class PreparedSegments {
  struct Node {
    float minX=0,maxX=0,minY=0,maxY=0;
    int begin=0,end=0,left=-1,right=-1;
  };
  std::vector<FastSegment> _segments;
  std::vector<Node> _nodes;
  zIntArray _tags;
  bool _hasDistances=false;
  int build(int begin,int end) {
    const int id=int(_nodes.size());_nodes.emplace_back();
    Node node;node.begin=begin;node.end=end;
    node.minX=node.maxX=_segments[begin].minX;
    node.minY=node.maxY=_segments[begin].minY;
    for(int i=begin;i<end;++i){const auto &e=_segments[i];
      node.minX=(std::min)(node.minX,e.minX);node.maxX=(std::max)(node.maxX,e.maxX);
      node.minY=(std::min)(node.minY,e.minY);node.maxY=(std::max)(node.maxY,e.maxY);}
    if(end-begin>8){const int middle=begin+(end-begin)/2;
      node.left=build(begin,middle);node.right=build(middle,end);}
    _nodes[id]=node;return id;
  }
  static float lowerSquared(const Node &node,float x,float y) {
    const float dx=(std::max)(0.f,(std::max)(node.minX-x,x-node.maxX));
    const float dy=(std::max)(0.f,(std::max)(node.minY-y,y-node.maxY));
    return dx*dx+dy*dy;
  }
  static bool ownerCannotWin(const Node &node,float x,float y,double best) {
    // Float projection can round outside an exact segment box. Retain a
    // conservative coordinate-scale allowance when pruning categorical owners.
    const double scale=1+std::abs(double(x))+std::abs(double(y))+
      (std::max)({std::abs(double(node.minX)),std::abs(double(node.maxX)),
                 std::abs(double(node.minY)),std::abs(double(node.maxY))})+
      std::abs(double(node.maxX)-node.minX)+std::abs(double(node.maxY)-node.minY);
    const double lower=lowerSquared(node,x,y);
    if(!std::isfinite(lower))return false;
    return std::sqrt(lower)>
      best+1e-8+16*std::numeric_limits<float>::epsilon()*scale;
  }
public:
  PreparedSegments()=default;
  explicit PreparedSegments(std::vector<FastSegment> edges,zIntArray tags={})
      :_segments(std::move(edges)),_tags(std::move(tags)) {
    _hasDistances=std::any_of(_segments.begin(),_segments.end(),[](const auto &e){return e.lenSqInv>0;});
    if(!_segments.empty()){_nodes.reserve(2*_segments.size());build(0,int(_segments.size()));}
  }
  bool empty()const{return _segments.empty();}
  zScalar distance(const zPoint &point)const {
    if(!_hasDistances)return (std::numeric_limits<zScalar>::max)();
    const float x=float(point.x),y=float(point.y);float best=(std::numeric_limits<float>::max)();
    int stack[64],count=0;stack[count++]=0;
    while(count){const auto &node=_nodes[stack[--count]];
      if(lowerSquared(node,x,y)>=best)continue;
      if(node.left>=0){stack[count++]=node.right;stack[count++]=node.left;continue;}
      for(int i=node.begin;i<node.end;++i){const auto &edge=_segments[i];if(edge.lenSqInv==0)continue;
        const float dx=(std::max)(0.f,(std::max)(edge.minX-x,x-edge.maxX));
        const float dy=(std::max)(0.f,(std::max)(edge.minY-y,y-edge.maxY));
        if(dx*dx+dy*dy>=best)continue;
        const float d=segmentDistanceSquared(x,y,edge);if(d<best)best=d;}
    }
    return static_cast<zScalar>(std::sqrt(best));
  }
  void distances(const zPoint *points,size_t count,zScalarArray &out)const {
    out.resize(count);for(size_t i=0;i<count;++i)out[i]=distance(points[i]);
  }
  bool inside(const zPoint &point)const {
    if(_nodes.empty())return false;
    const float x=float(point.x),y=float(point.y);int winding=0,stack[64],count=0;stack[count++]=0;
    while(count){const auto &node=_nodes[stack[--count]];
      if(y<node.minY||y>node.maxY)continue;
      if(node.left>=0){stack[count++]=node.right;stack[count++]=node.left;continue;}
      for(int i=node.begin;i<node.end;++i){const auto &edge=_segments[i];if(y<edge.minY||y>edge.maxY)continue;
        if(edge.y0<=y){if(edge.y1>y){const float cross=edge.dx*(y-edge.y0)-(x-edge.x0)*edge.dy;if(cross>0)++winding;}}
        else if(edge.y1<=y){const float cross=edge.dx*(y-edge.y0)-(x-edge.x0)*edge.dy;if(cross<0)--winding;}}
    }
    return winding!=0;
  }
  int owner(const zPoint &point)const {
    if(!_hasDistances)return 0;
    const float x=float(point.x),y=float(point.y);double best=1e30;int tag=0,stack[64],count=0;stack[count++]=0;
    while(count){const auto &node=_nodes[stack[--count]];if(ownerCannotWin(node,x,y,best))continue;
      if(node.left>=0){stack[count++]=node.right;stack[count++]=node.left;continue;}
      for(int i=node.begin;i<node.end;++i){const auto &edge=_segments[i];if(edge.lenSqInv==0)continue;
        const float dx=(std::max)(0.f,(std::max)(edge.minX-x,x-edge.maxX));
        const float dy=(std::max)(0.f,(std::max)(edge.minY-y,y-edge.maxY));
        float squared=(std::numeric_limits<float>::max)();
        if(dx*dx+dy*dy<squared)squared=(std::min)(squared,segmentDistanceSquared(x,y,edge));
        const double d=static_cast<zScalar>(std::sqrt(squared));const int candidate=_tags[i];
        if(d<best-1e-8||(std::abs(d-best)<=1e-8&&(!tag||candidate<tag))){best=d;tag=candidate;}}
    }
    return isSideFaceTag(tag)?tag:0;
  }
};
}

struct FaceTagQuery::Impl {
  struct Group {PreparedSegments distance,owners;double width=0,half=0;};
  struct Rib {PreparedSegments distance,slot;double offset=0,half=0;};
  PreparedSegments boundary,slot;
  std::vector<Group> groups;
  std::vector<Rib> ribs;
  bool usable=false,doubleBoundary=false;
  double slotHalf=0;
  explicit Impl(const FaceTagSource &source) {
    ScopedSlicingTiming timing("FaceTagQuery.prepare",true);
    zIntArray connects;zFnGraph(const_cast<zObjGraph&>(source.boundary)).getEdgeData(connects);
    if(source.edgeTags.size()!=connects.size()/2||source.edgeWidthIds.size()!=source.edgeTags.size())return;
    std::vector<FastSegment> edges,slotEdges;computeFastSegments(source.boundary,edges);computeFastSegments(source.slot,slotEdges);
    if(edges.size()!=source.edgeTags.size())return;
    boundary=PreparedSegments(edges);slot=PreparedSegments(std::move(slotEdges));doubleBoundary=source.doubleBoundary;
    if(source.slotWidth>0)slotHalf=.5*printedWidth(source.slotWidth,true,true,source.overlap);
    std::vector<FieldSegment> fields;loopWidthSegments(source.boundary,source.edgeWidthIds,source.widths,fields);
    zIntArray ids=source.edgeWidthIds;std::sort(ids.begin(),ids.end());ids.erase(std::unique(ids.begin(),ids.end()),ids.end());
    if(fields.size()!=ids.size())return;
    for(size_t g=0;g<fields.size();++g){Group group;std::vector<FastSegment> distanceEdges,ownerEdges;zIntArray tags;
      computeFastSegments(fields[g].graph,distanceEdges);
      for(size_t e=0;e<edges.size();++e)if(source.edgeWidthIds[e]==ids[g]){ownerEdges.push_back(edges[e]);tags.push_back(source.edgeTags[e]);}
      group.distance=PreparedSegments(std::move(distanceEdges));group.owners=PreparedSegments(std::move(ownerEdges),std::move(tags));
      group.width=fields[g].width;group.half=.5*printedWidth(fields[g].width,true,doubleBoundary,source.overlap);groups.push_back(std::move(group));}
    for(const auto &field:source.bracing){Rib rib;std::vector<FastSegment> distanceEdges,trimEdges;
      computeFastSegments(field.graph,distanceEdges);zObjGraph trim;
      if(field.trim>=0&&field.trim<=1&&createSlot(field.graph,trim,field.trim,field.width*1.5))computeFastSegments(trim,trimEdges);
      rib.distance=PreparedSegments(std::move(distanceEdges));rib.slot=PreparedSegments(std::move(trimEdges));
      rib.half=.5*printedWidth(field.width,false,true,source.overlap);rib.offset=.5*field.width-rib.half;ribs.push_back(std::move(rib));}
    usable=true;
  }
  int label(const zPoint &point)const {
    const bool within=boundary.inside(point);double outer=-1e30,inner=-1e30;int winner=-1;
    for(size_t g=0;g<groups.size();++g){const auto &group=groups[g];const double d=double(group.distance.distance(point))*(within?-1:1);
      if(d+group.half>outer){outer=d+group.half;winner=int(g);}inner=(std::max)(inner,d+group.width-group.half);}
    if(winner<0)return 0;
    if(doubleBoundary){double bracing=1e30;
      for(const auto &rib:ribs){double value=double(rib.distance.distance(point))-rib.offset;
        if(!rib.slot.empty())value=(std::max)(value,-(double(rib.slot.distance(point))-rib.half));bracing=(std::min)(bracing,value);}
      inner=(std::max)(inner,-bracing);if(-inner>outer+1e-8)return 0;
      if(!slot.empty()&&slotHalf>0&&-(double(slot.distance(point))-slotHalf)>outer+1e-8)return 0;}
    return groups[winner].owners.owner(point);
  }
  void query(const zPointArray &points,zIntArray &tags)const {
    // Shared distance buffers mirror computeWidthField's group-wise pass.
    // Original source-edge tags remain independent of equal numerical widths.
    constexpr size_t chunkSize=4096;
    struct Scratch {
      zScalarArray distances,slotDistances;
      std::vector<double> outer,inner,bracing;
      zIntArray winner;
      std::vector<char> within;
    };
    tags.assign(points.size(),0);
    const auto queryChunk=[&](size_t index,Scratch &scratch){
      auto &distances=scratch.distances,&slotDistances=scratch.slotDistances;
      auto &outer=scratch.outer,&inner=scratch.inner,&bracing=scratch.bracing;
      auto &winner=scratch.winner;auto &within=scratch.within;
      const size_t begin=index*chunkSize,n=(std::min)(chunkSize,points.size()-begin);const zPoint *chunk=points.data()+begin;
      outer.assign(n,-1e30);inner.assign(n,-1e30);winner.assign(n,-1);within.resize(n);
      for(size_t i=0;i<n;++i)within[i]=boundary.inside(chunk[i]);
      for(size_t g=0;g<groups.size();++g){const auto &group=groups[g];group.distance.distances(chunk,n,distances);
        for(size_t i=0;i<n;++i){const double d=double(distances[i])*(within[i]?-1:1);
          if(d+group.half>outer[i]){outer[i]=d+group.half;winner[i]=int(g);}inner[i]=(std::max)(inner[i],d+group.width-group.half);}}
      if(doubleBoundary){bracing.assign(n,1e30);
        for(const auto &rib:ribs){rib.distance.distances(chunk,n,distances);if(!rib.slot.empty())rib.slot.distances(chunk,n,slotDistances);
          for(size_t i=0;i<n;++i){double value=double(distances[i])-rib.offset;
            if(!rib.slot.empty())value=(std::max)(value,-(double(slotDistances[i])-rib.half));bracing[i]=(std::min)(bracing[i],value);}}
        for(size_t i=0;i<n;++i){inner[i]=(std::max)(inner[i],-bracing[i]);if(-inner[i]>outer[i]+1e-8)winner[i]=-1;}
        if(!slot.empty()&&slotHalf>0){slot.distances(chunk,n,slotDistances);for(size_t i=0;i<n;++i)if(-(double(slotDistances[i])-slotHalf)>outer[i]+1e-8)winner[i]=-1;}}
      for(size_t i=0;i<n;++i)if(winner[i]>=0)tags[begin+i]=groups[winner[i]].owners.owner(chunk[i]);
    };
    const size_t chunks=(points.size()+chunkSize-1)/chunkSize;
#ifdef _OPENMP
    // Avoid nested OpenMP teams when a caller already parallelizes layers.
    if(chunks>1&&!omp_in_parallel()&&omp_get_max_threads()>1){
      const int workers=(std::min)(int(chunks),omp_get_max_threads());
      zTsDetail::ParallelErrors failures(chunks);
#pragma omp parallel num_threads(workers)
      {
        Scratch scratch;
#pragma omp for schedule(static)
        for(int i=0;i<int(chunks);++i){try{queryChunk(size_t(i),scratch);}catch(...){failures.capture(size_t(i));}}
      }
      failures.rethrow();return;
    }
#endif
    Scratch scratch;for(size_t i=0;i<chunks;++i)queryChunk(i,scratch);
  }
};

FaceTagQuery::FaceTagQuery(const FaceTagSource &source):_impl(std::make_shared<Impl>(source)){}
bool FaceTagQuery::valid()const{return _impl&&_impl->usable;}
int FaceTagQuery::label(const zPoint &point)const{return valid()?_impl->label(point):0;}
bool FaceTagQuery::query(const zPointArray &points,zIntArray &tags)const {
  ScopedSlicingTiming timing("FaceTagQuery.query");tags.clear();if(!valid())return false;_impl->query(points,tags);return true;
}
bool queryFaceTags(const FaceTagSource &source,const zPointArray &points,zIntArray &tags) {
  return FaceTagQuery(source).query(points,tags);
}

bool labelFaceTagSegments(PrintPath &path,const FaceTagSource &source,bool splitTransitions) {
  return labelFaceTagSegments(path,FaceTagQuery(source),splitTransitions);
}

bool labelFaceTagSegments(PrintPath &path,const FaceTagQuery &query,bool splitTransitions) {
  zTsUtils::ScopedSlicingTiming methodTiming("FaceTagSampling.labelFaceTagSegments",true);
  if(!query.valid())return false;
  const auto label=[&](const zPoint &point){return query.label(point);};
  std::vector<GraphPolyline> curves;
  if (!graphPolylines(path.graph, curves) || curves.size()!=1) return false;
  auto ids = curves[0].Vertices;
  if (ids.size()>1 && ids.front()==ids.back()) ids.pop_back();
  const bool closed=curves[0].Closed;
  zPointArray points; zFnGraph(path.graph).getVertexPositions(points);
  const auto validSize=[&](size_t n){return n==0||n==points.size();};
  if (points.size()!=path.widths.size() || ids.size()<2 ||
      !validSize(path.heights.size()) || !validSize(path.normals.size()) ||
      !validSize(path.trackIds.size()) || !validSize(path.featureFlags.size())) return false;
  double coordinateScale=1;
  for(const auto &p:points)coordinateScale=(std::max)(coordinateScale,(std::max)(std::abs(double(p.x)),std::abs(double(p.y))));
  const double anchorTolerance=4*std::numeric_limits<float>::epsilon()*coordinateScale;
  PrintPath result; zPointArray output; zIntArray connects;
  const auto append = [&](int a, int b, double t) {
    const auto p=t==0?points[a]:t==1?points[b]:points[a]+(points[b]-points[a])*t;
    const auto value=[&](float a,float b){return t==0?a:t==1?b:float(a+(b-a)*t);};
    if (!output.empty() && vecDistSq(output.back(),p)<1e-18) return;
    output.push_back(p);
    result.featureFlags.push_back(path.featureFlags.empty()?0:t<=1e-8?path.featureFlags[a]:t>=1-1e-8?path.featureFlags[b]:0);
    result.widths.push_back(value(path.widths[a],path.widths[b]));
    if (!path.heights.empty()) result.heights.push_back(value(path.heights[a],path.heights[b]));
    if (!path.normals.empty()) { auto n=t==0?path.normals[a]:t==1?path.normals[b]:path.normals[a]+(path.normals[b]-path.normals[a])*t; if(n.length()>0) n.normalize(); result.normals.push_back(n); }
    if (!path.trackIds.empty()) result.trackIds.push_back((path.trackIds[a]+(path.trackIds[b]-path.trackIds[a])*t)>=.5 ? 1.f:0.f);
  };
  append(ids[0],ids[0],0);
  for (size_t e=0;e<ids.size()-(closed?0:1);++e) {
    const int a=ids[e], b=ids[(e+1)%ids.size()];
    const auto at=[&](double t){ return label(points[a]+(points[b]-points[a])*t); };
    // Raw contour edges are grid-scale; subdivision also catches two
    // transitions on one edge near a slot/corner. IDs are never interpolated.
    double previous=1e-7; int last=at(previous);
    if(splitTransitions) for(int q=1;q<=8;++q) {
      const double t=q==8 ? 1-1e-7 : double(q)/8;
      const int current=at(t);
      if(current!=last) {
        double lo=previous,hi=t;
        for(int it=0;it<30;++it) { const double mid=.5*(lo+hi); if(at(mid)==last) lo=mid; else hi=mid; }
        double cut=.5*(lo+hi);
        const double length=std::sqrt(vecDistSq(points[a],points[b]));
        // Reuse a geometric feature when a categorical crossing differs only
        // by float precision. Distinct recorded geometric features stay intact.
        if(path.featureFlags.size()==points.size()) {
          if(path.featureFlags[a] && cut*length<=anchorTolerance)cut=0;
          else if(path.featureFlags[b] && (1-cut)*length<=anchorTolerance)cut=1;
        }
        append(a,b,cut);
      }
      previous=t; last=current;
    }
    if (!closed || e+1<ids.size()) append(a,b,1);
  }
  result.faceTags.assign(output.size(),0);
  for(size_t i=0;i<output.size()-(closed?0:1);++i) {
    const int a=int(i),b=int((i+1)%output.size());
    connects.push_back(a);connects.push_back(b);
    const int tag=label((output[a]+output[b])*.5);
    result.edgeFaceTags.push_back(tag);
    for(int v:{a,b}) if(tag && (!result.faceTags[v] || tag<result.faceTags[v])) result.faceTags[v]=tag;
  }
  zFnGraph(result.graph).create(output,connects);
  if(!zPrintSyn::refreshJunctionFeatures(result))return false;
  path=std::move(result); return true;
}


} // namespace zTsField
} // namespace zSpace
