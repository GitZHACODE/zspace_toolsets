#pragma once

#include <zspace/zToolsets/fabrication/slicing/utils/Utils.h>
#include <algorithm>
#include <array>
#include <limits>
#include <numeric>

namespace zSpace { namespace zTsDetail {

// Read-only, per-section nearest-triangle cache for mapped path targets.
// Stored geometry is double precision; only final zSpace positions are float.
class SurfaceQuery {
public:
  using V=Eigen::Vector3d;
  struct Hit {int Triangle=-1;V Position=V::Zero(),Weights=V::Zero();double Distance=0;};
  struct Triangle {std::array<int,3> Ids;V P[3],Normal;};
private:
  struct Box {
    V Lo=V::Constant((std::numeric_limits<double>::infinity)()),Hi=-Lo;
    void add(const V &p){Lo=Lo.cwiseMin(p);Hi=Hi.cwiseMax(p);}
    double distance(const V &p)const{return ((Lo-p).cwiseMax(0)+(p-Hi).cwiseMax(0)).squaredNorm();}
  };
  struct Node {Box Bounds;int Begin=0,End=0,Left=-1,Right=-1;};
  std::vector<Triangle> triangles_;
  std::vector<Box> boxes_;
  std::vector<Node> nodes_;
  std::vector<int> order_;
  static V segment(const V &p,const V &a,const V &b){const V d=b-a;const double n=d.squaredNorm();return a+(n>0?std::clamp((p-a).dot(d)/n,0.0,1.0):0)*d;}
  static V weights(const V &p,const Triangle &t){const V a=t.P[1]-t.P[0],b=t.P[2]-t.P[0],d=p-t.P[0],n=a.cross(b);const double den=n.squaredNorm();const double v=d.cross(b).dot(n)/den,w=a.cross(d).dot(n)/den;return V(1-v-w,v,w);}
  static V closest(const V &p,const Triangle &t){const V projected=p-t.Normal*((p-t.P[0]).dot(t.Normal));
    if(weights(projected,t).minCoeff()>=0)return projected;
    V best=t.P[0];double d=(p-best).squaredNorm();for(int j=0;j<3;++j){const V q=segment(p,t.P[j],t.P[(j+1)%3]);const double next=(q-p).squaredNorm();if(next<d){d=next;best=q;}}return best;}
  int build(int begin,int end){const int id=int(nodes_.size());nodes_.push_back(Node());nodes_[id].Begin=begin;nodes_[id].End=end;
    for(int i=begin;i<end;++i){nodes_[id].Bounds.add(boxes_[order_[i]].Lo);nodes_[id].Bounds.add(boxes_[order_[i]].Hi);}
    if(end-begin>8){Eigen::Index axis;(nodes_[id].Bounds.Hi-nodes_[id].Bounds.Lo).maxCoeff(&axis);const int middle=(begin+end)/2;
      std::nth_element(order_.begin()+begin,order_.begin()+middle,order_.begin()+end,[&](int a,int b){return boxes_[a].Lo[axis]+boxes_[a].Hi[axis]<boxes_[b].Lo[axis]+boxes_[b].Hi[axis];});
      const int left=build(begin,middle),right=build(middle,end);nodes_[id].Left=left;nodes_[id].Right=right;}return id;}
public:
  explicit SurfaceQuery(zObjMesh &mesh){zPointArray p;zIntArray c,n;zFnMesh fn(mesh);fn.getVertexPositions(p);fn.getPolygonData(c,n);size_t offset=0;
    for(int count:n){for(int j=1;j+1<count;++j){Triangle t;t.Ids={c[offset],c[offset+j],c[offset+j+1]};Box box;
        for(int k=0;k<3;++k){t.P[k]=zTsUtils::toEigenVec(p.at(t.Ids[k]));box.add(t.P[k]);}
        t.Normal=(t.P[1]-t.P[0]).cross(t.P[2]-t.P[0]);if(t.Normal.squaredNorm()==0)continue;t.Normal.normalize();triangles_.push_back(t);boxes_.push_back(box);}
      offset+=count;}
    order_.resize(triangles_.size());std::iota(order_.begin(),order_.end(),0);if(!order_.empty())build(0,int(order_.size()));}
  bool nearest(const zPoint &position,Hit &hit)const{if(nodes_.empty())return false;const V p=zTsUtils::toEigenVec(position);double best=(std::numeric_limits<double>::infinity)();
    std::array<int,64> stack;int size=0;stack[size++]=0;
    while(size){const Node &node=nodes_[stack[--size]];if(node.Bounds.distance(p)>best)continue;
      if(node.Left>=0){const double a=nodes_[node.Left].Bounds.distance(p),b=nodes_[node.Right].Bounds.distance(p);
        if(a<b){if(b<=best)stack[size++]=node.Right;if(a<=best)stack[size++]=node.Left;}else{if(a<=best)stack[size++]=node.Left;if(b<=best)stack[size++]=node.Right;}}
      else for(int j=node.Begin;j<node.End;++j){const int id=order_[j];if(boxes_[id].distance(p)>best)continue;const V q=closest(p,triangles_[id]);const double d=(q-p).squaredNorm();
        if(d<best){best=d;hit.Triangle=id;hit.Position=q;hit.Weights=weights(q,triangles_[id]);}}
    }hit.Distance=std::sqrt(best);return hit.Triangle>=0;}
  const Triangle &triangle(int id)const{return triangles_[id];}
};
} }
