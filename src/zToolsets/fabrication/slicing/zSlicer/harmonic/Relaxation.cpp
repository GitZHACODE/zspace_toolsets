#include <zspace/zToolsets/fabrication/slicing/utils/IsoRelaxation.h>
#include <zspace/zToolsets/fabrication/slicing/utils/Parallel.h>
#include "../../utils/Interpolation.h"
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <Eigen/Eigenvalues>
#include <Eigen/IterativeLinearSolvers>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <cstdlib>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <queue>
#include <set>
#include <sstream>
#include <tuple>

namespace zSpace { namespace zTsRelax {
namespace {
using V = Eigen::Vector3d;
using Tri = std::array<int, 3>;
using Sparse = Eigen::SparseMatrix<double>;
using Key = unsigned long long;
using zTsUtils::toEigenVec;
using zTsUtils::makeEdgeKey;
const double inf = (std::numeric_limits<double>::infinity)();
zPoint point(const V &v) { return zPoint(v.x(), v.y(), v.z()); }

struct Face {
  V p[3];
  V normal() const { return (p[1]-p[0]).cross(p[2]-p[0]); }
};
V barycentric(const V &p, const Face &f) {
  const V a=f.p[1]-f.p[0], b=f.p[2]-f.p[0], n=a.cross(b), d=p-f.p[0];
  const double den=n.squaredNorm();
  if (!(den>0)) return V(-inf,-inf,-inf);
  const double v=d.cross(b).dot(n)/den, w=a.cross(d).dot(n)/den;
  return V(1-v-w,v,w);
}
V closestSegment(const V &p, const V &a, const V &b) {
  const V ab=b-a;
  const double den=ab.squaredNorm();
  return a+(den>0 ? std::clamp((p-a).dot(ab)/den,0.0,1.0) : 0.0)*ab;
}
V closest(const V &p, const Face &f) {
  const V n=f.normal();
  if (n.squaredNorm()>0) {
    V q=p-n*((p-f.p[0]).dot(n)/n.squaredNorm());
    if (barycentric(q,f).minCoeff()>=0) {
      // Preserve exact axis-aligned support planes. Roundoff such as -1e-16
      // on x=0 otherwise fails a strict, valid volume-containment callback.
      for(int c=0;c<3;++c)if(f.p[0][c]==f.p[1][c]&&f.p[0][c]==f.p[2][c])q[c]=f.p[0][c];
      return q;
    }
  }
  V best=f.p[0]; double distance=inf;
  for (int e=0;e<3;++e) {
    const V q=closestSegment(p,f.p[e],f.p[(e+1)%3]);
    const double d=(q-p).squaredNorm();
    if (d<distance) { best=q; distance=d; }
  }
  return best;
}
double segmentDistance(const V &a,const V &b,const V &c,const V &d) {
  const V u=b-a,v=d-c,w=a-c;
  const double aa=u.dot(u),bb=u.dot(v),cc=v.dot(v),dd=u.dot(w),ee=v.dot(w);
  const double den=aa*cc-bb*bb;
  double best=(closestSegment(a,c,d)-a).squaredNorm();
  best=(std::min)(best,(closestSegment(b,c,d)-b).squaredNorm());
  best=(std::min)(best,(closestSegment(c,a,b)-c).squaredNorm());
  best=(std::min)(best,(closestSegment(d,a,b)-d).squaredNorm());
  if (den>1e-15*aa*cc) {
    const double s=(bb*ee-cc*dd)/den,t=(aa*ee-bb*dd)/den;
    if (s>=0 && s<=1 && t>=0 && t<=1) best=(std::min)(best,(w+s*u-t*v).squaredNorm());
  }
  return best;
}
// Clip a coplanar segment against the triangle using its affine barycentrics.
bool segmentInterval(const V &a,const V &b,const Face &f,double tol,double &lo,double &hi,bool localPlane=false) {
  const V n=f.normal(); const double nn=n.norm();
  if (!(nn>0)) return false;
  const double distanceA=(a-f.p[0]).dot(n),distanceB=(b-f.p[0]).dot(n);
  if(!localPlane&&(std::abs(distanceA)>tol*nn||std::abs(distanceB)>tol*nn))return false;
  const V ba=barycentric(a,f),db=barycentric(b,f)-ba;
  lo=0; hi=1;
  for (int j=0;j<3;++j) {
    const double btol=tol*(f.p[(j+1)%3]-f.p[(j+2)%3]).norm()/nn;
    if (std::abs(db[j])<1e-15) { if (ba[j]<-btol) return false; }
    else if (db[j]>0) lo=(std::max)(lo,(-btol-ba[j])/db[j]);
    else hi=(std::min)(hi,(-btol-ba[j])/db[j]);
  }
  if(localPlane){const double delta=distanceB-distanceA,epsilon=tol*nn;
    if(std::abs(delta)<1e-30){if(std::abs(distanceA)>epsilon)return false;}
    else {double a=(-epsilon-distanceA)/delta,b=(epsilon-distanceA)/delta;if(a>b)std::swap(a,b);lo=(std::max)(lo,a);hi=(std::min)(hi,b);}}
  return lo<=hi;
}
struct Contacts {
  std::array<V,12> Points;int Count=0;
  void push_back(const V &p){Points[Count++]=p;}
  bool empty()const{return Count==0;}
  const V *begin()const{return Points.data();}const V *end()const{return Points.data()+Count;}
};
void segmentContacts(const V &a,const V &b,const Face &f,double tol,Contacts &hits) {
  const V n=f.normal(); const double nn=n.norm();
  if (!(nn>0)) return;
  const double da=(a-f.p[0]).dot(n)/nn,db=(b-f.p[0]).dot(n)/nn;
  if (std::abs(da)<=tol && std::abs(db)<=tol) {
    double lo,hi;
    if (segmentInterval(a,b,f,tol,lo,hi)) { hits.push_back(a+lo*(b-a));hits.push_back(a+hi*(b-a)); }
  } else if ((da<=tol && db>=-tol)||(db<=tol && da>=-tol)) {
    const double den=da-db;
    if (std::abs(den)<=1e-30) return;
    const double t=da/den;
    if (t<0 || t>1) return;
    const V q=a+t*(b-a),weights=barycentric(q,f);
    if (weights.minCoeff()>=-1e-8) hits.push_back(q);
  }
}
Contacts contacts(const Face &a,const Face &b,double tol) {
  Contacts hits;
  for (int e=0;e<3;++e) { segmentContacts(a.p[e],a.p[(e+1)%3],b,tol,hits);
    segmentContacts(b.p[e],b.p[(e+1)%3],a,tol,hits); }
  return hits;
}
double distance(const Face &a,const Face &b,double tol) {
  if (!contacts(a,b,tol).empty()) return 0;
  double best=inf;
  for (int i=0;i<3;++i) {
    best=(std::min)(best,(a.p[i]-closest(a.p[i],b)).squaredNorm());
    best=(std::min)(best,(b.p[i]-closest(b.p[i],a)).squaredNorm());
    for (int j=0;j<3;++j) best=(std::min)(best,segmentDistance(a.p[i],a.p[(i+1)%3],b.p[j],b.p[(j+1)%3]));
  }
  return std::sqrt(best);
}

struct Box {
  V lo=V::Constant(inf),hi=V::Constant(-inf);
  void add(const V &p) { lo=lo.cwiseMin(p);hi=hi.cwiseMax(p); }
  void add(const Box &b) { add(b.lo);add(b.hi); }
  double distanceSquared(const V &p) const {
    return ((lo-p).cwiseMax(0)+(p-hi).cwiseMax(0)).squaredNorm();
  }
  double distanceSquared(const Box &b) const {
    return ((lo-b.hi).cwiseMax(0)+(b.lo-hi).cwiseMax(0)).squaredNorm();
  }
  bool overlaps(const Box &b,double tol) const { return (lo.array()<=b.hi.array()+tol).all() && (b.lo.array()<=hi.array()+tol).all(); }
};
Box bounds(const Face &f) { Box b;for (const V &p:f.p) b.add(p);return b; }

// Small immutable BVH shared by projection, gap, ray and intersection queries.
class Tree {
  struct Node { Box box;int begin=0,end=0,left=-1,right=-1; };
  std::vector<Tri> tri_;
  std::vector<Face> faces_;
  std::vector<Box> boxes_;
  std::vector<int> order_; std::vector<Node> nodes_;
  int build(int begin,int end) {
    const int id=int(nodes_.size());nodes_.push_back(Node());
    nodes_[id].begin=begin;nodes_[id].end=end;
    for (int i=begin;i<end;++i) nodes_[id].box.add(boxes_[order_[i]]);
    if (end-begin>8) {
      Eigen::Index axis;(nodes_[id].box.hi-nodes_[id].box.lo).maxCoeff(&axis);
      const int middle=(begin+end)/2;
      std::nth_element(order_.begin()+begin,order_.begin()+middle,order_.begin()+end,[&](int a,int b){
        const Face fa=face(a),fb=face(b);return fa.p[0][axis]+fa.p[1][axis]+fa.p[2][axis]<fb.p[0][axis]+fb.p[1][axis]+fb.p[2][axis];});
      const int left=build(begin,middle),right=build(middle,end);
      nodes_[id].left=left;nodes_[id].right=right;
    }
    return id;
  }
public:
  Tree(const zPointArray &p,const std::vector<Tri> &t):tri_(t) {
    faces_.resize(t.size());boxes_.resize(t.size());
    for(size_t i=0;i<t.size();++i){for(int j=0;j<3;++j)faces_[i].p[j]=toEigenVec(p[t[i][j]]);boxes_[i]=bounds(faces_[i]);}
    order_.resize(t.size());std::iota(order_.begin(),order_.end(),0);
    if (!t.empty()) build(0,int(t.size()));
  }
  void refit(const zPointArray &p){
    for(size_t i=0;i<tri_.size();++i){boxes_[i]=Box();for(int j=0;j<3;++j){faces_[i].p[j]=toEigenVec(p[tri_[i][j]]);boxes_[i].add(faces_[i].p[j]);}}
    for(int i=int(nodes_.size())-1;i>=0;--i){auto &n=nodes_[i];n.box=Box();
      if(n.left>=0){n.box.add(nodes_[n.left].box);n.box.add(nodes_[n.right].box);}
      else for(int j=n.begin;j<n.end;++j)n.box.add(boxes_[order_[j]]);}
  }
  const Face &face(int i) const { return faces_[i]; }
  int size() const { return int(faces_.size()); }
  template<class F> bool visit(const Box &box,double tol,F callback) const {
    if (nodes_.empty()) return true;
    std::array<int,64> pending;int count=0;pending[count++]=0;
    while(count) {const Node &n=nodes_[pending[--count]];
      if(!n.box.overlaps(box,tol))continue;
      if(n.left>=0){pending[count++]=n.left;pending[count++]=n.right;}
      else for(int i=n.begin;i<n.end;++i)if(boxes_[order_[i]].overlaps(box,tol)&&!callback(order_[i]))return false;
    }return true;
  }
  double nearest(const V &p,int &triangle,V &q,int hint=-1) const {
    double best=inf;triangle=-1;std::array<int,64> pending;int count=0;
    if(hint>=0&&hint<int(tri_.size())){q=closest(p,face(hint));best=(q-p).squaredNorm();triangle=hint;}
    if(!nodes_.empty())pending[count++]=0;
    while(count){const Node &n=nodes_[pending[--count]];
      if(n.box.distanceSquared(p)>best)continue;
      if(n.left>=0){const double a=nodes_[n.left].box.distanceSquared(p),b=nodes_[n.right].box.distanceSquared(p);
        const int near=a<=b?n.left:n.right,far=a<=b?n.right:n.left;
        if((std::max)(a,b)<=best)pending[count++]=far;
        if((std::min)(a,b)<=best)pending[count++]=near;}
      else for(int i=n.begin;i<n.end;++i){const int id=order_[i];const V candidate=closest(p,face(id));
        const double d=(candidate-p).squaredNorm();if(d<best){best=d;q=candidate;triangle=id;}}
    }return std::sqrt(best);
  }
  double nearest(const Face &f,double tol) const {
    double best=inf;const Box b=bounds(f);std::array<int,64> pending;int count=0;
    if(!nodes_.empty())pending[count++]=0;
    while(count){const Node &n=nodes_[pending[--count]];
      if(n.box.distanceSquared(b)>best*best)continue;
      if(n.left>=0){pending[count++]=n.left;pending[count++]=n.right;}
      else for(int i=n.begin;i<n.end;++i)best=(std::min)(best,distance(f,face(order_[i]),tol));
    }return best;
  }
  // Global gap uses one best-first dual-tree traversal, rather than restarting
  // an unbounded nearest-triangle search for every triangle in the layer.
  double minimumDistance(const Tree &other,double tol) const {
    if(nodes_.empty()||other.nodes_.empty())return inf;
    struct Pair {double bound;int a,b;bool operator<(const Pair &x)const{return bound>x.bound;}};
    std::priority_queue<Pair> queue;queue.push({nodes_[0].box.distanceSquared(other.nodes_[0].box),0,0});double best=inf;
    const auto push=[&](int a,int b){const double bound=nodes_[a].box.distanceSquared(other.nodes_[b].box);if(bound<best*best)queue.push({bound,a,b});};
    while(!queue.empty()){const Pair item=queue.top();queue.pop();if(item.bound>=best*best)continue;
      const Node &a=nodes_[item.a],&b=other.nodes_[item.b];
      if(a.left<0&&b.left<0){for(int i=a.begin;i<a.end;++i)for(int j=b.begin;j<b.end;++j){const int u=order_[i],v=other.order_[j];
          if(boxes_[u].distanceSquared(other.boxes_[v])>=best*best)continue;
          best=(std::min)(best,distance(face(u),other.face(v),tol));if(best==0)return 0;}}
      else if(b.left<0||(a.left>=0&&a.end-a.begin>=b.end-b.begin)){push(a.left,item.b);push(a.right,item.b);}
      else {push(item.a,b.left);push(item.a,b.right);}
    }return best;
  }
  double ray(const V &p,const V &direction,double maximum,double tol,int *hitTriangle=nullptr,int hint=-1) const {
    double best=maximum;bool found=false;
    if(hint>=0&&hint<int(tri_.size())){const Face f=face(hint);const V n=f.normal();const double den=n.dot(direction);
      if(std::abs(den)>1e-12*n.norm()){const double d=n.dot(f.p[0]-p)/den;
        if(d>tol&&d<best&&barycentric(p+d*direction,f).minCoeff()>=-1e-8){best=d;found=true;if(hitTriangle)*hitTriangle=hint;}}}
    Box box;box.add(p);box.add(p+(best+tol)*direction);
    visit(box,tol,[&](int id){const Face f=face(id);const V n=f.normal();const double den=n.dot(direction);
      if(std::abs(den)<=1e-12*n.norm())return true;
      const double d=n.dot(f.p[0]-p)/den;
      if(d>tol&&d<best&&barycentric(p+d*direction,f).minCoeff()>=-1e-8){best=d;found=true;if(hitTriangle)*hitTriangle=id;}return true;});
    return found?best:inf;
  }
};

struct SidePatch { std::vector<int> triangles;int width=-1; };
struct BoundaryEdge { int a,b,patch; };
struct SideConstraints {
  const Reference &ref;Tree tree;double tol;bool robustPatches=false;
  std::vector<SidePatch> patches;
  std::map<Key,std::vector<std::pair<int,int>>> curves;
  SideConstraints(const Reference &r,double tolerance,bool robust=false):ref(r),tree(r.Positions,r.Triangles),tol(tolerance),robustPatches(robust) {}
  bool build(std::string &error) {
    const int n=int(ref.Triangles.size());std::vector<int> parent(n);std::iota(parent.begin(),parent.end(),0);
    const auto root=[&](int v){while(parent[v]!=v)v=parent[v];return v;};
    std::map<Key,std::vector<int>> edgeFaces;
    for(int t=0;t<n;++t)for(int j=0;j<3;++j)edgeFaces[makeEdgeKey(ref.Triangles[t][j],ref.Triangles[t][(j+1)%3])].push_back(t);
    for(const auto &entry:edgeFaces){const auto &ids=entry.second;if(ids.size()!=2)continue;
      const int a=ids[0],b=ids[1];if(ref.WidthIds[a]<0||ref.WidthIds[a]!=ref.WidthIds[b])continue;
      const Face fa=tree.face(a),fb=tree.face(b);const V na=fa.normal().normalized(),nb=fb.normal().normalized();
      if(std::abs(na.dot(nb))<1-(robustPatches?1e-6:1e-10))continue;
      const Face &plane=robustPatches&&fb.normal().squaredNorm()>fa.normal().squaredNorm()?fb:fa;
      const Face &other=&plane==&fa?fb:fa;const V pn=plane.normal().normalized();
      bool planar=true;for(const V &p:other.p)planar &= std::abs((p-plane.p[0]).dot(pn))<=tol;
      if(planar)parent[root(b)]=root(a);
    }
    std::map<int,int> ids;std::vector<int> patchId(n,-1);
    for(int t=0;t<n;++t)if(ref.WidthIds[t]>=0){const int r=root(t);auto inserted=ids.emplace(r,int(patches.size()));
      if(inserted.second){patches.push_back(SidePatch());patches.back().width=ref.WidthIds[t];}
      patchId[t]=inserted.first->second;patches[patchId[t]].triangles.push_back(t);}
    for(const auto &entry:edgeFaces){const auto &faces=entry.second;if(faces.size()!=2)continue;
      const int a=patchId[faces[0]],b=patchId[faces[1]];if(a<0||b<0||a==b)continue;
      const auto &ta=ref.Triangles[faces[0]];int u=-1,v=-1;
      for(int j=0;j<3;++j)if(makeEdgeKey(ta[j],ta[(j+1)%3])==entry.first){u=ta[j];v=ta[(j+1)%3];break;}
      if(u>=0)curves[makeEdgeKey(a,b)].push_back({u,v});
    }
    if(patches.empty()){error="Input has no tagged side patches.";return false;}return true;
  }
  bool covers(int patch,const V &a,const V &b) const {
    std::vector<std::pair<double,double>> intervals;
    for(int t:patches[patch].triangles){double lo,hi;if(segmentInterval(a,b,tree.face(t),tol,lo,hi,robustPatches))intervals.push_back({lo,hi});}
    std::sort(intervals.begin(),intervals.end());double end=0;
    for(const auto &i:intervals){if(i.first>end+1e-8)return false;end=(std::max)(end,i.second);if(end>=1-1e-8)return true;}return false;
  }
  int support(int width,const V &a,const V &b) const {
    for(size_t p=0;p<patches.size();++p)if(patches[p].width==width&&covers(int(p),a,b))return int(p);
    return -1;
  }
  int uniqueSupport(const V &a,const V &b) const {
    int found=-1;
    for(size_t p=0;p<patches.size();++p)if(covers(int(p),a,b)){
      if(found>=0)return -1;found=int(p);
    }return found;
  }
  V project(const V &candidate,const V &original,const std::set<int> &allowed) const {
    if(allowed.empty())return candidate;
    double best=inf;V result=original;
    if(allowed.size()==1){for(int t:patches[*allowed.begin()].triangles){const V q=closest(candidate,tree.face(t));
      const double d=(q-candidate).squaredNorm();if(d<best){best=d;result=q;}}}
    else if(allowed.size()==2){auto it=allowed.begin();const int a=*it++,b=*it;const auto found=curves.find(makeEdgeKey(a,b));
      if(found!=curves.end())for(const auto &edge:found->second){const V q=closestSegment(candidate,toEigenVec(ref.Positions[edge.first]),toEigenVec(ref.Positions[edge.second]));
        const double d=(q-candidate).squaredNorm();if(d<best){best=d;result=q;}}}
    // Three independent patches meet at a required geometric/semantic corner.
    return result;
  }
};

bool meshValid(const Layer &m,std::string &error) {
  if(m.Positions.empty()||m.Triangles.empty()){error="Empty slice mesh.";return false;}
  for(const auto &p:m.Positions)if(!toEigenVec(p).allFinite()){error="Non-finite slice position.";return false;}
  for(const auto &t:m.Triangles)for(int id:t)if(id<0||size_t(id)>=m.Positions.size()){error="Invalid triangle index.";return false;}
  return true;
}
double triangleQuality(const Layer &m,const Tri &t){const V p[3]={toEigenVec(m.Positions[t[0]]),toEigenVec(m.Positions[t[1]]),toEigenVec(m.Positions[t[2]])};
  double sum=0;for(int j=0;j<3;++j)sum+=(p[j]-p[(j+1)%3]).squaredNorm();return sum>0?(p[1]-p[0]).cross(p[2]-p[0]).norm()/sum:0;}
void repairSlivers(Layer &m,double length){
  // Link-condition collapses remove only tiny interior edges. Keep endpoint
  // coordinates/UVs, so the boundary and all semantic width tags are exact.
  const int n=int(m.Positions.size());std::vector<char> boundary(n,0);
  for(const auto &e:m.BoundaryWidths){boundary[int(e.first>>32)]=boundary[int(e.first&0xffffffffu)]=1;}
  for(int pass=0;pass<4;++pass){std::map<Key,std::vector<int>> edges;std::vector<std::set<int>> neighbours(n);
    for(int f=0;f<int(m.Triangles.size());++f){const auto &t=m.Triangles[f];for(int j=0;j<3;++j){edges[makeEdgeKey(t[j],t[(j+1)%3])].push_back(f);neighbours[t[j]].insert(t[(j+1)%3]);neighbours[t[j]].insert(t[(j+2)%3]);}}
    std::vector<std::pair<double,Key>> shortEdges;for(const auto &e:edges){const int a=int(e.first>>32),b=int(e.first&0xffffffffu);if(boundary[a]||boundary[b]||e.second.size()!=2)continue;
      const double d=(toEigenVec(m.Positions[a])-toEigenVec(m.Positions[b])).norm();if(d<length)shortEdges.push_back({d,e.first});}
    std::sort(shortEdges.begin(),shortEdges.end());std::vector<char> touched(n,0);std::vector<int> parent(n);std::iota(parent.begin(),parent.end(),0);bool changed=false;
    for(const auto &edge:shortEdges){const int a=int(edge.second>>32),b=int(edge.second&0xffffffffu);if(touched[a]||touched[b])continue;
      std::vector<int> shared;std::set_intersection(neighbours[a].begin(),neighbours[a].end(),neighbours[b].begin(),neighbours[b].end(),std::back_inserter(shared));if(shared.size()!=2)continue;
      bool valid=true;for(const auto &t:m.Triangles){if(std::find(t.begin(),t.end(),b)==t.end()||std::find(t.begin(),t.end(),a)!=t.end())continue;
        Tri next=t;for(int &id:next)if(id==b)id=a;const Face old{{toEigenVec(m.Positions[t[0]]),toEigenVec(m.Positions[t[1]]),toEigenVec(m.Positions[t[2]])}},candidate{{toEigenVec(m.Positions[next[0]]),toEigenVec(m.Positions[next[1]]),toEigenVec(m.Positions[next[2]])}};
        if(candidate.normal().dot(old.normal())<=0){valid=false;break;}}
      if(!valid)continue;parent[b]=a;changed=true;for(int id:neighbours[a])touched[id]=1;for(int id:neighbours[b])touched[id]=1;touched[a]=touched[b]=1;
    }
    if(!changed)break;std::vector<Tri> triangles;for(auto t:m.Triangles){for(int &id:t)id=parent[id];if(t[0]!=t[1]&&t[1]!=t[2]&&t[0]!=t[2])triangles.push_back(t);}m.Triangles.swap(triangles);
  }
  // Improve long skinny pairs by flipping their common diagonal. Surface
  // normals may change only slightly; this is a local quality preconditioner.
  for(int pass=0;pass<4;++pass){std::map<Key,std::vector<int>> edges;for(int f=0;f<int(m.Triangles.size());++f)for(int j=0;j<3;++j)edges[makeEdgeKey(m.Triangles[f][j],m.Triangles[f][(j+1)%3])].push_back(f);
    std::set<int> touched;bool changed=false;for(const auto &entry:edges){if(entry.second.size()!=2||m.BoundaryWidths.count(entry.first))continue;
      const int f=entry.second[0],g=entry.second[1];if(touched.count(f)||touched.count(g))continue;const auto oldA=m.Triangles[f],oldB=m.Triangles[g];
      const double quality=(std::min)(triangleQuality(m,oldA),triangleQuality(m,oldB));if(quality>=.01)continue;
      const int u=int(entry.first>>32),v=int(entry.first&0xffffffffu);int a=-1,b=-1;for(int id:oldA)if(id!=u&&id!=v)a=id;for(int id:oldB)if(id!=u&&id!=v)b=id;
      if(a<0||b<0||a==b||edges.count(makeEdgeKey(a,b)))continue;Tri nextA{a,b,u},nextB{b,a,v};
      const auto normal=[&](const Tri &t){return (toEigenVec(m.Positions[t[1]])-toEigenVec(m.Positions[t[0]])).cross(toEigenVec(m.Positions[t[2]])-toEigenVec(m.Positions[t[0]]));};
      V reference=normal(oldA)+normal(oldB);if(normal(nextA).dot(reference)<0){std::swap(nextA[1],nextA[2]);std::swap(nextB[1],nextB[2]);}
      if(normal(nextA).normalized().dot(reference.normalized())<.995||normal(nextB).normalized().dot(reference.normalized())<.995)continue;
      if((std::min)(triangleQuality(m,nextA),triangleQuality(m,nextB))<=1.5*quality)continue;
      m.Triangles[f]=nextA;m.Triangles[g]=nextB;touched.insert(f);touched.insert(g);changed=true;
    }if(!changed)break;
  }
  std::vector<int> map(n,-1);for(const auto &t:m.Triangles)for(int id:t)map[id]=0;int count=0;for(int &id:map)if(id==0)id=count++;
  const auto compact=[&](auto &values){if(values.size()!=size_t(n))return;auto next=values;next.resize(count);for(int i=0;i<n;++i)if(map[i]>=0)next[map[i]]=values[i];values.swap(next);};
  compact(m.Positions);compact(m.UV);compact(m.Heights);compact(m.Normals);for(auto &t:m.Triangles)for(int &id:t)id=map[id];
  auto widths=m.BoundaryWidths;m.BoundaryWidths.clear();for(const auto &e:widths){const int a=map[int(e.first>>32)],b=map[int(e.first&0xffffffffu)];if(a>=0&&b>=0&&a!=b)m.BoundaryWidths[makeEdgeKey(a,b)]=e.second;}
  // A corner triangle whose three vertices lie on the trimmed boundary has
  // no interior geometric freedom. Add one material point before fairing,
  // without moving or subdividing any boundary edge or fixed cap.
  std::set<int> outer;for(const auto &e:m.BoundaryWidths){outer.insert(int(e.first>>32));outer.insert(int(e.first&0xffffffffu));}
  std::vector<Tri> refined;for(const Tri &t:m.Triangles){if(!outer.count(t[0])||!outer.count(t[1])||!outer.count(t[2])){refined.push_back(t);continue;}
    const int id=int(m.Positions.size());const auto average=[&](const zPointArray &values){return point((toEigenVec(values[t[0]])+toEigenVec(values[t[1]])+toEigenVec(values[t[2]]))/3);};
    const zPoint p=average(m.Positions),uv=m.UV.empty()?zPoint():average(m.UV);m.Positions.push_back(p);if(!m.UV.empty())m.UV.push_back(uv);
    if(!m.Heights.empty())m.Heights.push_back((m.Heights[t[0]]+m.Heights[t[1]]+m.Heights[t[2]])/3);
    if(!m.Normals.empty())m.Normals.push_back(point((toEigenVec(m.Normals[t[0]])+toEigenVec(m.Normals[t[1]])+toEigenVec(m.Normals[t[2]])).normalized()));
    for(int j=0;j<3;++j)refined.push_back({t[j],t[(j+1)%3],id});
  }m.Triangles.swap(refined);
}
zVectorArray normals(const Layer &m) {
  std::vector<V> sums(m.Positions.size(),V::Zero());
  std::vector<double> scale(m.Positions.size(),0);
  std::vector<std::vector<int>> neighbours(m.Positions.size());
  double coordinateScale=1;for(const auto &p:m.Positions)coordinateScale=(std::max)(coordinateScale,toEigenVec(p).cwiseAbs().maxCoeff());
  for(const Tri &t:m.Triangles){const V n=(toEigenVec(m.Positions[t[1]])-toEigenVec(m.Positions[t[0]])).cross(toEigenVec(m.Positions[t[2]])-toEigenVec(m.Positions[t[0]]));
    double longest=0;for(int j=0;j<3;++j)longest=(std::max)(longest,(toEigenVec(m.Positions[t[j]])-toEigenVec(m.Positions[t[(j+1)%3]])).squaredNorm());
    for(int j=0;j<3;++j){const int id=t[j];sums[id]+=n;scale[id]=(std::max)(scale[id],longest);
      neighbours[id].push_back(t[(j+1)%3]);neighbours[id].push_back(t[(j+2)%3]);}}
  const auto initial=sums;
  const auto reliable=[&](int i){return initial[i].norm()>(std::max)(1e-6*scale[i],
      16*std::numeric_limits<float>::epsilon()*coordinateScale*std::sqrt(scale[i]));};
  for(size_t i=0;i<sums.size();++i)if(!reliable(int(i))){
    // A vertex supported only by a float-scale sliver has no trustworthy
    // cross-product normal. Recover from its nearest reliable mesh ring.
    std::set<int> visited{int(i)};std::vector<int> ring{int(i)};V recovered=V::Zero();
    while(!ring.empty()&&recovered.norm()==0){std::vector<int> next;
      for(int v:ring)for(int j:neighbours[v])if(visited.insert(j).second){
        if(reliable(j))recovered+=initial[j];else next.push_back(j);}
      ring.swap(next);}
    if(recovered.norm()>0)sums[i]=recovered;
  }
  zVectorArray out;out.reserve(sums.size());for(V &n:sums){if(n.norm()>0)n.normalize();out.push_back(point(n));}return out;
}

struct Work {
  Layer original;
  std::vector<BoundaryEdge> edges;
  std::vector<std::set<int>> allowed;
  std::vector<bool> fixed;
  Eigen::MatrixXd target;
  Eigen::VectorXd mass;
  Sparse bending;
  std::vector<Eigen::Matrix3d> Tangents;
  V planeOrigin=V::Zero(),planeNormal=V::Zero();
  bool planar=false;
  double limit=0;
};

bool prepare(Layer &m,const SideConstraints &side,double gap,const Settings &settings,Work &w,Report &report,std::string &error,bool buildTarget=true) {
  w.original=m;w.allowed.resize(m.Positions.size());w.limit=0.2*gap;
  if(settings.MaxDeviation>0)w.limit=(std::min)(w.limit,settings.MaxDeviation);
  std::map<Key,std::pair<int,std::pair<int,int>>> edges;
  for(const Tri &t:m.Triangles)for(int j=0;j<3;++j){auto &e=edges[makeEdgeKey(t[j],t[(j+1)%3])];++e.first;e.second={t[j],t[(j+1)%3]};}
  for(const auto &item:edges){if(item.second.first>2){error="Non-manifold slice edge.";return false;}
    if(item.second.first!=1)continue;const int a=item.second.second.first,b=item.second.second.second;
    auto width=m.BoundaryWidths.find(item.first);
    int p=-1;
    if(width==m.BoundaryWidths.end()){
      // Float iso coordinates can expose a tiny edge after a degenerate face
      // is discarded. Recover provenance only from unique full-segment cover.
      p=side.uniqueSupport(toEigenVec(m.Positions[a]),toEigenVec(m.Positions[b]));
      if(p<0){error="Slice boundary lacks an unambiguous input-side patch.";return false;}
      width=m.BoundaryWidths.emplace(item.first,side.patches[p].width).first;++report.RecoveredBoundaryTags;
    }else p=side.support(width->second,toEigenVec(m.Positions[a]),toEigenVec(m.Positions[b]));
    if(p<0){std::ostringstream message;message<<"Slice boundary has no conforming input-side patch: width="<<width->second
        <<" ("<<toEigenVec(m.Positions[a]).transpose()<<") -> ("<<toEigenVec(m.Positions[b]).transpose()<<").";error=message.str();return false;}
    w.edges.push_back({a,b,p});w.allowed[a].insert(p);w.allowed[b].insert(p);
  }
  const int n=int(m.Positions.size());Eigen::VectorXd mass=Eigen::VectorXd::Zero(n);
  w.Tangents.assign(n,Eigen::Matrix3d::Identity());
  for(int i=0;i<n;++i){if(w.allowed[i].size()>2){w.Tangents[i].setZero();continue;}
    for(int patch:w.allowed[i]){V normal=V::Zero();double area=0;
      for(int id:side.patches[patch].triangles){const V next=side.tree.face(id).normal();if(next.norm()>area){area=next.norm();normal=next.normalized();}}
      V d=w.Tangents[i]*normal;if(d.norm()>1e-8){d.normalize();w.Tangents[i]-=d*d.transpose();}}
  }
  for(const auto &p:m.Positions)w.planeOrigin+=toEigenVec(p);w.planeOrigin/=n;
  Eigen::Matrix3d covariance=Eigen::Matrix3d::Zero();for(const auto &p:m.Positions){const V d=toEigenVec(p)-w.planeOrigin;covariance+=d*d.transpose();}
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eig(covariance);w.planeNormal=eig.eigenvectors().col(0);w.planar=true;
  for(const auto &p:m.Positions)if(std::abs((toEigenVec(p)-w.planeOrigin).dot(w.planeNormal))>side.tol){w.planar=false;break;}
  w.fixed.assign(n,false);
  std::vector<Eigen::Triplet<double>> coefficients;double sumEdge=0;size_t edgeCount=0;
  for(const Tri &t:m.Triangles){const V p[3]={toEigenVec(m.Positions[t[0]]),toEigenVec(m.Positions[t[1]]),toEigenVec(m.Positions[t[2]])};
    const double area2=(p[1]-p[0]).cross(p[2]-p[0]).norm();if(!(area2>0)){error="Degenerate slice triangle.";return false;}
    double shortest=inf,longest=0;for(int j=0;j<3;++j){const double d=(p[(j+1)%3]-p[j]).squaredNorm();shortest=(std::min)(shortest,d);longest=(std::max)(longest,d);}
    if(shortest<64*side.tol*side.tol||area2<1e-5*longest)for(int id:t)w.fixed[id]=true;
    for(int c=0;c<3;++c){mass[t[c]]+=area2/6;const int a=t[(c+1)%3],b=t[(c+2)%3];
      const double cot=0.5*(p[(c+1)%3]-p[c]).dot(p[(c+2)%3]-p[c])/area2;
      coefficients.emplace_back(a,a,cot);coefficients.emplace_back(b,b,cot);coefficients.emplace_back(a,b,-cot);coefficients.emplace_back(b,a,-cot);
      sumEdge+=(p[(c+1)%3]-p[c]).norm();++edgeCount;}
  }
  if((mass.array()<=0).any()){error="Slice contains an unsupported vertex.";return false;}
  report.SliverConstrainedVertices+=int(std::count(w.fixed.begin(),w.fixed.end(),true));
  Sparse lap(n,n);lap.setFromTriplets(coefficients.begin(),coefficients.end());
  Sparse inv(n,n),lumped(n,n);inv.reserve(n);lumped.reserve(n);
  for(int i=0;i<n;++i){inv.insert(i,i)=1.0/mass[i];lumped.insert(i,i)=mass[i];}
  const double length=sumEdge/edgeCount,lambda=settings.Strength*std::pow(length,4);
  w.mass=mass;w.bending=lambda*(lap.transpose()*inv*lap);
  if(!buildTarget)return true;
  const Sparse matrix=lumped+w.bending;
  Eigen::SimplicialLDLT<Sparse> solver;solver.compute(matrix);
  if(solver.info()!=Eigen::Success){error="Relaxation matrix factorisation failed.";return false;}
  Eigen::MatrixXd rhs(n,3);for(int i=0;i<n;++i)rhs.row(i)=mass[i]*(toEigenVec(m.Positions[i])-w.planeOrigin).transpose();
  w.target=solver.solve(rhs);
  w.target.rowwise()+=w.planeOrigin.transpose();
  if(solver.info()!=Eigen::Success||!w.target.allFinite()){error="Relaxation solve returned invalid positions.";return false;}return true;
}

bool validTrial(const Layer &m,const Work &w,const SideConstraints &side,const SampleVolume &sample,double tol,std::string &reason,std::set<int> &blocked,
                const std::vector<char> *changedVertices=nullptr,const std::vector<char> *changedFaces=nullptr,Tree *preparedTree=nullptr) {
  blocked.clear();
  // Projection/containment need a float-coordinate tolerance. Reusing that
  // expanded slab for intersection predicates invents overlaps between nearly
  // coplanar neighbouring triangles. Intersect the stored geometry in double.
  const double intersectionTol=tol*1e-3;
  for(size_t i=0;i<m.Positions.size();++i)if((!changedVertices||(*changedVertices)[i])&&!sample(m.Positions[i],nullptr)){reason="volume containment";blocked.insert(int(i));return false;}
  for(const auto &edge:w.edges)if((!changedVertices||(*changedVertices)[edge.a]||(*changedVertices)[edge.b])&&!side.covers(edge.patch,toEigenVec(m.Positions[edge.a]),toEigenVec(m.Positions[edge.b]))){double da=inf,db=inf;for(int id:side.patches[edge.patch].triangles){da=(std::min)(da,(closest(toEigenVec(m.Positions[edge.a]),side.tree.face(id))-toEigenVec(m.Positions[edge.a])).norm());db=(std::min)(db,(closest(toEigenVec(m.Positions[edge.b]),side.tree.face(id))-toEigenVec(m.Positions[edge.b])).norm());}
    std::ostringstream msg;msg<<"boundary support "<<edge.a<<","<<edge.b<<" patch="<<edge.patch<<" triangles="<<side.patches[edge.patch].triangles.size()<<" distances="<<da<<","<<db<<" tolerance="<<tol;
    reason=msg.str();blocked.insert(edge.a);blocked.insert(edge.b);return false;}
  std::unique_ptr<Tree> ownTree;if(!preparedTree)ownTree.reset(new Tree(m.Positions,m.Triangles));Tree &tree=preparedTree?*preparedTree:*ownTree;
  for(size_t i=0;i<m.Triangles.size();++i){if(changedFaces&&!(*changedFaces)[i])continue;const Tri &t=m.Triangles[i];const Face f=tree.face(int(i));
    const Face old={ {toEigenVec(w.original.Positions[t[0]]),toEigenVec(w.original.Positions[t[1]]),toEigenVec(w.original.Positions[t[2]])} };
    if(f.normal().squaredNorm()<1e-12*old.normal().squaredNorm()||f.normal().dot(old.normal())<=0){reason="triangle collapse/reversal";blocked.insert(t.begin(),t.end());return false;}
    // Check the complete triangle against the immutable volume boundary. Contacts
    // are allowed only at its genuine side-boundary vertices/edges.
    if(!side.tree.visit(bounds(f),tol,[&](int r){
      if(side.robustPatches){const Face &boundary=side.tree.face(r);const V n=boundary.normal().normalized();double low=inf,high=-inf;
        for(const V &p:f.p){const double d=(p-boundary.p[0]).dot(n);low=(std::min)(low,d);high=(std::max)(high,d);}
        // Contact within coordinate precision is grazing, not a crossing of
        // the volume boundary. True crossings straddle its plane beyond tol.
        if(low>=-tol||high<=tol)return true;
      }
      for(const V &hit:contacts(f,side.tree.face(r),intersectionTol)){
        bool boundary=false;
        for(int j=0;j<3;++j){const int a=t[j],b=t[(j+1)%3];
          if(!w.allowed[a].empty()&&(hit-f.p[j]).norm()<=4*tol)boundary=true;
          if(m.BoundaryWidths.count(makeEdgeKey(a,b))&&(closestSegment(hit,f.p[j],f.p[(j+1)%3])-hit).norm()<=4*tol)boundary=true;}
        if(!boundary)return false;}return true;})){reason="input-surface contact";blocked.insert(t.begin(),t.end());return false;}
    if(!tree.visit(bounds(f),tol,[&](int j){if(j==int(i)||(j<int(i)&&(!changedFaces||(*changedFaces)[j])))return true;const Tri &b=m.Triangles[j];
      const auto hits=contacts(f,tree.face(j),intersectionTol);if(hits.empty())return true;
      std::vector<V> shared;
      for(int a:t)for(int v:b)if(a==v||
          (toEigenVec(w.original.Positions[a])-toEigenVec(w.original.Positions[v])).norm()<=2*tol){
        const V pa=toEigenVec(m.Positions[a]),pb=toEigenVec(m.Positions[v]);
        if((pa-pb).norm()<=2*tol)shared.push_back((pa+pb)*0.5);}
      // Float marching-tet vertices can coincide geometrically while retaining
      // different edge provenance. Permit their existing local contact only;
      // also check adjacent triangles for overlap beyond their shared simplex.
      for(const V &hit:hits){double d=inf;for(const V &p:shared)d=(std::min)(d,(hit-p).norm());
        for(size_t a=0;a<shared.size();++a)for(size_t b=a+1;b<shared.size();++b)
          d=(std::min)(d,(hit-closestSegment(hit,shared[a],shared[b])).norm());
        if(d>4*tol){std::ostringstream message;message<<"self intersection triangles "<<i<<","<<j<<" distance="<<d;reason=message.str();blocked.insert(t.begin(),t.end());blocked.insert(b.begin(),b.end());return false;}}
      return true;}))return false;
  }return true;
}

bool separated(const Layer &a,const Layer &b,double tol) {
  Tree tree(b.Positions,b.Triangles);Tree other(a.Positions,a.Triangles);
  for(size_t i=0;i<a.Triangles.size();++i){const Face f=other.face(int(i));
    if(!tree.visit(bounds(f),tol,[&](int j){return contacts(f,tree.face(j),tol).empty();}))return false;}
  return true;
}
bool separatedTrees(const Tree &a,const Tree &b,double tol,std::pair<int,int> *blocked=nullptr){
  for(int i=0;i<a.size();++i){const Face &f=a.face(i);
    if(!b.visit(bounds(f),tol,[&](int j){const bool clear=contacts(f,b.face(j),tol).empty();if(!clear&&blocked)*blocked={i,j};return clear;}))return false;}
  return true;
}

bool advanceLayer(const Layer &current,Layer &out,Work &w,const SideConstraints &side,
                  const SampleVolume &sample,double tol,Report &report,std::set<int> &uvConstrained) {
  Layer candidate=current;
  for(double step=0.5;step>=1.0/1024;step*=0.5){
    for(int repair=0;repair<32;++repair){candidate.Positions=current.Positions;candidate.UV=current.UV;double moved=0;
      for(size_t i=0;i<candidate.Positions.size();++i){if(w.fixed[i])continue;
        const V old=toEigenVec(current.Positions[i]),origin=toEigenVec(w.original.Positions[i]);
        V desired=old+step*(w.target.row(int(i)).transpose()-old);const V delta=desired-origin;
        if(delta.norm()>w.limit)desired=origin+delta.normalized()*w.limit;
        if(w.planar)desired-=w.planeNormal*((desired-w.planeOrigin).dot(w.planeNormal));
        desired=side.project(desired,origin,w.allowed[i]);
        if((desired-origin).norm()>w.limit+tol){w.fixed[i]=true;++report.GeometryConstrainedVertices;continue;}
        const zPoint p=point(desired);zPoint uv;
        if((desired-old).norm()>tol){if(!sample(p,&uv)){uvConstrained.insert(int(i));continue;}
          if(candidate.UV.size()==candidate.Positions.size())candidate.UV[i]=uv;}
        candidate.Positions[i]=p;moved=(std::max)(moved,(toEigenVec(p)-old).norm());
      }
      std::set<int> blocked;
      if(!validTrial(candidate,w,side,sample,tol,report.LastRejectedReason,blocked)){
        bool constrained=false;for(int id:blocked)if(!w.fixed[id]){w.fixed[id]=true;++report.GeometryConstrainedVertices;constrained=true;}
        if(constrained)continue;break;
      }
      out=std::move(candidate);return moved>tol;
    }
  }return false;
}

struct HeightLink {
  int LayerId=-1;Tri Ids{};V Weights=V::Zero(),PlaneNormal=V::Zero(),Direction=V::Zero(),Position=V::Zero();
  double Height=0,NormalLength=1;bool Fallback=false;
  int Triangle=-1;
};
using HeightLinks=std::vector<std::vector<std::vector<HeightLink>>>;
struct HeightQueryFailure {int LayerId=-1,Vertex=-1,Adjacent=-1;};
bool refreshHeights(std::vector<Layer> &layers,double tol,Report &report,std::string &error,bool parallel,
                    HeightLinks *links=nullptr,const std::vector<char> *refreshMask=nullptr,HeightQueryFailure *queryFailure=nullptr,
                    std::vector<std::unique_ptr<Tree>> *cache=nullptr,bool preparedTrees=false) {
  zTsUtils::ScopedSlicingTiming methodTiming("Height.refreshHeights");
  std::vector<HeightQueryFailure> queryFailures(layers.size());
  if(links&&links->size()!=layers.size())links->resize(layers.size());
  std::vector<std::unique_ptr<Tree>> localTrees;
  if(cache){if(cache->size()!=layers.size())cache->resize(layers.size());}else localTrees.resize(layers.size());
  auto &trees=cache?*cache:localTrees;std::vector<zVectorArray> n(layers.size());
  zTsDetail::ParallelErrors preparation(layers.size());
  {
    zTsUtils::ScopedSlicingTiming geometryTiming("Height.refreshGeometry");
    #pragma omp parallel for schedule(dynamic) if(parallel)
    for(int k=0;k<int(layers.size());++k){try{
      const bool needed=!refreshMask||(*refreshMask)[k]||(k>0&&(*refreshMask)[k-1])||(k+1<int(layers.size())&&(*refreshMask)[k+1]);
      if(needed){if(!trees[k])trees[k].reset(new Tree(layers[k].Positions,layers[k].Triangles));
        else if(!preparedTrees)trees[k]->refit(layers[k].Positions);
        n[k]=normals(layers[k]);}
    }catch(...){preparation.capture(k);}}
    preparation.rethrow();
  }
  std::vector<int> planeEstimates(layers.size(),0),unavailableAuxiliary(layers.size(),0);std::vector<std::string> errors(layers.size());
  const auto refresh=[&](size_t k){std::string &error=errors[k];Layer &m=layers[k];m.Heights.resize(m.Positions.size());
    std::vector<std::vector<HeightLink>> previous;
    if(links){previous=std::move((*links)[k]);(*links)[k].resize(m.Positions.size());}
    if(k>0)m.Normals=n[k];
    for(size_t i=0;i<m.Positions.size();++i){const V p=toEigenVec(m.Positions[i]),normal=toEigenVec(n[k][i]);
      if(queryFailure)queryFailures[k]={int(k),int(i),-1};
      if(normal.norm()<0.5){error="Cannot recover a finite layer normal.";return false;}
      double primaryHeight=0;
      for(int offset:{-1,1}){const int adjacent=int(k)+offset;if(adjacent<0||adjacent>=int(layers.size()))continue;
        const bool primary=k>0&&adjacent==int(k)-1;
        if(queryFailure)queryFailures[k].Adjacent=adjacent;
        int hint=-1;if(i<previous.size())for(const auto &old:previous[i])if(old.LayerId==adjacent){hint=old.Triangle;break;}
        int tri=-1;V q;const double nearest=trees[adjacent]->nearest(p,tri,q,hint);
        if(!(nearest>tol)||tri<0){error="Relaxed layers have zero or invalid separation.";return false;}
        V direction=offset*normal;
        // The bottom cap's input winding is outward (opposite the slicing
        // direction); an imported cap may also use the opposite convention.
        if((k==0||k+1==layers.size())&&direction.dot(q-p)<0)direction=-direction;
        int hitTriangle=-1;
        double h=trees[adjacent]->ray(p,direction,4*nearest+tol,tol,links?&hitTriangle:nullptr,hint);
        V planeNormal=V::Zero(),weights=V::Zero();double normalLength=1;bool fallback=false;
        if(!std::isfinite(h)){
          // A side-boundary normal ray can leave the solid before hitting the
          // adjacent trimmed layer. Extrapolate its nearest local face plane;
          // this is computed from relaxed geometry, not the old field height.
          const Face face=trees[adjacent]->face(tri);
          const V localWeights=barycentric(q,face);
          V faceNormal=V::Zero();for(int j=0;j<3;++j)faceNormal+=localWeights[j]*toEigenVec(n[adjacent][layers[adjacent].Triangles[tri][j]]);
          if(links)normalLength=faceNormal.norm();faceNormal.normalize();
          const double den=faceNormal.dot(direction);
          if(!std::isfinite(den)||std::abs(den)<0.25){
            if(!primary){++unavailableAuxiliary[k];continue;}
            error="Previous layer is nearly tangent to the TCP height ray on layer "+std::to_string(k)+", vertex "+std::to_string(i)+", neighbour "+std::to_string(adjacent)+".";return false;}
          h=faceNormal.dot(q-p)/den;++planeEstimates[k];
          if(links){hitTriangle=tri;planeNormal=faceNormal;weights=barycentric(q,face);fallback=true;}
        }
        if(!std::isfinite(h)||h<=tol){if(!primary){++unavailableAuxiliary[k];continue;}
          error="Cannot estimate positive previous-layer TCP spacing on layer "+std::to_string(k)+", vertex "+std::to_string(i)+", neighbour "+std::to_string(adjacent)+".";return false;}
        // TCP height is the directed clearance to the previously printed
        // layer along -normal. S0 is an unprinted substrate reference; its
        // optional next-gap value is diagnostic and never a height constraint.
        if(primary||k==0)primaryHeight=h;
        if(links){if(planeNormal.squaredNorm()==0){const Face &face=trees[adjacent]->face(hitTriangle);planeNormal=face.normal().normalized();weights=barycentric(p+h*direction,face);}
          weights=weights.cwiseMax(0);weights/=weights.sum();
          (*links)[k][i].push_back({adjacent,layers[adjacent].Triangles[hitTriangle],weights,planeNormal,direction,fallback?q:p+h*direction,h,normalLength,fallback,hitTriangle});}
      }
      if(k>0&&!(primaryHeight>tol)){error="Missing previous-layer TCP height.";return false;}
      m.Heights[i]=float(primaryHeight);
    }
    return true;
  };
  zTsUtils::ScopedSlicingTiming correspondenceTiming("Height.correspondences");
  zTsDetail::ParallelErrors failures(layers.size());
  #pragma omp parallel for schedule(dynamic) if(parallel)
  for(int k=0;k<int(layers.size());++k){try{if(!refreshMask||(*refreshMask)[k])refresh(size_t(k));}catch(...){failures.capture(k);}}
  failures.rethrow();
  for(size_t k=0;k<layers.size();++k){report.PlaneHeightEstimates+=planeEstimates[k];report.UnavailableAuxiliaryGaps+=unavailableAuxiliary[k];
    if(!errors[k].empty()){error=errors[k];if(queryFailure)*queryFailure=queryFailures[k];return false;}}
  return true;
}
} // namespace

bool relax(const Reference &reference,std::vector<Layer> &layers,const Settings &settings,
           const SampleVolume &sample,Report &report,std::string &error) {
  report=Report();error.clear();
  using Clock=std::chrono::steady_clock;const auto start=Clock::now();auto last=start;
  const auto mark=[&](int stage){const auto now=Clock::now();report.StageMilliseconds[stage]=std::chrono::duration<double,std::milli>(now-last).count();last=now;};
  if(settings.Iterations<0||settings.Iterations>100||!std::isfinite(settings.Strength)||settings.Strength<=0||settings.Strength>1||
     !std::isfinite(settings.MaxDeviation)||settings.MaxDeviation<0||!std::isfinite(settings.MinEdgeLength)||settings.MinEdgeLength<0){error="Invalid iso relaxation settings (iterations 0..100, strength (0,1], deviation/edge length >=0).";return false;}
  if(settings.Iterations==0)return true;
  if(layers.size()<2||!sample||reference.Triangles.size()!=reference.WidthIds.size()){error="Missing layer, reference or volume-query data.";return false;}
  Layer ref;ref.Positions=reference.Positions;ref.Triangles=reference.Triangles;if(!meshValid(ref,error))return false;
  for(const Layer &m:layers)if(!meshValid(m,error))return false;
  Box box;for(const auto &p:reference.Positions)box.add(toEigenVec(p));
  const double tol=(std::max)(1e-9,2e-7*(box.hi-box.lo).norm());
  SideConstraints side(reference,tol,settings.MinEdgeLength>0);if(!side.build(error))return false;
  if(settings.MinEdgeLength>0)for(size_t k=1;k+1<layers.size();++k){Layer candidate=layers[k];repairSlivers(candidate,settings.MinEdgeLength);
    Work check;Report checkReport;std::string reason;std::set<int> blocked;
    if(prepare(candidate,side,1,settings,check,checkReport,reason)&&validTrial(candidate,check,side,sample,tol,reason,blocked))layers[k]=std::move(candidate);
    else zTsUtils::printDebug("[heightOptimizer] Retained original triangulation on layer "+std::to_string(k)+": "+reason);
  }
  std::vector<std::unique_ptr<Tree>> initial;for(const Layer &m:layers)initial.emplace_back(new Tree(m.Positions,m.Triangles));
  mark(0);
  std::vector<double> gap(layers.size(),inf);report.MinimumInitialGap=inf;
  std::vector<double> pairGap(layers.size()-1,inf);
  zTsDetail::ParallelErrors gapFailures(pairGap.size());
  #pragma omp parallel for schedule(dynamic) if(settings.Parallel)
  for(int k=0;k<int(pairGap.size());++k){try{pairGap[k]=initial[k]->minimumDistance(*initial[k+1],tol);}catch(...){gapFailures.capture(k);}}
  gapFailures.rethrow();
  for(size_t k=0;k<pairGap.size();++k){const double d=pairGap[k];
    if(!std::isfinite(d)||d<=4*tol){error="Initial neighbouring layers touch or have insufficient separation.";return false;}
    gap[k]=(std::min)(gap[k],d);gap[k+1]=(std::min)(gap[k+1],d);report.MinimumInitialGap=(std::min)(report.MinimumInitialGap,d);
  }
  std::vector<Work> work(layers.size());
  std::vector<Report> localReports(layers.size());std::vector<std::string> prepareErrors(layers.size());
  mark(1);
  zTsDetail::ParallelErrors prepareFailures(layers.size());
  #pragma omp parallel for schedule(dynamic) if(settings.Parallel)
  for(int k=1;k<int(layers.size())-1;++k){try{prepare(layers[k],side,gap[k],settings,work[k],localReports[k],prepareErrors[k]);}catch(...){prepareFailures.capture(k);}}
  prepareFailures.rethrow();
  for(size_t k=1;k+1<layers.size();++k)if(!prepareErrors[k].empty()){error="Layer "+std::to_string(k)+": "+prepareErrors[k];return false;}
  std::vector<std::set<int>> uvConstrained(layers.size());
  mark(2);
  // Geometrically separated meshes can still have unusable directed height
  // rays at a sharp boundary. Remember whether the seed is measurable, so a
  // fairing step cannot discard the last usable complete stack.
  bool measurable=false;
  {auto measured=layers;Report ignored;std::string reason;
    measurable=refreshHeights(measured,tol,ignored,reason,settings.Parallel);}
  for(int iteration=0;iteration<settings.Iterations;++iteration){std::vector<Layer> trial=layers;std::vector<char> moved(layers.size(),0),crossed(layers.size()-1,0);
    zTsDetail::ParallelErrors trialFailures(layers.size());
    #pragma omp parallel for schedule(dynamic) if(settings.Parallel)
    for(int k=1;k<int(layers.size())-1;++k){try{moved[k]=advanceLayer(layers[k],trial[k],work[k],side,sample,tol,localReports[k],uvConstrained[k]);}catch(...){trialFailures.capture(k);}}
    trialFailures.rethrow();
    zTsDetail::ParallelErrors separationFailures(crossed.size());
    #pragma omp parallel for schedule(dynamic) if(settings.Parallel)
    for(int k=0;k<int(crossed.size());++k){try{crossed[k]=!separated(trial[k],trial[k+1],tol);}catch(...){separationFailures.capture(k);}}
    separationFailures.rethrow();
    const bool accepted=std::find(moved.begin(),moved.end(),char(1))!=moved.end();
    const bool valid=std::find(crossed.begin(),crossed.end(),char(1))==crossed.end();
    if(!valid)report.LastRejectedReason="adjacent layer intersection";
    if(!valid||!accepted)break;
    Report heightReport;std::string heightError;bool trialMeasurable=false;
    {zTsUtils::ScopedSlicingTiming feasibilityTiming("Relaxation.heightFeasibility");
      trialMeasurable=refreshHeights(trial,tol,heightReport,heightError,settings.Parallel);}
    if(!trialMeasurable&&measurable){report.LastRejectedReason="Rejected unmeasurable fairing step: "+heightError;break;}
    measurable=trialMeasurable;
    layers.swap(trial);++report.AcceptedIterations;
  }
  for(size_t k=1;k+1<layers.size();++k){report.UVConstrainedVertices+=int(uvConstrained[k].size());
    report.GeometryConstrainedVertices+=localReports[k].GeometryConstrainedVertices;report.SliverConstrainedVertices+=localReports[k].SliverConstrainedVertices;
    report.RecoveredBoundaryTags+=localReports[k].RecoveredBoundaryTags;
    if(!localReports[k].LastRejectedReason.empty())report.LastRejectedReason=localReports[k].LastRejectedReason+" on layer "+std::to_string(k);}
  mark(3);
  for(size_t k=1;k+1<layers.size();++k)for(size_t i=0;i<layers[k].Positions.size();++i){const V p=toEigenVec(layers[k].Positions[i]);
    const double moved=(p-toEigenVec(work[k].original.Positions[i])).norm();report.MaxDisplacement=(std::max)(report.MaxDisplacement,moved);
    if(!work[k].allowed[i].empty()){if(moved>tol)++report.MovedBoundaryVertices;int t;V q;
      report.MaxBoundaryError=(std::max)(report.MaxBoundaryError,side.tree.nearest(p,t,q));}}
  const bool ok=refreshHeights(layers,tol,report,error,settings.Parallel);mark(4);
  report.StageMilliseconds[5]=std::chrono::duration<double,std::milli>(Clock::now()-start).count();return ok;
}

namespace {
bool validHeightSettings(const PrintHeightOptimizationSettings &s,std::string &error) {
  if(!std::isfinite(s.Minimum)||!std::isfinite(s.Maximum)||s.Minimum<=0||s.Maximum<=s.Minimum||
     s.Iterations<1||s.Iterations>32||!std::isfinite(s.MaxDisplacement)||s.MaxDisplacement<0||
     !std::isfinite(s.Tolerance)||s.Tolerance<0||s.Tolerance>=s.Minimum){
    error="Height domain requires finite 0 < minimum < maximum, iterations 1..32, nonnegative displacement and 0 <= tolerance < minimum.";return false;}
  return true;
}
double bandViolation(double h,const PrintHeightOptimizationSettings &s) {
  return (std::max)({0.0,s.Minimum-h,h-s.Maximum});
}
struct HeightSample {
  Tri Ids{};V Weights=V::Zero();bool Primary=true,Exact=false;
  HeightLink Hit;V Position=V::Zero(),Normal=V::Zero();double NormalLength=1;
};
using HeightSamples=std::vector<std::vector<HeightSample>>;
double sampleHeight(const Layer &layer,const HeightSample &sample){return sample.Exact?sample.Hit.Height:zTsDetail::interpolateScalar(layer.Heights[sample.Ids[0]],layer.Heights[sample.Ids[1]],layer.Heights[sample.Ids[2]],sample.Weights[0],sample.Weights[1],sample.Weights[2]);}

bool refreshExactSamples(const std::vector<Layer> &layers,HeightSamples &samples,
    const std::vector<std::unique_ptr<Tree>> &trees,double maximum,double tol,
    std::string &error,bool parallel,const std::vector<char> *mask=nullptr,HeightQueryFailure *failure=nullptr) {
  zTsUtils::ScopedSlicingTiming timing("Height.refreshExactSamples");
  std::vector<std::string> errors(layers.size());zTsDetail::ParallelErrors failures(layers.size());
  #pragma omp parallel for schedule(dynamic) if(parallel)
  for(int k=1;k<int(layers.size());++k){try{if(mask&&!(*mask)[k])continue;
    for(auto &sample:samples[k]){if(!sample.Exact)continue;sample.Position.setZero();V sum=V::Zero();
      for(int j=0;j<3;++j){sample.Position+=sample.Weights[j]*toEigenVec(layers[k].Positions[sample.Ids[j]]);
        sum+=sample.Weights[j]*toEigenVec(layers[k].Normals[sample.Ids[j]]);}
      sample.NormalLength=sum.norm();if(!sum.allFinite()||sample.NormalLength<1e-8){errors[k]="Invalid interpolated TCP normal.";break;}
      sample.Normal=sum/sample.NormalLength;const V direction=-sample.Normal;
      int nearestTri=-1;V nearest;const double clearance=trees[k-1]->nearest(sample.Position,nearestTri,nearest,sample.Hit.Triangle);
      if(!(clearance>tol)){errors[k]="TCP touches its previously printed surface.";break;}
      int triangle=-1;const double height=trees[k-1]->ray(sample.Position,direction,maximum,tol,&triangle,sample.Hit.Triangle);
      if(!std::isfinite(height)||triangle<0){errors[k]="TCP normal ray misses its previously printed surface on section "+std::to_string(k)+".";break;}
      const Face &face=trees[k-1]->face(triangle);V weights=barycentric(sample.Position+height*direction,face).cwiseMax(0);weights/=weights.sum();
      sample.Hit={k-1,layers[k-1].Triangles[triangle],weights,face.normal().normalized(),direction,
        sample.Position+height*direction,height,1,false,triangle};
    }
  }catch(...){failures.capture(k);}}
  failures.rethrow();for(int k=1;k<int(errors.size());++k)if(!errors[k].empty()){
    error=errors[k];if(failure)*failure={k,-1,k-1};return false;}
  return true;
}
double heightScore(const std::vector<Layer> &layers,const HeightSamples &samples,const PrintHeightOptimizationSettings &s,HeightReport &r) {
  r.Minimum=inf;r.Maximum=0;r.WorstViolation=0;r.ViolatingVertices=0;double score=0;
  // TCPs occupy S1..Slast. S0 is the fixed substrate reference; the fixed top
  // contact is the final printed TCP surface and has a previous-layer gap.
  for(size_t k=1;k<layers.size();++k)for(const auto &point:samples[k]){const double h=sampleHeight(layers[k],point);
    const double v=bandViolation(h,s),residual=(std::max)(0.0,v-.5*s.Tolerance);
    score+=residual*residual;if(point.Primary){r.Minimum=(std::min)(r.Minimum,double(h));r.Maximum=(std::max)(r.Maximum,double(h));
      r.WorstViolation=(std::max)(r.WorstViolation,v);if(v>s.Tolerance)++r.ViolatingVertices;}}
  r.Converged=r.ViolatingVertices==0;return score;
}
}

bool redistributeIsoValues(std::vector<Layer> layers,const PrintHeightOptimizationSettings &s,
                           std::vector<double> &levels,bool parallel,std::string &error) {
  if(!s.Enabled)return true;
  if(!validHeightSettings(s,error)||layers.size()<3){if(error.empty())error="Height scheduling needs at least three sections.";return false;}
  for(const auto &m:layers)if(!meshValid(m,error))return false;
  Box bounds;for(const auto &m:layers)for(const auto &p:m.Positions)bounds.add(toEigenVec(p));
  const double tol=(std::max)(1e-9,2e-7*(bounds.hi-bounds.lo).norm());Report ignored;
  if(!refreshHeights(layers,tol,ignored,error,parallel))return false;
  const int n=int(layers.size())-1;Eigen::VectorXd scale=Eigen::VectorXd::Ones(n);
  const double margin=(std::min)(.02*(s.Maximum-s.Minimum),.1*s.Minimum);
  std::vector<std::pair<double,double>> ranges(n);
  for(int k=0;k<n;++k){auto h=layers[k+1].Heights;std::sort(h.begin(),h.end());
    // Robust initializer only: isolated nodal slivers must not displace an
    // entire layer. Every vertex/path is still checked by the geometric stage.
    ranges[k]={h[size_t(.01*(h.size()-1))],h[size_t(.99*(h.size()-1))]};}
  // Tiny active-set least-squares system. Only the iso schedule is changed;
  // the scalar field and the cap values are immutable.
  for(int it=0;it<20;++it){Eigen::MatrixXd a=Eigen::MatrixXd::Identity(n+1,n+1)*1e-3;Eigen::VectorXd b=Eigen::VectorXd::Zero(n+1);a(n,n)=0;
    for(int k=0;k<n;++k){a(k,n)=a(n,k)=1;b[k]=1e-3;}b[n]=n;
    for(int k=0;k<n;++k)for(double h:{ranges[k].first,ranges[k].second}){
      Eigen::VectorXd row=Eigen::VectorXd::Zero(n);row[k]=h;
      const double value=row.dot(scale),epsilon=1e-5*s.Maximum;
      // This approximate initializer must leave an already feasible schedule
      // alone, including exact bounds and float-scale field-gradient noise.
      const double target=value<s.Minimum-epsilon?s.Minimum+margin:(value>s.Maximum+epsilon?s.Maximum-margin:value);
      if(value!=target){a.topLeftCorner(n,n)+=row*row.transpose();b.head(n)+=row*target;}
    }
    for(int k=1;k<n;++k){a(k,k)+=2;a(k-1,k-1)+=2;a(k,k-1)-=2;a(k-1,k)-=2;}
    Eigen::VectorXd next=a.ldlt().solve(b).head(n);if(!next.allFinite()){error="Iso schedule solve failed.";return false;}
    // Projection to positive intervals with their sum fixed. Bisection is
    // deterministic and keeps ordered iso values even for infeasible bounds.
    double low=-10,high=10;for(int j=0;j<60;++j){const double shift=(low+high)*.5;
      const double sum=(next.array()-shift).max(.9).min(1.1).sum();if(sum>n)low=shift;else high=shift;}
    next=(next.array()-(low+high)*.5).max(.9).min(1.1);const double change=(next-scale).norm();scale=next;if(change<1e-6)break;
  }
  levels.assign(n+1,0);for(int k=0;k<n;++k)levels[k+1]=levels[k]+scale[k]/n;levels.back()=1;return true;
}


namespace {
struct SurfaceSystem {
  int Generation=0;
  std::vector<int> Index;
  std::vector<Eigen::Matrix3d> Basis;
  std::vector<int> Dimensions;
  std::unique_ptr<Eigen::SimplicialLDLT<Sparse>> Solver;
  Eigen::MatrixXd FixedRhs;
};
constexpr double seedWeight=.001, coherenceWeight=.1, consensusPenalty=8;
bool surfaceSystem(const Layer &m,const Work &w,SurfaceSystem &system,std::string &error){
  zTsUtils::ScopedSlicingTiming methodTiming("Height.factorSurface",true);
  ++system.Generation;
  const int n=int(m.Positions.size());system.Index.assign(n,-1);system.Basis.assign(n,Eigen::Matrix3d::Zero());system.Dimensions.assign(n,0);int count=0;
  for(int i=0;i<n;++i)if(!w.fixed[i]){Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eig(w.Tangents[i]);
    for(int d=0;d<3;++d)if(eig.eigenvalues()[d]>.5)system.Basis[i].col(system.Dimensions[i]++)=eig.eigenvectors().col(d);
    if(system.Dimensions[i]){system.Index[i]=count;count+=system.Dimensions[i];}}
  system.Solver.reset();system.FixedRhs=Eigen::MatrixXd::Zero(count,1);if(count==0)return true;
  // The unknown is a displacement from the seed. Absolute surface fairing
  // contributes -B*X_seed; omitting this term smooths only the displacement
  // and can preserve or amplify a rough seed during height correction.
  Eigen::MatrixXd seed(n,3);for(int i=0;i<n;++i)seed.row(i)=(toEigenVec(w.original.Positions[i])-w.planeOrigin).transpose();
  const Eigen::MatrixXd fairingForce=-w.bending*seed;
  for(int i=0;i<n;++i)if(system.Index[i]>=0)for(int d=0;d<system.Dimensions[i];++d)
    system.FixedRhs(system.Index[i]+d,0)=system.Basis[i].col(d).dot(fairingForce.row(i).transpose());
  Sparse matrix=w.bending;for(int i=0;i<n;++i)matrix.coeffRef(i,i)+=(consensusPenalty+seedWeight+coherenceWeight)*w.mass[i];
  std::vector<Eigen::Triplet<double>> entries;
  for(int c=0;c<matrix.outerSize();++c)for(Sparse::InnerIterator it(matrix,c);it;++it){
    const int row=system.Index[it.row()],col=system.Index[it.col()];if(row<0)continue;
    for(int a=0;a<system.Dimensions[it.row()];++a){const V direction=system.Basis[it.row()].col(a);
      if(col>=0){for(int b=0;b<system.Dimensions[it.col()];++b){const double value=it.value()*direction.dot(system.Basis[it.col()].col(b));
        if(std::abs(value)>1e-20)entries.emplace_back(row+a,col+b,value);}}
      else system.FixedRhs(row+a,0)-=it.value()*direction.dot(toEigenVec(m.Positions[it.col()])-toEigenVec(w.original.Positions[it.col()]));}
  }
  Sparse reduced(count,count);reduced.setFromTriplets(entries.begin(),entries.end());
  system.Solver.reset(new Eigen::SimplicialLDLT<Sparse>());system.Solver->compute(reduced);
  if(system.Solver->info()!=Eigen::Success){error="Surface consensus factorization failed.";return false;}return true;
}
struct ProjectionTerm {int Layer,Vertex;V A=V::Zero(),Z=V::Zero(),Dual=V::Zero();double Weight=0;};
struct HeightProjection {std::vector<ProjectionTerm> Terms;double Bias=0,Lower=0,Upper=inf;int Layer=-1,Sample=0;};
struct ContactGuard {
  int Layer=0,FaceA=0,FaceB=0,VersionA=-1,VersionB=-1;
  V Normal=V::Zero();double Minimum=0;
  std::array<Eigen::MatrixXd,3> ResponseA,ResponseB;
};
using Positions=std::vector<Eigen::MatrixXd>;
Positions positions(const std::vector<Layer> &layers){Positions x(layers.size());
  for(size_t k=0;k<layers.size();++k){x[k].resize(layers[k].Positions.size(),3);for(int i=0;i<x[k].rows();++i)x[k].row(i)=toEigenVec(layers[k].Positions[i]);}return x;}

void projectContacts(Positions &x,const std::vector<Layer> &layers,const std::vector<SurfaceSystem> &systems,
                     std::vector<ContactGuard> &guards,double tolerance){
  zTsUtils::ScopedSlicingTiming methodTiming("Height.projectContacts");
  for(auto &guard:guards){const auto prepare=[&](int layer,int face,int &version,std::array<Eigen::MatrixXd,3> &responses,double sign){
      if(version==systems[layer].Generation)return;version=systems[layer].Generation;const auto &system=systems[layer];
      for(int j=0;j<3;++j){responses[j]=Eigen::MatrixXd::Zero(x[layer].rows(),3);if(!system.Solver)continue;
        const int vertex=layers[layer].Triangles[face][j],row=system.Index[vertex];if(row<0)continue;
        Eigen::VectorXd rhs=Eigen::VectorXd::Zero(system.FixedRhs.rows());for(int d=0;d<system.Dimensions[vertex];++d)
          rhs[row+d]=sign*system.Basis[vertex].col(d).dot(guard.Normal);
        const Eigen::VectorXd u=system.Solver->solve(rhs);
        for(int i=0;i<x[layer].rows();++i)if(system.Index[i]>=0)for(int d=0;d<system.Dimensions[i];++d)
          responses[j].row(i)+=u[system.Index[i]+d]*system.Basis[i].col(d).transpose();
      }};
    prepare(guard.Layer,guard.FaceA,guard.VersionA,guard.ResponseA,-1);
    prepare(guard.Layer+1,guard.FaceB,guard.VersionB,guard.ResponseB,1);
  }
  // Each response is Q^-1 times a contact gradient. This is the minimum
  // fairing-metric correction, spread smoothly across each independent layer,
  // rather than a point displacement or a permanently frozen contact patch.
  for(int pass=0;pass<12;++pass){double violation=0;
    for(const auto &guard:guards){const int layer=guard.Layer;const auto &a=layers[layer].Triangles[guard.FaceA],&b=layers[layer+1].Triangles[guard.FaceB];
      for(int i=0;i<3;++i)for(int j=0;j<3;++j){const double gap=guard.Normal.dot((x[layer+1].row(b[j])-x[layer].row(a[i])).transpose()),missing=guard.Minimum-gap;
        violation=(std::max)(violation,missing);if(missing<=tolerance)continue;
        const double den=guard.Normal.dot((guard.ResponseB[j].row(b[j])-guard.ResponseA[i].row(a[i])).transpose());
        if(den>1e-20){const double step=missing/den;x[layer]+=step*guard.ResponseA[i];x[layer+1]+=step*guard.ResponseB[j];}}}
    if(violation<=tolerance)break;
  }
}
V tangentCoefficient(const V &a,int layer,int vertex,const std::vector<Work> &work,const SideConstraints &side){
  return work[layer].Tangents[vertex]*a;
}

using Matrix3=Eigen::Matrix3d;
Matrix3 crossMatrix(const V &v){Matrix3 m;m<<0,-v.z(),v.y(),v.z(),0,-v.x(),-v.y(),v.x(),0;return m;}
// Sparse one-ring differentials of the EXISTING area-weighted/recovered normal.
// They are evaluated lazily, only for active height supports, and frozen for
// a consensus batch. No cubic-control derivative arrays or Newton matrix.
class NormalStencils {
  const Layer &layer_;
  std::vector<V> sums_;
  std::vector<double> scale_;
  std::vector<std::vector<int>> faces_,neighbours_;
  std::map<int,std::vector<std::pair<int,Matrix3>>> cache_;
  double coordinateScale_=1;
  bool reliable(int i)const{return sums_[i].norm()>(std::max)(1e-6*scale_[i],16*std::numeric_limits<float>::epsilon()*coordinateScale_*std::sqrt(scale_[i]));}
public:
  explicit NormalStencils(const Layer &layer):layer_(layer),sums_(layer.Positions.size(),V::Zero()),scale_(layer.Positions.size(),0),faces_(layer.Positions.size()),neighbours_(layer.Positions.size()){
    for(const auto &p:layer.Positions)coordinateScale_=(std::max)(coordinateScale_,toEigenVec(p).cwiseAbs().maxCoeff());
    for(int f=0;f<int(layer.Triangles.size());++f){const auto &t=layer.Triangles[f];V p[3];for(int j=0;j<3;++j)p[j]=toEigenVec(layer.Positions[t[j]]);
      const V n=(p[1]-p[0]).cross(p[2]-p[0]);double longest=0;for(int j=0;j<3;++j)longest=(std::max)(longest,(p[j]-p[(j+1)%3]).squaredNorm());
      for(int j=0;j<3;++j){sums_[t[j]]+=n;scale_[t[j]]=(std::max)(scale_[t[j]],longest);faces_[t[j]].push_back(f);neighbours_[t[j]].push_back(t[(j+1)%3]);neighbours_[t[j]].push_back(t[(j+2)%3]);}}
  }
  const std::vector<std::pair<int,Matrix3>> &get(int id){auto found=cache_.find(id);if(found!=cache_.end())return found->second;
    std::vector<int> supports;if(reliable(id))supports.push_back(id);else{std::set<int> visited{id};std::vector<int> ring{id};
      while(!ring.empty()&&supports.empty()){std::vector<int> next;for(int i:ring)for(int j:neighbours_[i])if(visited.insert(j).second){if(reliable(j))supports.push_back(j);else next.push_back(j);}ring.swap(next);}}
    std::map<int,Matrix3> derivative;V sum=V::Zero();
    for(int anchor:supports){sum+=sums_[anchor];for(int f:faces_[anchor]){const auto &tri=layer_.Triangles[f];V p[3];for(int j=0;j<3;++j)p[j]=toEigenVec(layer_.Positions[tri[j]]);
      for(int j=0;j<3;++j){const Matrix3 d=crossMatrix(p[(j+2)%3]-p[(j+1)%3]);auto entry=derivative.find(tri[j]);if(entry==derivative.end())derivative.emplace(tri[j],d);else entry->second+=d;}}}
    auto &result=cache_[id];if(sum.norm()>0){const V n=sum.normalized();const Matrix3 normalize=(Matrix3::Identity()-n*n.transpose())/sum.norm();
      for(const auto &entry:derivative)result.emplace_back(entry.first,normalize*entry.second);}return result;
  }
};
V closestWeightDerivative(const Face &face,const V &weights,const V &p,const V &dp,const std::array<V,3> &d){
  const V a=face.p[1]-face.p[0],b=face.p[2]-face.p[0],q=face.p[0]+weights[1]*a+weights[2]*b;
  const V material=weights[0]*d[0]+weights[1]*d[1]+weights[2]*d[2],residual=p-q;
  int active[3],count=0;for(int j=0;j<3;++j)if(weights[j]>1e-8)active[count++]=j;
  if(count==3){const double aa=a.squaredNorm(),ab=a.dot(b),bb=b.squaredNorm(),det=aa*bb-ab*ab;
    if(det<=1e-18*aa*bb)return V::Zero();const double x=(d[1]-d[0]).dot(residual)+a.dot(dp-material),y=(d[2]-d[0]).dot(residual)+b.dot(dp-material);
    const double dv=(bb*x-ab*y)/det,dw=(aa*y-ab*x)/det;return V(-dv-dw,dv,dw);}
  if(count==2){const int i=active[0],j=active[1];const V edge=face.p[j]-face.p[i];if(edge.squaredNorm()==0)return V::Zero();
    const double dt=((d[j]-d[i]).dot(residual)+edge.dot(dp-material))/edge.squaredNorm();V derivative=V::Zero();derivative[i]=-dt;derivative[j]=dt;return derivative;}
  return V::Zero();
}
std::vector<HeightProjection> heightProjections(const std::vector<Layer> &layers,const HeightLinks &links,
    const HeightSamples &samples,const std::vector<Work> &work,const SideConstraints &side,
    std::vector<Eigen::VectorXd> &degrees,const PrintHeightOptimizationSettings &settings,bool orientation,
    std::vector<std::vector<double>> &multipliers,const std::vector<ContactGuard> &contacts,HeightReport &report){
  zTsUtils::ScopedSlicingTiming methodTiming("Height.buildProjections");
  report.UncontrollableSamples=0;
  std::vector<HeightProjection> rows;degrees.resize(layers.size());
  std::vector<std::unique_ptr<NormalStencils>> normalStencils(layers.size());
  std::vector<zVectorArray> capNormals(layers.size());
  for(size_t k=0;k<layers.size();++k)degrees[k]=Eigen::VectorXd::Zero(layers[k].Positions.size());
  for(int k=1;k<int(samples.size());++k)for(int sampleId=0;sampleId<int(samples[k].size());++sampleId){const auto &sample=samples[k][sampleId];
    const double height=sampleHeight(layers[k],sample),margin=.05*(settings.Maximum-settings.Minimum);
    if(height>settings.Minimum+margin&&height<settings.Maximum-margin){multipliers[k][sampleId]=0;continue;}
    std::map<std::pair<int,int>,V> coefficients;
    const auto add=[&](int layer,int id,const V &a){if(layer<=0||layer+1>=int(layers.size())||work[layer].fixed[id])return;
      const auto key=std::make_pair(layer,id);auto found=coefficients.find(key);
      if(found==coefficients.end())coefficients.emplace(key,a);else found->second+=a;};

    const auto normalContribution=[&](int layer,int vertex,const V &gradient){if(!orientation||layer<=0||layer+1>=int(layers.size()))return;
      if(!normalStencils[layer])normalStencils[layer].reset(new NormalStencils(layers[layer]));
      for(const auto &d:normalStencils[layer]->get(vertex))add(layer,d.first,d.second.transpose()*gradient);
    };
    if(sample.Exact){const auto &hit=sample.Hit;const double den=hit.PlaneNormal.dot(hit.Direction);
      if(std::abs(den)>1e-10){const V a=hit.PlaneNormal/den;
        for(int j=0;j<3;++j)add(k,sample.Ids[j],-sample.Weights[j]*a);
        for(int j=0;j<3;++j)add(hit.LayerId,hit.Ids[j],hit.Weights[j]*a);
        // N = normalize(sum alpha_j N_j). Differentiate that normalization
        // before the existing recovered-normal one-ring stencils.
        const V normalGradient=(hit.Height/sample.NormalLength)*
          (a-sample.Normal*sample.Normal.dot(a));
        for(int j=0;j<3;++j)if(sample.Weights[j])
          normalContribution(k,sample.Ids[j],sample.Weights[j]*normalGradient);
      }
    }else for(int j=0;j<3;++j){const double alpha=sample.Weights[j];if(alpha==0)continue;const int id=sample.Ids[j];const auto &hits=links[k][id];
      for(const auto &hit:hits){if(hit.LayerId!=(k==0?1:k-1))continue;
        const double den=hit.PlaneNormal.dot(hit.Direction);if(std::abs(den)<1e-10)continue;
        const V a=alpha*hit.PlaneNormal/den;add(k,id,-a);
        for(int v=0;v<3;++v)add(hit.LayerId,hit.Ids[v],hit.Weights[v]*a);
        if(k>0&&k+1<int(layers.size())){const double sign=hit.Direction.dot(toEigenVec(layers[k].Normals[id]))>=0?1:-1;
          normalContribution(k,id,-hit.Height*sign*a);}
        if(hit.Fallback){const V p=toEigenVec(layers[k].Positions[id]),residual=hit.Position-p-hit.Height*hit.Direction;
          const V normalGradient=alpha*(residual-hit.PlaneNormal*hit.PlaneNormal.dot(residual))/(den*hit.NormalLength);
          for(int v=0;v<3;++v)normalContribution(hit.LayerId,hit.Ids[v],hit.Weights[v]*normalGradient);
          Face face;V sensitivity;if(layers[hit.LayerId].Normals.empty()&&capNormals[hit.LayerId].empty())capNormals[hit.LayerId]=normals(layers[hit.LayerId]);
          const auto &adjacentNormals=layers[hit.LayerId].Normals.empty()?capNormals[hit.LayerId]:layers[hit.LayerId].Normals;
          for(int v=0;v<3;++v){face.p[v]=toEigenVec(layers[hit.LayerId].Positions[hit.Ids[v]]);
            sensitivity[v]=alpha*hit.PlaneNormal.dot(face.p[v])/den+normalGradient.dot(toEigenVec(adjacentNormals[hit.Ids[v]]));}
          V sourceGradient=V::Zero();std::array<V,3> zero{V::Zero(),V::Zero(),V::Zero()};
          for(int c=0;c<3;++c)sourceGradient[c]=sensitivity.dot(closestWeightDerivative(face,hit.Weights,p,V::Unit(c),zero));add(k,id,sourceGradient);
          for(int v=0;v<3;++v){V gradient;for(int c=0;c<3;++c){auto d=zero;d[v]=V::Unit(c);gradient[c]=sensitivity.dot(closestWeightDerivative(face,hit.Weights,p,V::Zero(),d));}add(hit.LayerId,hit.Ids[v],gradient);}
        }
      }}
    HeightProjection row;row.Layer=k;row.Sample=sampleId;row.Bias=sampleHeight(layers[k],sample);
    const double guard=(std::min)(.02*(settings.Maximum-settings.Minimum),(std::max)(1e-5*settings.Minimum,2*settings.Tolerance));
    row.Lower=settings.Minimum+guard;row.Upper=settings.Maximum-guard;
    for(const auto &entry:coefficients){const int l=entry.first.first,i=entry.first.second;const V a=tangentCoefficient(entry.second,l,i,work,side);
      if(a.squaredNorm()<1e-20)continue;ProjectionTerm term;term.Layer=l;term.Vertex=i;term.A=a;
      row.Bias-=a.dot(toEigenVec(layers[l].Positions[i]));row.Terms.push_back(term);}

    double available=0;for(const auto &term:row.Terms)available+=2*work[term.Layer].limit*term.A.norm();
    if(!row.Terms.empty()&&available>=bandViolation(height,settings)){
      for(const auto &term:row.Terms)degrees[term.Layer][term.Vertex]+=1;rows.push_back(std::move(row));}
    else if(sample.Primary&&bandViolation(height,settings)>settings.Tolerance)++report.UncontrollableSamples;

  }
  // Both directions still supply contact guards. Constraining the deposited
  // previous-layer TCP gap does not by itself prevent contact with the next.
  const double gapFloor=.1*settings.Minimum;
  for(int k=0;k<int(layers.size());++k)for(int i=0;i<int(links[k].size());++i)for(const auto &hit:links[k][i]){
    if(hit.Height>3*gapFloor)continue;const double den=hit.PlaneNormal.dot(hit.Direction);if(std::abs(den)<1e-10)continue;
    HeightProjection row;row.Bias=hit.Height;row.Lower=gapFloor;const V a=hit.PlaneNormal/den;
    const auto add=[&](int layer,int vertex,const V &gradient){if(layer<=0||layer+1>=int(layers.size())||work[layer].fixed[vertex])return;
      ProjectionTerm t;t.Layer=layer;t.Vertex=vertex;t.A=work[layer].Tangents[vertex]*gradient;if(t.A.squaredNorm()<1e-20)return;
      row.Bias-=t.A.dot(toEigenVec(layers[layer].Positions[vertex]));row.Terms.push_back(t);};
    add(k,i,-a);for(int j=0;j<3;++j)add(hit.LayerId,hit.Ids[j],hit.Weights[j]*a);
    if(!row.Terms.empty()){for(const auto &t:row.Terms)degrees[t.Layer][t.Vertex]+=1;rows.push_back(std::move(row));}
  }
  for(const auto &guard:contacts)for(int a=0;a<3;++a)for(int b=0;b<3;++b){HeightProjection row;row.Lower=guard.Minimum;
    const auto add=[&](int layer,int face,const V &weights,double sign){const auto &ids=layers[layer].Triangles[face];
      for(int j=0;j<3;++j){const V gradient=sign*weights[j]*guard.Normal;const int vertex=ids[j];
        if(layer<=0||layer+1>=int(layers.size())||work[layer].fixed[vertex]){row.Bias+=gradient.dot(toEigenVec(layers[layer].Positions[vertex]));continue;}
        ProjectionTerm term;term.Layer=layer;term.Vertex=vertex;term.A=work[layer].Tangents[vertex]*gradient;
        row.Bias+=(gradient-term.A).dot(toEigenVec(layers[layer].Positions[vertex]));if(term.A.squaredNorm()>1e-20)row.Terms.push_back(term);}};
    add(guard.Layer,guard.FaceA,V::Unit(a),-1);add(guard.Layer+1,guard.FaceB,V::Unit(b),1);
    if(!row.Terms.empty()){for(const auto &t:row.Terms)degrees[t.Layer][t.Vertex]+=1;rows.push_back(std::move(row));}
  }
  // Incidence-normalized replica masses sum to M at every participating vertex.
  // The global matrix therefore stays M + bending, independent of normals,
  // triangle-hit changes and which height constraints are currently active.
  for(auto &row:rows)for(auto &term:row.Terms){term.Weight=consensusPenalty*work[term.Layer].mass[term.Vertex]/degrees[term.Layer][term.Vertex];
    if(row.Layer>=0)term.Dual=multipliers[row.Layer][row.Sample]*term.A/term.Weight;}
  return rows;
}
double bendingEnergy(const std::vector<Layer> &layers,const std::vector<Work> &work){double value=0,area=0;
  for(size_t k=1;k+1<layers.size();++k){Eigen::MatrixXd x(layers[k].Positions.size(),3);
    for(int i=0;i<x.rows();++i)x.row(i)=(toEigenVec(layers[k].Positions[i])-work[k].planeOrigin).transpose();
    value+=(x.array()*(work[k].bending*x).array()).sum();area+=work[k].mass.sum();}
  return value/(std::max)(area,1e-30);
}
}

bool measureTargetHeights(const std::vector<Layer> &layers,
    const std::vector<zPointArray> &targets,const std::vector<zVectorArray> &targetNormals,
    std::vector<zFloatArray> &outHeights,std::string &error,bool parallel) {
  zTsUtils::ScopedSlicingTiming timing("Height.measureTargetHeights");outHeights.clear();error.clear();
  if(layers.size()<2||targets.size()+1!=layers.size()||targetNormals.size()!=targets.size()){
    error="TCP targets/normals must match sections S1..Slast.";return false;}
  Box box;for(const auto &layer:layers){if(!meshValid(layer,error))return false;for(const auto &p:layer.Positions)box.add(toEigenVec(p));}
  const double maximum=2*(box.hi-box.lo).norm(),tol=(std::max)(1e-9,2e-7*(box.hi-box.lo).norm());
  std::vector<std::unique_ptr<Tree>> trees(layers.size());zTsDetail::ParallelErrors preparation(layers.size());
  #pragma omp parallel for schedule(dynamic) if(parallel)
  for(int k=0;k<int(layers.size());++k){try{trees[k].reset(new Tree(layers[k].Positions,layers[k].Triangles));}catch(...){preparation.capture(k);}}
  preparation.rethrow();std::vector<zFloatArray> measured(targets.size());std::vector<std::string> errors(targets.size());
  zTsDetail::ParallelErrors queries(targets.size());
  #pragma omp parallel for schedule(dynamic) if(parallel)
  for(int p=0;p<int(targets.size());++p){try{
    if(targetNormals[p].size()!=targets[p].size()){errors[p]="TCP normal count differs from its targets.";continue;}
    measured[p].reserve(targets[p].size());int sourceHint=-1,previousHint=-1;
    for(size_t i=0;i<targets[p].size();++i){const V position=toEigenVec(targets[p][i]);V normal=toEigenVec(targetNormals[p][i]);
      if(!position.allFinite()||!normal.allFinite()||normal.norm()<1e-12){errors[p]="Invalid TCP position or normal.";break;}normal.normalize();
      int source=-1,previous=-1;V nearest;
      if(trees[p+1]->nearest(position,source,nearest,sourceHint)>10*tol){errors[p]="TCP target is off its indexed section.";break;}sourceHint=source;
      if(trees[p]->nearest(position,previous,nearest,previousHint)<=tol){errors[p]="TCP target touches its previous section.";break;}previousHint=previous;
      int hit=-1;const double height=trees[p]->ray(position,-normal,maximum,tol,&hit,previousHint);
      if(!std::isfinite(height)||hit<0){errors[p]="TCP normal ray has no previous-layer support on printed layer "+std::to_string(p)+", target "+std::to_string(i)+".";break;}
      previousHint=hit;measured[p].push_back(float(height));
    }
  }catch(...){queries.capture(p);}}
  queries.rethrow();for(const auto &e:errors)if(!e.empty()){error=e;return false;}
  outHeights.swap(measured);return true;
}

bool optimizeHeightDomain(const Reference &reference,std::vector<Layer> &layers,
    const PrintHeightOptimizationSettings &s,const SampleVolume &sample,bool parallel,
    HeightReport &report,std::string &error,const std::vector<zPointArray> *printTargets,double fairingStrength) {
  zTsUtils::ScopedSlicingTiming methodTiming("Height.optimizeHeightDomain");
  report=HeightReport();error.clear();const auto start=std::chrono::steady_clock::now();
  const auto finish=[&](){report.Milliseconds=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-start).count();};
  if(!s.Enabled||!s.OptimizeSurfaces)return true;
  if(!validHeightSettings(s,error)||!std::isfinite(fairingStrength)||fairingStrength<=0)return false;
  if(layers.size()<3||!sample||reference.Triangles.size()!=reference.WidthIds.size()){error="Missing surface optimizer geometry or volume query.";return false;}
  Layer ref;ref.Positions=reference.Positions;ref.Triangles=reference.Triangles;if(!meshValid(ref,error))return false;
  for(const auto &m:layers)if(!meshValid(m,error))return false;
  report.InitialSections=report.FinalSections=int(layers.size());
  Box box;for(const auto &p:reference.Positions)box.add(toEigenVec(p));const double tol=(std::max)(1e-9,2e-7*(box.hi-box.lo).norm());
  SideConstraints side(reference,tol,true);if(!side.build(error))return false;
  HeightLinks links;Report ignored;std::vector<std::unique_ptr<Tree>> trees;
  if(!refreshHeights(layers,tol,ignored,error,parallel,&links,nullptr,nullptr,&trees))return false;
  const auto directedGaps=[&](){report.MinimumDirectedGap=inf;report.MaximumDirectedGap=0;report.UnavailableAuxiliaryGaps=0;
    for(size_t k=0;k<links.size();++k)for(const auto &vertex:links[k]){
      report.UnavailableAuxiliaryGaps+=(k>0&&k+1<links.size()?2:1)-int(vertex.size());
      for(const auto &hit:vertex){
      report.MinimumDirectedGap=(std::min)(report.MinimumDirectedGap,hit.Height);
      report.MaximumDirectedGap=(std::max)(report.MaximumDirectedGap,hit.Height);}}};
  directedGaps();
  HeightSamples samples(layers.size());
  if(printTargets&&printTargets->size()+1!=layers.size()){error="Height targets must match printed TCP sections S1..Slast.";return false;}
  for(size_t k=1;k<layers.size();++k){if(printTargets){
      for(const auto &p:(*printTargets)[k-1]){int tri;V q;const double d=trees[k]->nearest(toEigenVec(p),tri,q);if(tri<0||d>10*tol){error="Print-height target is off its TCP section.";return false;}
        V w=barycentric(q,trees[k]->face(tri)).cwiseMax(0);w/=w.sum();samples[k].push_back({layers[k].Triangles[tri],w});samples[k].back().Exact=true;
      }
    }else {
      // Unprinted trim-boundary slivers must not dominate the initial stack.
      // Cover the surface interior first; the authoritative final footprint
      // is checked after mapping and may add barycentric feedback constraints.
      std::set<int> boundary;for(const auto &e:layers[k].BoundaryWidths){boundary.insert(int(e.first>>32));boundary.insert(int(e.first&0xffffffffu));}
      for(int i=0;i<int(layers[k].Positions.size());++i)if(!boundary.count(i))samples[k].push_back({{i,i,i},V(1,0,0)});
      if(samples[k].empty())for(const auto &tri:layers[k].Triangles)samples[k].push_back({tri,V::Constant(1.0/3)});
    }}
  size_t primaryCount=0;for(const auto &layer:samples)for(const auto &point:layer)primaryCount+=point.Primary;
  if(primaryCount==0){error="No material print-curve height supports were supplied.";return false;}
  const double maximumRay=2*(box.hi-box.lo).norm();
  if(!refreshExactSamples(layers,samples,trees,maximumRay,tol,error,parallel))return false;
  double score=heightScore(layers,samples,s,report);report.InitialMinimum=report.Minimum;report.InitialMaximum=report.Maximum;
  if(report.Converged){finish();return true;}
  std::vector<Work> work(layers.size());std::vector<SurfaceSystem> systems(layers.size());std::vector<std::string> errors(layers.size());
  std::vector<std::unique_ptr<Tree>> trialTrees(layers.size());std::vector<char> dirtyTrees(layers.size(),0);
  std::vector<int> geometryPriority(layers.size(),0);
  std::vector<int> initialSlivers(layers.size(),0);
  Settings fair;fair.Iterations=1;fair.Strength=fairingStrength;fair.Parallel=parallel;
  zTsDetail::ParallelErrors setup(layers.size());
  zTsUtils::ScopedSlicingTiming setupTiming("Height.prepareSystems");
  #pragma omp parallel for schedule(dynamic) if(parallel)
  for(int k=1;k<int(layers.size())-1;++k){try{Report local;
    if(!prepare(layers[k],side,s.Maximum,fair,work[k],local,errors[k],false))continue;
    initialSlivers[k]=local.SliverConstrainedVertices;
    // A thin triangle does not by itself require an immobile patch. The
    // scalar normal system moves its supports coherently; hard trial guards
    // still reject collapse, reversal, contact and missing transported UV.
    work[k].fixed.assign(layers[k].Positions.size(),false);
    work[k].limit=s.MaxDisplacement>0?s.MaxDisplacement:2*s.Maximum;
    // Height correction needs a smooth scalar offset, not free tangential
    // mesh motion. Project each seed normal into the exact side/crease tangent
    // space and retain that one direction throughout the correction solve.
    for(int i=0;i<int(layers[k].Positions.size());++i){V direction=work[k].Tangents[i]*toEigenVec(layers[k].Normals[i]);
      if(direction.norm()>1e-8){direction.normalize();work[k].Tangents[i]=direction*direction.transpose();}
      else {work[k].Tangents[i].setZero();work[k].fixed[i]=true;}}
    std::vector<double> masses(work[k].mass.data(),work[k].mass.data()+work[k].mass.size());
    std::nth_element(masses.begin(),masses.begin()+masses.size()/2,masses.end());const double floor=.1*masses[masses.size()/2];
    for(int i=0;i<work[k].mass.size();++i)work[k].mass[i]=(std::max)(work[k].mass[i],floor);

    std::set<int> blocked;if(!validTrial(layers[k],work[k],side,sample,tol,errors[k],blocked))continue;
    surfaceSystem(layers[k],work[k],systems[k],errors[k]);
  }catch(...){setup.capture(k);}}
  setup.rethrow();for(const auto &e:errors)if(!e.empty()){error=e;return false;}
  setupTiming.stop();
  for(size_t k=1;k+1<layers.size();++k){report.Factorizations+=systems[k].Solver?1:0;
    report.InitialSliverVertices+=initialSlivers[k];
    for(int dimensions:systems[k].Dimensions)report.DegreesOfFreedom+=dimensions;}
  report.InitialBendingEnergy=bendingEnergy(layers,work);
  const auto seed=positions(layers);int stalled=0;
  std::vector<ContactGuard> contactGuards;std::set<std::tuple<int,int,int>> contactKeys;
  std::vector<std::vector<double>> multipliers(layers.size());for(size_t k=0;k<samples.size();++k)multipliers[k].assign(samples[k].size(),0);
  // Parallel planar stacks have an exact constant-offset subspace. Preserve
  // it explicitly instead of amplifying float-scale ray-hit changes into
  // spatial variation during successive correspondence refreshes.
  const V planeAxis=toEigenVec(normals(layers.front()).front());bool parallelPlanes=true;
  for(const auto &layer:layers){const double level=planeAxis.dot(toEigenVec(layer.Positions.front()));
    for(const auto &p:layer.Positions)if(std::abs(planeAxis.dot(toEigenVec(p))-level)>tol){parallelPlanes=false;break;}}
  for(int outer=0;outer<s.Iterations;++outer){++report.OuterIterations;
    const auto outerStart=std::chrono::steady_clock::now();
    const size_t initialContacts=contactGuards.size();
    std::vector<Eigen::VectorXd> degrees;auto rows=heightProjections(layers,links,samples,work,side,degrees,s,!parallelPlanes,multipliers,contactGuards,report);if(rows.empty())break;
    std::vector<char> active(layers.size(),0);
    for(size_t k=1;k+1<layers.size();++k)active[k]=degrees[k].sum()>0;
    if(outer==0&&printTargets&&std::getenv("ZSPACE_CHECK_HEIGHT_GRADIENT")){
      std::vector<int> check;for(int r=0;r<int(rows.size());++r)if(rows[r].Layer>=0)check.push_back(r);
      std::sort(check.begin(),check.end(),[&](int a,int b){return bandViolation(sampleHeight(layers[rows[a].Layer],samples[rows[a].Layer][rows[a].Sample]),s)>
        bandViolation(sampleHeight(layers[rows[b].Layer],samples[rows[b].Layer][rows[b].Sample]),s);});
      for(int r=0;r<(std::min)(3,int(check.size()));++r){const auto &row=rows[check[r]];
        const auto term=std::max_element(row.Terms.begin(),row.Terms.end(),[](const ProjectionTerm &a,const ProjectionTerm &b){return a.A.squaredNorm()<b.A.squaredNorm();});
        const V direction=term->A.normalized();const double epsilon=.001*s.Minimum;auto plus=layers,minus=layers;
        plus[term->Layer].Positions[term->Vertex]=point(toEigenVec(layers[term->Layer].Positions[term->Vertex])+epsilon*direction);
        minus[term->Layer].Positions[term->Vertex]=point(toEigenVec(layers[term->Layer].Positions[term->Vertex])-epsilon*direction);
        std::vector<char> mask(layers.size(),0);mask[row.Layer]=1;Report temp;std::string message;
        std::vector<std::unique_ptr<Tree>> plusTrees,minusTrees;auto plusSamples=samples,minusSamples=samples;
        if(refreshHeights(plus,tol,temp,message,false,nullptr,&mask,nullptr,&plusTrees)&&
            refreshHeights(minus,tol,temp,message,false,nullptr,&mask,nullptr,&minusTrees)&&
            refreshExactSamples(plus,plusSamples,plusTrees,maximumRay,tol,message,false,&mask)&&
            refreshExactSamples(minus,minusSamples,minusTrees,maximumRay,tol,message,false,&mask)){
          const double finite=(sampleHeight(plus[row.Layer],plusSamples[row.Layer][row.Sample])-sampleHeight(minus[row.Layer],minusSamples[row.Layer][row.Sample]))/(2*epsilon);
          std::ostringstream log;log<<"[surfaceGradientCheck] layer="<<row.Layer<<" predicted="<<term->A.norm()<<" finiteDifference="<<finite;zTsUtils::printDebug(log.str());}
      }
    }
    Positions x=positions(layers);
    // Neighbour consensus is lagged once per outer batch, just like the hit
    // supports. Its right-hand side does not need to be assembled again in
    // every inexpensive local/global projection.
    Positions baseRhs(layers.size());
    for(size_t k=1;k+1<layers.size();++k){auto &system=systems[k];baseRhs[k]=system.FixedRhs;
      for(int i=0;i<x[k].rows();++i){const int row=system.Index[i];if(row<0)continue;V neighbour=V::Zero();int count=0;
        for(const auto &hit:links[k][i]){V delta=V::Zero();for(int j=0;j<3;++j)
          delta+=hit.Weights[j]*(x[hit.LayerId].row(hit.Ids[j])-seed[hit.LayerId].row(hit.Ids[j])).transpose();neighbour+=delta;++count;}
        if(count)neighbour/=count;const V target=work[k].mass[i]*(coherenceWeight*neighbour-consensusPenalty*seed[k].row(i).transpose());
        for(int d=0;d<system.Dimensions[i];++d)baseRhs[k](row+d,0)+=system.Basis[i].col(d).dot(target);}}
    // Consensus ADMM: local slab projections, independent cached surface solves,
    // then dual updates. Sparse support/orientation coefficients are frozen;
    // there is no global Jacobian matrix or Newton factorization.
    // Coefficients and hit supports are frozen only for a bounded batch.
    // Persistent duals carry progress across refreshes; contact appearance
    // must not silently multiply every later batch's sparse-solve budget.
    constexpr int innerBudget=32;
    zTsUtils::ScopedSlicingTiming consensusTiming("Height.consensus");
    for(int inner=0;inner<innerBudget;++inner){
      ++report.InnerIterations;
      #pragma omp parallel for schedule(static) if(parallel)
      for(int r=0;r<int(rows.size());++r){auto &row=rows[r];double value=row.Bias,den=0;
        for(const auto &t:row.Terms){value+=t.A.dot(x[t.Layer].row(t.Vertex).transpose()+t.Dual);den+=t.A.squaredNorm()/t.Weight;}
        const double target=std::clamp(value,row.Lower,row.Upper),scale=(target-value)/den;
        for(auto &t:row.Terms){const V previous=x[t.Layer].row(t.Vertex).transpose();
          const V projected=previous+t.Dual+scale*t.A/t.Weight;
          t.Z=1.6*projected-.6*previous;}
      }
      Positions rhs(layers.size());for(size_t k=1;k+1<layers.size();++k)rhs[k]=Eigen::MatrixXd::Zero(x[k].rows(),3);
      for(const auto &row:rows)for(const auto &t:row.Terms)rhs[t.Layer].row(t.Vertex)+=t.Weight*(t.Z-t.Dual).transpose();
      const auto previous=x;zTsDetail::ParallelErrors solves(layers.size());
      #pragma omp parallel for schedule(dynamic) if(parallel)
      for(int k=1;k<int(layers.size())-1;++k){try{auto &system=systems[k];
        // A height batch has no authority to refair a layer without active
        // height/contact supports. Keep it exact, and avoid its sparse solves.
        if(!active[k]||!system.Solver)continue;
        Eigen::MatrixXd b=baseRhs[k];
        for(int i=0;i<x[k].rows();++i){const int row=system.Index[i];if(row<0)continue;const double mass=work[k].mass[i];
          V target=rhs[k].row(i).transpose();if(degrees[k][i]==0)target+=consensusPenalty*mass*previous[k].row(i).transpose();
          for(int d=0;d<system.Dimensions[i];++d)b(row+d,0)+=system.Basis[i].col(d).dot(target);
        }
        const Eigen::MatrixXd solution=system.Solver->solve(b);if(!solution.allFinite())throw std::runtime_error("Surface consensus solve returned nonfinite positions.");
        for(int i=0;i<x[k].rows();++i)if(system.Index[i]>=0){V p=seed[k].row(i).transpose();
          for(int d=0;d<system.Dimensions[i];++d)p+=solution(system.Index[i]+d,0)*system.Basis[i].col(d);x[k].row(i)=p.transpose();}
      }catch(...){solves.capture(k);}}
      solves.rethrow();
      if(parallelPlanes)for(size_t k=1;k+1<layers.size();++k){double offset=0;
        for(int i=0;i<x[k].rows();++i)offset+=planeAxis.dot(x[k].row(i).transpose()-toEigenVec(layers[k].Positions[i]));offset/=x[k].rows();
        if(std::find(work[k].fixed.begin(),work[k].fixed.end(),true)!=work[k].fixed.end())offset=0;
        for(int i=0;i<x[k].rows();++i)x[k].row(i)=(toEigenVec(layers[k].Positions[i])+offset*planeAxis).transpose();}
      for(auto &row:rows)for(auto &t:row.Terms)t.Dual+=x[t.Layer].row(t.Vertex).transpose()-t.Z;
      if(inner>=7&&inner%4==3){double violation=0,change=0;
        for(const auto &row:rows){double value=row.Bias;for(const auto &t:row.Terms)value+=t.A.dot(x[t.Layer].row(t.Vertex).transpose());
          violation=(std::max)({violation,row.Lower-value,value-row.Upper});}
        for(size_t k=1;k+1<layers.size();++k)change=(std::max)(change,(x[k]-previous[k]).rowwise().norm().maxCoeff());
        if(violation<.25*s.Tolerance+1e-6*s.Minimum&&change<1e-3*s.Minimum)break;}
    }
    consensusTiming.stop();
    projectContacts(x,layers,systems,contactGuards,1e-6*s.Minimum);
    // Trust limits belong to each complete layer. One ill-conditioned local
    // proposal must not throttle every unrelated layer in the stack.
    std::vector<double> scales(layers.size(),1);
    for(size_t k=1;k+1<layers.size();++k)for(int i=0;i<x[k].rows();++i){const V d=x[k].row(i).transpose()-toEigenVec(layers[k].Positions[i]);
      if(d.norm()>0)scales[k]=(std::min)(scales[k],.5*s.Minimum/d.norm());}
    const double minimumTrust=*std::min_element(scales.begin(),scales.end());
    bool accepted=false;double acceptedLine=0;std::string rejected;std::vector<std::set<int>> blocked(layers.size()),blockedUnion(layers.size());HeightQueryFailure failedQuery;
    std::vector<double> lineSteps(layers.size(),1),validatedSteps(layers.size(),-1);
    const auto reduceLayer=[&](int k){if(k<=0||k+1>=int(layers.size())||lineSteps[k]==0)return false;
      lineSteps[k]*=.5;if(lineSteps[k]<1.0/128)lineSteps[k]=0;return true;};
    for(int search=0;search<16;++search){++report.LineSearchTrials;
      for(auto &ids:blockedUnion)ids.clear();failedQuery=HeightQueryFailure();
      zTsUtils::ScopedSlicingTiming trialTiming("Height.lineSearchTrial");
      auto trial=layers;std::vector<char> valid(layers.size(),1),moved(layers.size(),0);zTsDetail::ParallelErrors geometry(layers.size());
      std::vector<std::vector<char>> changedVertices(layers.size()),changedFaces(layers.size());
      #pragma omp parallel for schedule(dynamic) if(parallel)
      for(int k=1;k<int(layers.size())-1;++k){try{
        blocked[k].clear();changedVertices[k].assign(layers[k].Positions.size(),0);changedFaces[k].assign(layers[k].Triangles.size(),0);
        for(int i=0;i<x[k].rows();++i){if(work[k].fixed[i])continue;const V old=toEigenVec(layers[k].Positions[i]),origin=toEigenVec(work[k].original.Positions[i]);
          const V p=side.project(old+lineSteps[k]*scales[k]*(x[k].row(i).transpose()-old),origin,work[k].allowed[i]);
          if((p-origin).norm()>work[k].limit+tol){valid[k]=0;errors[k]="Surface displacement limit reached.";break;}
          if((p-old).norm()<=tol)continue;
          trial[k].Positions[i]=point(p);zPoint uv;
          if(trial[k].UV.size()==trial[k].Positions.size()){
            if(!sample(trial[k].Positions[i],&uv)){valid[k]=0;errors[k]="Transported surface UV is unavailable.";blocked[k].insert(i);break;}
            trial[k].UV[i]=uv;
          }
          moved[k]=1;changedVertices[k][i]=1;
        }

        if(valid[k]&&moved[k])for(int f=0;f<int(trial[k].Triangles.size());++f){const auto &t=trial[k].Triangles[f];
          if(!changedVertices[k][t[0]]&&!changedVertices[k][t[1]]&&!changedVertices[k][t[2]])continue;changedFaces[k][f]=1;
          const V p[3]={toEigenVec(trial[k].Positions[t[0]]),toEigenVec(trial[k].Positions[t[1]]),toEigenVec(trial[k].Positions[t[2]])};
          const V original[3]={toEigenVec(work[k].original.Positions[t[0]]),toEigenVec(work[k].original.Positions[t[1]]),toEigenVec(work[k].original.Positions[t[2]])};
          const V n=(p[1]-p[0]).cross(p[2]-p[0]),old=(original[1]-original[0]).cross(original[2]-original[0]);
          if(n.squaredNorm()<1e-12*old.squaredNorm()||n.dot(old)<=0){valid[k]=0;errors[k]="triangle collapse/reversal";blocked[k].insert(t.begin(),t.end());break;}}

      }catch(...){geometry.capture(k);}}
      geometry.rethrow();if(std::find(valid.begin(),valid.end(),char(0))!=valid.end()){
        bool reduced=false;for(size_t k=1;k+1<layers.size();++k)if(!valid[k]){blockedUnion[k].insert(blocked[k].begin(),blocked[k].end());
          rejected=errors[k]+" on layer "+std::to_string(k);reduced=reduceLayer(int(k))||reduced;}
        if(reduced)continue;break;}

      zTsDetail::ParallelErrors refits(layers.size());
      #pragma omp parallel for schedule(dynamic) if(parallel)
      for(int k=0;k<int(layers.size());++k){try{const bool needed=moved[k]||(k>0&&moved[k-1])||(k+1<int(layers.size())&&moved[k+1]);
        if(needed||dirtyTrees[k]){if(trialTrees[k])trialTrees[k]->refit(trial[k].Positions);else trialTrees[k].reset(new Tree(trial[k].Positions,trial[k].Triangles));dirtyTrees[k]=moved[k];}
      }catch(...){refits.capture(k);}}
      refits.rethrow();std::vector<int> order;for(int k=1;k<int(layers.size())-1;++k)if(moved[k])order.push_back(k);
      std::sort(order.begin(),order.end(),[&](int a,int b){return geometryPriority[a]!=geometryPriority[b]?geometryPriority[a]>geometryPriority[b]:a<b;});
      // A source/self-contact test depends only on this layer's candidate.
      // Reuse its successful validation while another layer backtracks.
      #pragma omp parallel for schedule(dynamic) if(parallel)
      for(int j=0;j<int(order.size());++j){const int k=order[j];try{
        if(validatedSteps[k]==lineSteps[k])continue;
        valid[k]=validTrial(trial[k],work[k],side,sample,tol,errors[k],blocked[k],&changedVertices[k],&changedFaces[k],trialTrees[k].get());
        if(valid[k])validatedSteps[k]=lineSteps[k];
      }catch(...){geometry.capture(k);}}
      geometry.rethrow();bool rejectedGeometry=false,reducedGeometry=false;
      for(int k:order)if(!valid[k]){++geometryPriority[k];blockedUnion[k].insert(blocked[k].begin(),blocked[k].end());
        rejected=errors[k]+" on layer "+std::to_string(k);rejectedGeometry=true;reducedGeometry=reduceLayer(k)||reducedGeometry;}
      if(rejectedGeometry){if(reducedGeometry)continue;break;}
      std::vector<char> crossed(layers.size()-1,0);std::vector<std::pair<int,int>> contactFaces(layers.size()-1,{-1,-1});zTsDetail::ParallelErrors separation(layers.size()-1);
      #pragma omp parallel for schedule(dynamic) if(parallel)
      for(int k=0;k<int(layers.size())-1;++k){try{if(moved[k]||moved[k+1])crossed[k]=!separatedTrees(*trialTrees[k],*trialTrees[k+1],tol,&contactFaces[k]);}catch(...){separation.capture(k);}}
      separation.rethrow();if(std::find(crossed.begin(),crossed.end(),char(1))!=crossed.end()){
        for(int k=0;k<int(crossed.size());++k)if(crossed[k]){const auto pair=contactFaces[k];
          if(pair.first<0||pair.second<0||!contactKeys.emplace(k,pair.first,pair.second).second)continue;
          Face a,b;for(int j=0;j<3;++j){a.p[j]=toEigenVec(layers[k].Positions[layers[k].Triangles[pair.first][j]]);
            b.p[j]=toEigenVec(layers[k+1].Positions[layers[k+1].Triangles[pair.second][j]]);}
          V normal=V::Zero();double gap=0;
          const auto axis=[&](V n){if(n.norm()<1e-12)return;n.normalize();double alo=inf,ahi=-inf,blo=inf,bhi=-inf;
            for(int j=0;j<3;++j){const double pa=n.dot(a.p[j]),pb=n.dot(b.p[j]);alo=(std::min)(alo,pa);ahi=(std::max)(ahi,pa);blo=(std::min)(blo,pb);bhi=(std::max)(bhi,pb);}
            if(blo-ahi>gap){gap=blo-ahi;normal=n;}if(alo-bhi>gap){gap=alo-bhi;normal=-n;}};
          axis(a.normal());axis(b.normal());for(int i=0;i<3;++i)for(int j=0;j<3;++j)axis((a.p[(i+1)%3]-a.p[i]).cross(b.p[(j+1)%3]-b.p[j]));
          if(gap>tol){ContactGuard guard;guard.Layer=k;guard.FaceA=pair.first;guard.FaceB=pair.second;
            guard.Normal=normal;guard.Minimum=(std::max)(4*tol,(std::min)(.1*s.Minimum,.25*gap));contactGuards.push_back(guard);}}
        bool reduced=false;std::vector<char> collisionLayers(layers.size(),0);
        for(int k=0;k<int(crossed.size());++k)if(crossed[k])collisionLayers[k]=collisionLayers[k+1]=1;
        for(int k=1;k<int(collisionLayers.size())-1;++k)if(collisionLayers[k])reduced=reduceLayer(k)||reduced;
        rejected="Neighbouring surfaces intersect.";if(reduced)continue;break;}
      std::vector<char> affected(layers.size(),0);for(int k=0;k<int(layers.size());++k)if(moved[k])
        for(int j=(std::max)(0,k-1);j<=(std::min)(int(layers.size())-1,k+1);++j)affected[j]=1;
      HeightLinks nextLinks=links;Report measured;std::string queryError;
      if(!refreshHeights(trial,tol,measured,queryError,parallel,&nextLinks,&affected,&failedQuery,&trialTrees,true)){
        rejected=queryError;const bool a=reduceLayer(failedQuery.LayerId),b=reduceLayer(failedQuery.Adjacent);
        if(a||b)continue;break;}
      auto nextSamples=samples;
      if(!refreshExactSamples(trial,nextSamples,trialTrees,maximumRay,tol,queryError,parallel,&affected,&failedQuery)){
        rejected=queryError;const bool a=reduceLayer(failedQuery.LayerId),b=reduceLayer(failedQuery.Adjacent);
        if(a||b)continue;break;}
      HeightReport next=report;const double nextScore=heightScore(trial,nextSamples,s,next);
      // Never buy aggregate least-squares progress by worsening the largest
      // physical TCP violation or creating a new near-contact gap.
      double previousGap=inf;for(const auto &layer:links)for(const auto &vertex:layer)for(const auto &hit:vertex)previousGap=(std::min)(previousGap,hit.Height);
      const double gapFloor=(std::min)(previousGap,.1*s.Minimum);bool nearContact=false,reducedContact=false;std::vector<char> contactLayers(layers.size(),0);
      for(int k=0;k<int(nextLinks.size());++k)for(const auto &vertex:nextLinks[k])for(const auto &hit:vertex)
        if(hit.Height+tol<gapFloor){nearContact=true;contactLayers[k]=contactLayers[hit.LayerId]=1;}
      const double targetFloor=(std::min)(report.Minimum,.1*s.Minimum);
      for(int k=1;k<int(nextSamples.size());++k)for(const auto &sample:nextSamples[k])
        if(sampleHeight(trial[k],sample)+tol<targetFloor){nearContact=true;contactLayers[k]=contactLayers[k-1]=1;}
      for(int k=1;k<int(contactLayers.size())-1;++k)if(contactLayers[k])reducedContact=reduceLayer(k)||reducedContact;
      if(nearContact){rejected="Surface update creates a new near-contact gap.";if(reducedContact)continue;break;}
      const bool preservesWorst=next.WorstViolation<=report.WorstViolation+1e-7*s.Maximum;
      if(preservesWorst&&(next.Converged||nextScore<score*(1-1e-7))){
        const double improvement=(score-nextScore)/(std::max)(score,1e-30);stalled=improvement<1e-4?stalled+1:0;
        layers.swap(trial);links.swap(nextLinks);samples.swap(nextSamples);score=nextScore;report=next;++report.AcceptedIterations;accepted=true;
        acceptedLine=1;for(size_t k=1;k+1<lineSteps.size();++k)if(lineSteps[k]>0)acceptedLine=(std::min)(acceptedLine,lineSteps[k]);
        directedGaps();
        // Transport scalar constraint forces across frozen-coefficient batches.
        // Resetting them every outer iteration turns a hard bound into a weak
        // finite-iteration penalty and stalls before the requested domain.
        for(const auto &row:rows){if(row.Layer<0)continue;double force=0,norm=0;for(const auto &t:row.Terms){force+=t.Weight*t.A.dot(t.Dual);norm+=t.A.squaredNorm();}
          if(norm>0){const double nextMultiplier=force/norm;auto &old=multipliers[row.Layer][row.Sample];
            double trust=1;for(const auto &term:row.Terms)trust=(std::min)(trust,lineSteps[term.Layer]*scales[term.Layer]);
            old+=trust*(nextMultiplier-old);}}
        break;}
      rejected="Frozen-normal surface update did not reduce measured height error.";
      bool reduced=false;for(int k=1;k<int(layers.size())-1;++k)reduced=reduceLayer(k)||reduced;
      if(!reduced)break;
    }
    if(!accepted){bool changed=false;
      for(size_t k=1;k+1<layers.size();++k)blocked[k].insert(blockedUnion[k].begin(),blockedUnion[k].end());
      const auto freezeRing=[&](int layer,int id){if(layer<=0||layer+1>=int(layers.size())||id<0)return;
        blocked[layer].insert(id);for(const auto &tri:layers[layer].Triangles)if(std::find(tri.begin(),tri.end(),id)!=tri.end())blocked[layer].insert(tri.begin(),tri.end());};
      if(failedQuery.LayerId>=0&&failedQuery.Vertex>=0){freezeRing(failedQuery.LayerId,failedQuery.Vertex);
        for(const auto &hit:links[failedQuery.LayerId][failedQuery.Vertex])if(hit.LayerId==failedQuery.Adjacent)
          for(int id:hit.Ids)freezeRing(hit.LayerId,id);}
      for(size_t k=1;k+1<layers.size();++k){bool rebuild=false;for(int id:blocked[k])if(!work[k].fixed[id]){work[k].fixed[id]=true;rebuild=changed=true;}
        if(rebuild){if(!surfaceSystem(layers[k],work[k],systems[k],error))return false;++report.Factorizations;}}
      if(changed&&!accepted){if(!refreshHeights(layers,tol,ignored,error,parallel,&links,nullptr,nullptr,&trees))return false;continue;}
      if(!accepted&&contactGuards.size()>initialContacts)continue;
      if(!accepted){report.LastRejectedReason=rejected.empty()?"Surface optimization stagnated.":rejected;break;}}
    {std::ostringstream message;message<<"[surfaceIteration] outer="<<report.OuterIterations<<" constraints="<<rows.size()<<" domain=["<<report.Minimum<<","<<report.Maximum
      <<"] violations="<<report.ViolatingVertices<<" uncontrollable="<<report.UncontrollableSamples
      <<" step="<<acceptedLine<<" minimumLayerTrust="<<minimumTrust<<" ms="<<std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-outerStart).count();zTsUtils::printDebug(message.str());}
    if(report.Converged)break;if(stalled>=5){report.LastRejectedReason="Surface optimization stagnated.";break;}
  }
  for(size_t k=1;k+1<layers.size();++k)for(size_t i=0;i<layers[k].Positions.size();++i){const V p=toEigenVec(layers[k].Positions[i]);
    report.MaxDisplacement=(std::max)(report.MaxDisplacement,(p-toEigenVec(work[k].original.Positions[i])).norm());
    if(!work[k].allowed[i].empty()){int t;V q;report.MaxBoundaryError=(std::max)(report.MaxBoundaryError,side.tree.nearest(p,t,q));}}
  report.FinalBendingEnergy=bendingEnergy(layers,work);
  for(size_t k=1;k+1<work.size();++k)report.FixedVertices+=int(std::count(work[k].fixed.begin(),work[k].fixed.end(),true));
  if(!report.Converged&&report.LastRejectedReason.empty())report.LastRejectedReason="Surface iteration budget exhausted; domain not yet achieved.";
  if(report.Converged)report.LastRejectedReason.clear();finish();return true;
}

} } // namespace zSpace::zTsRelax
