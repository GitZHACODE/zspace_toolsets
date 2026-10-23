#include <zspace/zToolsets/fabrication/slicing/utils/IsoRelaxation.h>
#include <zspace/zToolsets/fabrication/slicing/utils/Parallel.h>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <queue>
#include <set>
#include <sstream>

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
    const V q=p-n*((p-f.p[0]).dot(n)/n.squaredNorm());
    if (barycentric(q,f).minCoeff()>=0) return q;
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
bool segmentInterval(const V &a,const V &b,const Face &f,double tol,double &lo,double &hi) {
  const V n=f.normal(); const double nn=n.norm();
  if (!(nn>0) || std::abs((a-f.p[0]).dot(n))>tol*nn ||
      std::abs((b-f.p[0]).dot(n))>tol*nn) return false;
  const V ba=barycentric(a,f),db=barycentric(b,f)-ba;
  lo=0; hi=1;
  for (int j=0;j<3;++j) {
    const double btol=tol*(f.p[(j+1)%3]-f.p[(j+2)%3]).norm()/nn;
    if (std::abs(db[j])<1e-15) { if (ba[j]<-btol) return false; }
    else if (db[j]>0) lo=(std::max)(lo,(-btol-ba[j])/db[j]);
    else hi=(std::min)(hi,(-btol-ba[j])/db[j]);
  }
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
  const zPointArray &p_; const std::vector<Tri> &tri_;
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
  Tree(const zPointArray &p,const std::vector<Tri> &t):p_(p),tri_(t) {
    faces_.resize(t.size());boxes_.resize(t.size());
    for(size_t i=0;i<t.size();++i){for(int j=0;j<3;++j)faces_[i].p[j]=toEigenVec(p[t[i][j]]);boxes_[i]=bounds(faces_[i]);}
    order_.resize(t.size());std::iota(order_.begin(),order_.end(),0);
    if (!t.empty()) build(0,int(t.size()));
  }
  const Face &face(int i) const { return faces_[i]; }
  template<class F> bool visit(const Box &box,double tol,F callback) const {
    if (nodes_.empty()) return true;
    std::array<int,64> pending;int count=0;pending[count++]=0;
    while(count) {const Node &n=nodes_[pending[--count]];
      if(!n.box.overlaps(box,tol))continue;
      if(n.left>=0){pending[count++]=n.left;pending[count++]=n.right;}
      else for(int i=n.begin;i<n.end;++i)if(boxes_[order_[i]].overlaps(box,tol)&&!callback(order_[i]))return false;
    }return true;
  }
  double nearest(const V &p,int &triangle,V &q) const {
    double best=inf;triangle=-1;std::array<int,64> pending;int count=0;
    if(!nodes_.empty())pending[count++]=0;
    while(count){const Node &n=nodes_[pending[--count]];
      if(n.box.distanceSquared(p)>best)continue;
      if(n.left>=0){pending[count++]=n.left;pending[count++]=n.right;}
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
  double ray(const V &p,const V &direction,double maximum,double tol) const {
    double best=maximum;bool found=false;Box box;box.add(p);box.add(p+maximum*direction);
    visit(box,tol,[&](int id){const Face f=face(id);const V n=f.normal();const double den=n.dot(direction);
      if(std::abs(den)<=1e-12*n.norm())return true;
      const double d=n.dot(f.p[0]-p)/den;
      if(d>tol&&d<best&&barycentric(p+d*direction,f).minCoeff()>=-1e-8){best=d;found=true;}return true;});
    return found?best:inf;
  }
};

struct SidePatch { std::vector<int> triangles;int width=-1; };
struct BoundaryEdge { int a,b,patch; };
struct SideConstraints {
  const Reference &ref;Tree tree;double tol;
  std::vector<SidePatch> patches;
  std::map<Key,std::vector<std::pair<int,int>>> curves;
  SideConstraints(const Reference &r,double tolerance):ref(r),tree(r.Positions,r.Triangles),tol(tolerance) {}
  bool build(std::string &error) {
    const int n=int(ref.Triangles.size());std::vector<int> parent(n);std::iota(parent.begin(),parent.end(),0);
    const auto root=[&](int v){while(parent[v]!=v)v=parent[v];return v;};
    std::map<Key,std::vector<int>> edgeFaces;
    for(int t=0;t<n;++t)for(int j=0;j<3;++j)edgeFaces[makeEdgeKey(ref.Triangles[t][j],ref.Triangles[t][(j+1)%3])].push_back(t);
    for(const auto &entry:edgeFaces){const auto &ids=entry.second;if(ids.size()!=2)continue;
      const int a=ids[0],b=ids[1];if(ref.WidthIds[a]<0||ref.WidthIds[a]!=ref.WidthIds[b])continue;
      const Face fa=tree.face(a),fb=tree.face(b);const V na=fa.normal().normalized(),nb=fb.normal().normalized();
      if(std::abs(na.dot(nb))<1-1e-10)continue;
      bool planar=true;for(const V &p:fb.p)planar &= std::abs((p-fa.p[0]).dot(na))<=tol;
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
    for(int t:patches[patch].triangles){double lo,hi;if(segmentInterval(a,b,tree.face(t),tol,lo,hi))intervals.push_back({lo,hi});}
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
  V planeOrigin=V::Zero(),planeNormal=V::Zero();
  bool planar=false;
  double limit=0;
};

bool prepare(Layer &m,const SideConstraints &side,double gap,const Settings &settings,Work &w,Report &report,std::string &error) {
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
  const Sparse matrix=lumped+lambda*(lap.transpose()*inv*lap);
  Eigen::SimplicialLDLT<Sparse> solver;solver.compute(matrix);
  if(solver.info()!=Eigen::Success){error="Relaxation matrix factorisation failed.";return false;}
  Eigen::MatrixXd rhs(n,3);for(int i=0;i<n;++i)rhs.row(i)=mass[i]*(toEigenVec(m.Positions[i])-w.planeOrigin).transpose();
  w.target=solver.solve(rhs);
  w.target.rowwise()+=w.planeOrigin.transpose();
  if(solver.info()!=Eigen::Success||!w.target.allFinite()){error="Relaxation solve returned invalid positions.";return false;}return true;
}

bool validTrial(const Layer &m,const Work &w,const SideConstraints &side,const SampleVolume &sample,double tol,std::string &reason,std::set<int> &blocked) {
  blocked.clear();
  // Projection/containment need a float-coordinate tolerance. Reusing that
  // expanded slab for intersection predicates invents overlaps between nearly
  // coplanar neighbouring triangles. Intersect the stored geometry in double.
  const double intersectionTol=tol*1e-3;
  for(size_t i=0;i<m.Positions.size();++i)if(!sample(m.Positions[i],nullptr)){reason="volume containment";blocked.insert(int(i));return false;}
  for(const auto &edge:w.edges)if(!side.covers(edge.patch,toEigenVec(m.Positions[edge.a]),toEigenVec(m.Positions[edge.b]))){reason="boundary support";blocked.insert(edge.a);blocked.insert(edge.b);return false;}
  Tree tree(m.Positions,m.Triangles);
  for(size_t i=0;i<m.Triangles.size();++i){const Tri &t=m.Triangles[i];const Face f=tree.face(int(i));
    const Face old={ {toEigenVec(w.original.Positions[t[0]]),toEigenVec(w.original.Positions[t[1]]),toEigenVec(w.original.Positions[t[2]])} };
    if(f.normal().squaredNorm()<1e-12*old.normal().squaredNorm()||f.normal().dot(old.normal())<=0){reason="triangle collapse/reversal";blocked.insert(t.begin(),t.end());return false;}
    // Check the complete triangle against the immutable volume boundary. Contacts
    // are allowed only at its genuine side-boundary vertices/edges.
    if(!side.tree.visit(bounds(f),tol,[&](int r){for(const V &hit:contacts(f,side.tree.face(r),intersectionTol)){
        bool boundary=false;
        for(int j=0;j<3;++j){const int a=t[j],b=t[(j+1)%3];
          if(!w.allowed[a].empty()&&(hit-f.p[j]).norm()<=4*tol)boundary=true;
          if(m.BoundaryWidths.count(makeEdgeKey(a,b))&&(closestSegment(hit,f.p[j],f.p[(j+1)%3])-hit).norm()<=4*tol)boundary=true;}
        if(!boundary)return false;}return true;})){reason="input-surface contact";blocked.insert(t.begin(),t.end());return false;}
    if(!tree.visit(bounds(f),tol,[&](int j){if(j<=int(i))return true;const Tri &b=m.Triangles[j];
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

bool refreshHeights(std::vector<Layer> &layers,double tol,Report &report,std::string &error,bool parallel) {
  std::vector<std::unique_ptr<Tree>> trees;for(const Layer &m:layers)trees.emplace_back(new Tree(m.Positions,m.Triangles));
  std::vector<zVectorArray> n;for(const Layer &m:layers)n.push_back(normals(m));
  std::vector<int> planeEstimates(layers.size(),0);std::vector<std::string> errors(layers.size());
  const auto refresh=[&](size_t k){std::string &error=errors[k];Layer &m=layers[k];m.Heights.resize(m.Positions.size());
    if(k>0&&k+1<layers.size())m.Normals=n[k];
    for(size_t i=0;i<m.Positions.size();++i){const V p=toEigenVec(m.Positions[i]),normal=toEigenVec(n[k][i]);
      if(normal.norm()<0.5){error="Cannot recover a finite layer normal.";return false;}
      double total=0;int count=0;
      for(int offset:{-1,1}){const int adjacent=int(k)+offset;if(adjacent<0||adjacent>=int(layers.size()))continue;
        int tri=-1;V q;const double nearest=trees[adjacent]->nearest(p,tri,q);
        if(!(nearest>tol)||tri<0){error="Relaxed layers have zero or invalid separation.";return false;}
        V direction=offset*normal;
        // The bottom cap's input winding is outward (opposite the slicing
        // direction); an imported cap may also use the opposite convention.
        if((k==0||k+1==layers.size())&&direction.dot(q-p)<0)direction=-direction;
        double h=trees[adjacent]->ray(p,direction,4*nearest+tol,tol);
        if(!std::isfinite(h)){
          // A side-boundary normal ray can leave the solid before hitting the
          // adjacent trimmed layer. Extrapolate its nearest local face plane;
          // this is computed from relaxed geometry, not the old field height.
          const Face face=trees[adjacent]->face(tri);
          const V weights=barycentric(q,face);
          V faceNormal=V::Zero();for(int j=0;j<3;++j)faceNormal+=weights[j]*toEigenVec(n[adjacent][layers[adjacent].Triangles[tri][j]]);
          faceNormal.normalize();
          const double den=faceNormal.dot(direction);
          if(std::abs(den)<0.25){error="Adjacent layer is nearly tangent to the print-height ray.";return false;}
          h=faceNormal.dot(q-p)/den;++planeEstimates[k];
        }
        if(!std::isfinite(h)||h<=tol){error="Cannot estimate positive spacing on layer "+std::to_string(k)+", vertex "+std::to_string(i)+", neighbour "+std::to_string(adjacent)+".";return false;}
        total+=h;++count;
      }
      m.Heights[i]=float(total/count);
    }
    return true;
  };
  zTsDetail::ParallelErrors failures(layers.size());
  #pragma omp parallel for schedule(dynamic) if(parallel)
  for(int k=0;k<int(layers.size());++k){try{refresh(size_t(k));}catch(...){failures.capture(k);}}
  failures.rethrow();
  for(size_t k=0;k<layers.size();++k){report.PlaneHeightEstimates+=planeEstimates[k];if(!errors[k].empty()){error=errors[k];return false;}}
  return true;
}
} // namespace

bool relax(const Reference &reference,std::vector<Layer> &layers,const Settings &settings,
           const SampleVolume &sample,Report &report,std::string &error) {
  report=Report();error.clear();
  using Clock=std::chrono::steady_clock;const auto start=Clock::now();auto last=start;
  const auto mark=[&](int stage){const auto now=Clock::now();report.StageMilliseconds[stage]=std::chrono::duration<double,std::milli>(now-last).count();last=now;};
  if(settings.Iterations<0||settings.Iterations>100||!std::isfinite(settings.Strength)||settings.Strength<=0||settings.Strength>1||
     !std::isfinite(settings.MaxDeviation)||settings.MaxDeviation<0){error="Invalid iso relaxation settings (iterations 0..100, strength (0,1], deviation >=0).";return false;}
  if(settings.Iterations==0)return true;
  if(layers.size()<2||!sample||reference.Triangles.size()!=reference.WidthIds.size()){error="Missing layer, reference or volume-query data.";return false;}
  Layer ref;ref.Positions=reference.Positions;ref.Triangles=reference.Triangles;if(!meshValid(ref,error))return false;
  for(const Layer &m:layers)if(!meshValid(m,error))return false;
  Box box;for(const auto &p:reference.Positions)box.add(toEigenVec(p));
  const double tol=(std::max)(1e-9,2e-7*(box.hi-box.lo).norm());
  SideConstraints side(reference,tol);if(!side.build(error))return false;
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

} } // namespace zSpace::zTsRelax
