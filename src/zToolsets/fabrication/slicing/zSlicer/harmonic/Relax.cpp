#include <zspace/zToolsets/fabrication/slicing/zSlicer/HarmonicSlicer.h>
#include <algorithm>
#include <numeric>
#include <sstream>

namespace zSpace {
namespace {
// Immutable point locator for the existing tetrahedral coordinate field. Its
// purpose is containment and transported UV; it does not reimpose the old iso.
class TetSampler {
  using V=Eigen::Vector3d;
  struct Bounds { V lo=V::Constant(INFINITY),hi=V::Constant(-INFINITY); };
  struct Node { Bounds box;int first=0,last=0,left=-1,right=-1; };
  const HarmonicField &tet_;
  std::vector<Bounds> boxes_;
  std::vector<int> order_;
  std::vector<Node> nodes_;
  double tolerance_=1e-9;
  int build(int first,int last) {
    const int id=int(nodes_.size());nodes_.push_back(Node());nodes_[id].first=first;nodes_[id].last=last;
    for(int i=first;i<last;++i){nodes_[id].box.lo=nodes_[id].box.lo.cwiseMin(boxes_[order_[i]].lo);nodes_[id].box.hi=nodes_[id].box.hi.cwiseMax(boxes_[order_[i]].hi);}
    if(last-first>8){Eigen::Index axis;(nodes_[id].box.hi-nodes_[id].box.lo).maxCoeff(&axis);const int mid=(first+last)/2;
      std::nth_element(order_.begin()+first,order_.begin()+mid,order_.begin()+last,[&](int a,int b){return boxes_[a].lo[axis]+boxes_[a].hi[axis]<boxes_[b].lo[axis]+boxes_[b].hi[axis];});
      const int left=build(first,mid),right=build(mid,last);nodes_[id].left=left;nodes_[id].right=right;}
    return id;
  }
public:
  explicit TetSampler(const HarmonicField &tet):tet_(tet),boxes_(tet.Volumes.size()) {
    for(size_t t=0;t<tet.Volumes.size();++t)if(tet.Volumes[t].Valid){for(int id:tet.Volumes[t].VertexIds){const V p=zTsUtils::toEigenVec(tet.Vertices[id].Position);boxes_[t].lo=boxes_[t].lo.cwiseMin(p);boxes_[t].hi=boxes_[t].hi.cwiseMax(p);}order_.push_back(int(t));}
    if(!order_.empty()){build(0,int(order_.size()));tolerance_=(std::max)(1e-9,2e-7*(nodes_[0].box.hi-nodes_[0].box.lo).norm());}
  }
  bool sample(const zPoint &position,zPoint *uv) const {
    if(nodes_.empty())return false;const V p=zTsUtils::toEigenVec(position);std::array<int,64> pending;int count=0;pending[count++]=0;
    while(count){const Node &node=nodes_[pending[--count]];
      const double tol=tolerance_;
      if((p.array()<node.box.lo.array()-tol).any()||(p.array()>node.box.hi.array()+tol).any())continue;
      if(node.left>=0){pending[count++]=node.left;pending[count++]=node.right;continue;}
      for(int j=node.first;j<node.last;++j){const int t=order_[j];const auto &ids=tet_.Volumes[t].VertexIds;
        const V delta=p-zTsUtils::toEigenVec(tet_.Vertices[ids[0]].Position);
        Eigen::Vector4d b=tet_.VolumeBasisGradients[t].transpose()*delta;b[0]+=1;
        bool inside=true;for(int i=0;i<4;++i){const double btol=2e-6+tol*tet_.VolumeBasisGradients[t].col(i).norm();
          if(b[i]<-btol||b[i]>1+btol){inside=false;break;}}
        if(!inside)continue;
        if(!uv)return true;
        b=b.cwiseMax(0);b/=b.sum();V value=V::Zero();bool valid=true;
        for(int i=0;i<4;++i)if(b[i]>1e-10){const V v=zTsUtils::toEigenVec(tet_.UVW[ids[i]]);if(!v.allFinite()){valid=false;break;}value+=b[i]*v;}
        if(valid&&value.allFinite()){*uv=zPoint(value.x(),value.y(),value.z());return true;}
      }
    }return false;
  }
};

// Cap triangles retain the solid's outward winding for all solver operations.
// Fabrication starts at the bottom and deposits into the solid, toward higher
// harmonic layers. Supply that inward direction only to path projection.
// Rebuild from geometry, never negate a previous override (repeated runs are
// idempotent). relaxationLayers() intentionally does not import this table.
bool prepareBottomPrintNormals(zObjMeshArray &meshes,
    std::vector<zVectorArray> &normals,std::string &error) {
  zTsUtils::ScopedSlicingTiming methodTiming("HarmonicSlicer.preparePrintCapNormals");
  if(meshes.empty()){error="No bottom cap for fabrication normals.";return false;}
  zPointArray points;zIntArray counts,connects;zFnMesh cap(meshes.front());
  cap.getVertexPositions(points);cap.getPolygonData(connects,counts);
  std::vector<Eigen::Vector3d> inward(points.size(),Eigen::Vector3d::Zero());size_t offset=0;
  for(int count:counts){for(int i=1;i+1<count;++i){const int ids[]={connects[offset],connects[offset+i],connects[offset+i+1]};
      const auto a=zTsUtils::toEigenVec(points[ids[0]]),b=zTsUtils::toEigenVec(points[ids[1]]),c=zTsUtils::toEigenVec(points[ids[2]]);
      const Eigen::Vector3d normal=(c-a).cross(b-a);
      if(!normal.allFinite()){error="Nonfinite bottom-cap fabrication normal.";return false;}
      for(int id:ids)inward[id]+=normal;
    }offset+=count;
  }
  zVectorArray result;result.reserve(points.size());
  for(auto &normal:inward){const double length=normal.norm();
    if(!std::isfinite(length)||length<=0){error="Unsupported bottom-cap fabrication normal.";return false;}
    normal/=length;result.emplace_back(normal.x(),normal.y(),normal.z());
  }
  if(result.empty()){error="Empty bottom-cap fabrication normal table.";return false;}
  normals.resize(meshes.size());normals.front()=std::move(result);return true;
}
}

std::vector<zTsRelax::Layer> HarmonicSlicer::relaxationLayers() const {
  std::vector<zTsRelax::Layer> layers(o_IsoMeshes.size());
  for(size_t k=0;k<layers.size();++k){auto &layer=layers[k];zFnMesh fn(const_cast<zObjMesh&>(o_IsoMeshes[k]));fn.getVertexPositions(layer.Positions);
    zIntArray connects,counts;fn.getPolygonData(connects,counts);size_t offset=0;
    for(int count:counts){if(count!=3)throw std::runtime_error("Iso relaxation requires triangulated layers.");
      layer.Triangles.push_back({connects[offset],connects[offset+1],connects[offset+2]});offset+=count;}
    layer.UV=o_IsoMeshUV[k];layer.Heights=o_VertexLayerHeights[k];layer.BoundaryWidths=_HarmonicIsoEdgeWidthIds[k];
  }
  return layers;
}

bool HarmonicSlicer::relaxContours(std::string &error, bool useParallel) {
  zTsUtils::ScopedSlicingTiming methodTiming("HarmonicSlicer.relaxContours");
  o_VertexLayerNormals.assign(o_IsoMeshes.size(),zVectorArray());
  zTsRelax::Settings settings;
  settings.Iterations=i_Settings._HIsoRelaxIterations;
  settings.Strength=i_Settings._HIsoRelaxStrength;
  settings.MaxDeviation=i_Settings._HIsoRelaxMaxDeviation;
  // Sliver preconditioning belongs to the optional surface-height solver.
  // Iso-value redistribution alone keeps the ordinary relaxation mesh.
  if(i_Settings.HPrintHeight.Enabled&&i_Settings.HPrintHeight.OptimizeSurfaces)
    settings.MinEdgeLength=.05*i_Settings.HPrintHeight.Minimum;
  settings.Parallel=useParallel;
  zTsRelax::Reference reference;reference.Positions=_Tet.Surface.Positions;reference.Triangles=_Tet.Surface.Triangles;
  for(int face:_Tet.SurfaceFaceIds)reference.WidthIds.push_back(_Tet.FaceTags[face]==HarmonicField::TetFaceTag::Side?_Tet.FaceGroupIds[face]:-1);
  auto layers=relaxationLayers();
  TetSampler sampler(_Tet);
  if(!zTsRelax::relax(reference,layers,settings,[&](const zPoint &p,zPoint *uv){return sampler.sample(p,uv);},_IsoRelaxationReport,error))return false;
  if(settings.Iterations==0)return prepareBottomPrintNormals(o_IsoMeshes,o_VertexLayerNormals,error);
  for(size_t k=0;k<layers.size();++k){if(k>0&&k+1<layers.size()){
      if(settings.MinEdgeLength>0){zIntArray counts(layers[k].Triangles.size(),3),connects;for(const auto &t:layers[k].Triangles)connects.insert(connects.end(),t.begin(),t.end());zFnMesh(o_IsoMeshes[k]).create(layers[k].Positions,counts,connects);}
      else zFnMesh(o_IsoMeshes[k]).setVertexPositions(layers[k].Positions);}
    _HarmonicIsoEdgeWidthIds[k].swap(layers[k].BoundaryWidths);
    o_IsoMeshUV[k].swap(layers[k].UV);o_VertexLayerHeights[k].swap(layers[k].Heights);o_VertexLayerNormals[k].swap(layers[k].Normals);}
  std::ostringstream report;report<<"[relaxContours] iterations="<<_IsoRelaxationReport.AcceptedIterations
      <<" boundaryVerticesMoved="<<_IsoRelaxationReport.MovedBoundaryVertices
      <<" maxDisplacement="<<_IsoRelaxationReport.MaxDisplacement
      <<" maxBoundaryError="<<_IsoRelaxationReport.MaxBoundaryError
      <<" minimumInitialGap="<<_IsoRelaxationReport.MinimumInitialGap
      <<" recoveredBoundaryTags="<<_IsoRelaxationReport.RecoveredBoundaryTags
      <<" uvConstrainedVertices="<<_IsoRelaxationReport.UVConstrainedVertices
      <<" sliverConstrainedVertices="<<_IsoRelaxationReport.SliverConstrainedVertices
      <<" geometryConstrainedVertices="<<_IsoRelaxationReport.GeometryConstrainedVertices
      <<" lastRejectedStep="<<_IsoRelaxationReport.LastRejectedReason
      <<" boundaryHeightPlaneEstimates="<<_IsoRelaxationReport.PlaneHeightEstimates;
  zTsUtils::printDebug(report.str());return prepareBottomPrintNormals(o_IsoMeshes,o_VertexLayerNormals,error);
}

bool HarmonicSlicer::optimizeSurfaceHeightDomain(const std::vector<zPointArray> *targets,std::string &error,bool parallel){
  zTsUtils::ScopedSlicingTiming methodTiming("HarmonicSlicer.optimizeSurfaceHeightDomain");
  zTsRelax::Reference reference;reference.Positions=_Tet.Surface.Positions;reference.Triangles=_Tet.Surface.Triangles;
  for(int face:_Tet.SurfaceFaceIds)reference.WidthIds.push_back(_Tet.FaceTags[face]==HarmonicField::TetFaceTag::Side?_Tet.FaceGroupIds[face]:-1);
  auto layers=relaxationLayers();TetSampler sampler(_Tet);
  auto settings=i_Settings.HPrintHeight;
  if(!targets)settings.Iterations=(std::min)(settings.Iterations,8);
  else if(settings.Tolerance>0){
    // Keep a small guide margin for float storage and a rebuilt chart. Final
    // acceptance still uses the original user domain, on the complete curves.
    const double inset=(std::min)(.05*(settings.Maximum-settings.Minimum),.5*settings.Tolerance);
    settings.Minimum+=inset;settings.Maximum-=inset;
  }
  if(!zTsRelax::optimizeHeightDomain(reference,layers,settings,[&](const zPoint &p,zPoint *uv){return sampler.sample(p,uv);},parallel,_PrintHeightOptimizationReport,error,targets,i_Settings._HIsoRelaxStrength))return false;
  for(size_t k=0;k<layers.size();++k){if(k>0&&k+1<layers.size())zFnMesh(o_IsoMeshes[k]).setVertexPositions(layers[k].Positions);
    o_IsoMeshUV[k].swap(layers[k].UV);o_VertexLayerHeights[k].swap(layers[k].Heights);o_VertexLayerNormals[k].swap(layers[k].Normals);}
  const auto &h=_PrintHeightOptimizationReport;std::ostringstream message;message<<"[surfaceHeightOptimizer] iterations="<<h.AcceptedIterations<<" factors="<<h.Factorizations<<" domain=["<<h.Minimum<<","<<h.Maximum
    <<"] violations="<<h.ViolatingVertices<<" worst="<<h.WorstViolation<<" timeMs="<<h.Milliseconds<<" constrained="<<h.LastRejectedReason;
  zTsUtils::printDebug(message.str());return prepareBottomPrintNormals(o_IsoMeshes,o_VertexLayerNormals,error);
}
} // namespace zSpace
