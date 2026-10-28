// Bindings for the plain data structures shared across the slicer toolset:
// zTsSlicerMode (zTsSlicing.h), SlicerSettings (Utils.h), PrintSynSetting
// (zPrintSyn.h), ParameterGraph
// and PrintPath (Utils.h), plus a minimal opaque binding for
// zObjMeshScalarField (owned by zTsSlicing::getScalarField() but otherwise
// outside the requested binding scope - it is a zSpace Core type, not one
// of zObjMesh/zFnMesh/zObjGraph/zFnGraph).

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <zspace/zToolsets/fabrication/slicing/zTsSlicing.h>
#include <zspace/zInterface/objects/zObjMeshField.h>

namespace py = pybind11;
using namespace zSpace;

namespace zspace_python {

void bind_types(py::module_ &m) {
  py::enum_<zTsSlicerMode>(m, "SlicerMode")
      .value("Topology", zTsSlicerMode::Topology)
      .value("Harmonic", zTsSlicerMode::Harmonic);

  // One settings set for both slicers: shared entries plus T* (Topology
  // only) and H* (Harmonic only) - see SlicerSettings in Utils.h.
  py::class_<SlicerSettings>(m, "SlicerSettings")
      .def(py::init<>())
      .def_readwrite("BoundaryWidth", &SlicerSettings::BoundaryWidth)
      .def_readwrite("BracingWidth", &SlicerSettings::BracingWidth)
      .def_readwrite("BracingCurveWidths", &SlicerSettings::BracingCurveWidths)
      .def_readwrite("OverlapWidth", &SlicerSettings::OverlapWidth)
      .def_readwrite("DoubleBoundary", &SlicerSettings::DoubleBoundary)
      .def_readwrite("SlotOffset", &SlicerSettings::SlotOffset)
      .def_readwrite("BracingSlot", &SlicerSettings::BracingSlot)
      .def_readwrite("FieldScaleUp", &SlicerSettings::FieldScaleUp)
      .def_readwrite("FieldCellSize", &SlicerSettings::FieldCellSize)
      .def_readwrite("LayerSpacing", &SlicerSettings::LayerSpacing)
      .def_readwrite("_HTraceStepFraction", &SlicerSettings::_HTraceStepFraction)
      .def_readwrite("_HTraceMaxSteps", &SlicerSettings::_HTraceMaxSteps)
      .def_readwrite("_HEqualiseIterations", &SlicerSettings::_HEqualiseIterations)
      .def_readwrite("_HArapIterations", &SlicerSettings::_HArapIterations)
      .def_readwrite("_HUnrollStart", &SlicerSettings::_HUnrollStart);

  // Values for zTsSlicing.setInputMesh(mesh, face_tags) (one int per face).
  py::enum_<VertexTag>(m, "VertexTag", py::arithmetic())
      .value("Interior", VertexTag::Interior)
      .value("Bottom", VertexTag::Bottom)
      .value("Top", VertexTag::Top)
      .value("Side", VertexTag::Side);

  py::class_<PrintSynSetting>(m, "PrintSynSetting")
      .def(py::init<>())
      .def_readwrite("ResampleLength", &PrintSynSetting::ResampleLength)
      .def_readwrite("FeatureAngleTol", &PrintSynSetting::FeatureAngleTol)
      .def_readwrite("SeamOffset", &PrintSynSetting::SeamOffset)
      .def_readwrite("SeamSize", &PrintSynSetting::SeamSize);

  py::class_<ParameterGraph::Vertex>(m, "ParameterGraphVertex")
      .def(py::init<>())
      .def_readwrite("v0", &ParameterGraph::Vertex::v0)
      .def_readwrite("v1", &ParameterGraph::Vertex::v1)
      .def_readwrite("v2", &ParameterGraph::Vertex::v2)
      .def_readwrite("u", &ParameterGraph::Vertex::u)
      .def_readwrite("v", &ParameterGraph::Vertex::v)
      .def_readwrite("w", &ParameterGraph::Vertex::w);

  py::class_<ParameterGraph>(m, "ParameterGraph")
      .def(py::init<>())
      .def_readwrite("vertices", &ParameterGraph::vertices)
      .def_readwrite("edgeConnects", &ParameterGraph::edgeConnects)
      .def("empty", &ParameterGraph::empty)
      .def("clear", &ParameterGraph::clear);

  py::class_<PrintPath>(m, "PrintPath")
      .def(py::init<>())
      .def_readwrite("graph", &PrintPath::graph)
      .def_readwrite("widths", &PrintPath::widths)
      .def_readwrite("heights", &PrintPath::heights)
      .def_readwrite("normals", &PrintPath::normals)
      // Per-vertex inner(0)/outer(1) tag, parallel to widths - empty for a
      // single-perimeter path with no inner/outer distinction.
      .def_readwrite("trackIds", &PrintPath::trackIds)
      .def_readwrite("faceTags", &PrintPath::faceTags)
      .def_readwrite("edgeFaceTags", &PrintPath::edgeFaceTags);

}

} // namespace zspace_python
