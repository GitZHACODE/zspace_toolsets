# Tagged slicing and PrintSyn

The tagged slicing engine is available through `<zspace/slicing.h>` as
`zSpace::zTsSlicing`. It ports the topology/harmonic slicers and PrintSyn from
the desktop `slicer` branch (`2db24962`) into the maintained v3 source layout.
The earlier v3 3DP prototype has been removed. `zTsSlicing` is the maintained C++ slicing API.

Public headers live in `include/zspace/zToolsets/fabrication/slicing/` and
implementations in `src/zToolsets/fabrication/slicing/`. Headers no longer include
implementation `.cpp` files; link against `zspace::toolsets`.

## Input contract

- A polygon mesh, copied into the session. Topology mode requires all-quad faces,
  one connected floor patch and regular columns connecting the floor to the top.
- One tag per **face**, in polygon order: 1 Bottom, 2 Top, or 3–6 side
  regions. Operators derive Bottom/Top/Side vertex roles while preserving
  each original side-region ID independently of its width ID.
- Optional per-face width-group IDs, indexing `SlicerSettings::BoundaryWidth`.
- An optional bracing graph in the input mesh's coordinate system.
- An explicit slot anchor in that coordinate system, required by both operators.
- Slicing and PrintSyn settings. Coordinates, layer spacing and widths use the
  same model units; the compute engine does not convert units.

`zTsSlicing` owns copies of the mesh, face tags, width IDs, bracing graph and
anchor. The lower-level operators borrow their inputs. Keep those inputs alive
when using `TopologySlicer` or `HarmonicSlicer` directly.

## Minimal usage

```cpp
#include <zspace/slicing.h>
#include <stdexcept>

zSpace::zTsSlicing session;
session.setInputMesh(mesh, faceTags, faceWidthIds);
session.setBracingGraph(bracing); // optional
session.setSlotAnchor(anchor);
session.getSettings().LayerSpacing = 1.0f;
session.getSettings().BoundaryWidth = {0.5f};
session.getSettings().FieldCellSize = 0.05f;
session.getPrintSynSetting().ResampleLength = 0.1f;

if (!session.computeSlicer(false))
    throw std::runtime_error(session.getStatus());
if (!session.computePostProcess(false))
    throw std::runtime_error(session.getStatus());

const auto& layers = session.get3DPrintPaths();
const auto& combined = session.getResultantPrintPath();
```

`computeBasicStructure()` validates/builds the initial structure.
`computeSlicer()` runs through layer meshes, unrolled charts and fields.
`computePrintPath()` extracts/maps paths without PrintSyn.
`computePostProcess()` extracts paths and applies PrintSyn resampling, seam
alignment, projection and linking. Failures retain diagnostic partial results
where the desktop implementation did so; check the return value and status
before consuming the combined path.

Per-path data includes the graph, widths, heights, normals and track IDs.
Inputs containing side tags 4–6 also enable FaceTag sampling. PrintSyn prepares
corners and orientation, aligns/cuts the seam, then samples the tagged outer
segments before projection/linking. Legacy inputs with only side tag 3 retain
their existing sampling order and behavior.

For each FaceTag, the target interval count comes from its mean **flattened 2D**
printable outer length across all layers divided by `ResampleLength`, rounded to
the nearest integer. Open-piece endpoint counts and required corner anchors set
the point-budget minimum. Every layer uses the same total point budget for that
tag. A staggered slot's disconnected pieces share that budget in proportion to
their lengths on that layer; their individual counts can change. Spacing is
uniform between protected anchors, subject to integer allocation. No outer
segment bridges the slot; its returns are excluded from the outer budget.

The shape/width fields remain unchanged. A compact continuous field provenance
query uses the same width-group boundary distances and active CSG branches to
identify regions, including bracing and slot returns. Category transitions
become protected contour anchors; integer IDs are never linearly interpolated.
`PrintPath::faceTags` carries a side ID on each outer sample, or 0 for another
role. `edgeFaceTags` labels graph edges. `faceTagSegments` records ordered,
endpoint-inclusive region membership. Shared corners are stored once and can
belong to both incident segments; their single point label is the smaller
nonzero incident tag. Therefore segment budgets include shared endpoints while
per-point label counts assign each corner once. Both counts stay consistent
across layers with the same region topology. Missing tags fail explicitly.
Raw contour labeling preserves geometry. Category anchors are inserted after
corner/seam preparation; consecutive anchors that project to the same physical
point merge before budgeting. This prevents zero-length 3D edges while retaining
the shared point count and coincident feature identities.
Projection preserves these indices and labels; linking assigns travel edges tag
0. USD `/TDP_gcode.printFaceTags`, WASM result slot 6 and Python `faceTags`
expose the per-point labels. `trackIds` continues to mean inner 0 / outer 1.
Topology intermediate outputs are available through `getTopologySlicer()`;
harmonic outputs through `getHarmonicSlicer()`. Changing mode replaces the
operator and invalidates references to its results.

Harmonic heights use the local normal-spacing estimate `deltaW / |grad(w)|`,
where `deltaW = 1 / (numberOfLayers - 1)` and `w` is the final harmonic field
(including any field equalisation). Tetrahedron volumes weight the nodal average;
heights are then interpolated on tetrahedral edges and barycentrically onto each
layer's paths. This is a local linear estimate, not an exact distance to the
previous curved section or a bound enforced by `LayerSpacing`. Both cap layers
use a full adjacent iso interval, matching topology's nonzero cap-height contract.
Zero/nonfinite gradients or unsupported vertices fail contour generation instead
of silently producing zero heights. `HarmonicSlicer::getVertexLayerHeights()`
returns one table per iso mesh, in that mesh's vertex order. Raw paths and PrintSyn
projection use these tables; linked paths retain the heights, so the existing
height-domain and print-volume statistics work in harmonic mode. Settings/input
invalidation clears the tables and statistics. Target-domain optimization remains
a separate operation.

`zTsUtils::buildPrintMesh()` lofts one centered rectangle per print-path vertex,
using its width, height and normalized print normal. Each graph component and
layer stays separate; open paths receive end caps. Zero widths/heights retain
collapsed sections, and invalid attributes or connectivity clear the mesh and
return a diagnostic with the layer/path location.

USD export preserves `/TDP_gcode` and also writes `/TDP_printmesh` (Mesh, display
name `print mesh`) and `/TDP_printpath` (linear BasisCurves, display name `print
path`). Pass the processed per-layer paths through the optional `layerPaths`
argument of `zTsIO::exportUSD()` to avoid lofting travel connectors
between layers. Without it, the supplied path alone is rendered. Geometry and
widths are converted from internal millimetres to the source stage units;
re-export replaces these output prims while preserving the input stage.

## Build configurations

The default native build includes topology slicing and PrintSyn without TetGen,
OpenUSD or OpenMP. Enable tests with `ZSPACE_TOOLSETS_BUILD_TESTS=ON`.

```powershell
cmake -S . -B build/slicing -DZSPACE_TOOLSETS_BUILD_TESTS=ON
cmake --build build/slicing --config Release --parallel 4
ctest --test-dir build/slicing -C Release --output-on-failure
```

Run from a compiler-enabled shell, using the same generator/toolchain as the
rest of your v3 build. The native smoke suite includes topology/PrintSyn,
input validation, width groups/bracing/boundary slots, serial/parallel parity,
repeat computation, and moved-session reuse, alongside the
existing v3 3DP regression.

For harmonic slicing, build/provide TetGen for the same platform and runtime:

```powershell
cmake -S . -B build/slicing -DZSPACE_TOOLSETS_WITH_TETGEN=ON `
  -DZSPACE_TETGEN_ROOT=C:/SDKs/tetgen
```

The root should contain `include/tetgen.h` and `lib/tetgen.lib` (or the platform
equivalent). Separate installations can set `ZSPACE_TETGEN_INCLUDE_DIR` and
`ZSPACE_TETGEN_LIBRARY` explicitly. Nothing is downloaded or vendored by this
option. Harmonic sources and their regression are included only when enabled.
Query `zTsSlicing::supportsSlicerMode()` before selecting a mode. Selecting an
unavailable or invalid mode throws `std::invalid_argument` and preserves the
current session.

`ZSPACE_TOOLSETS_SLICING_OPENMP=ON` enables existing parallel loops. Sequential
calls (`useParallel=false`) remain available in either build. Eigen stays sequential
to avoid nested parallel scheduling.

USD IO belongs to Core's `zUSDsdk`, enabled with `ZSPACE_IO_WITH_OPENUSD`
and `ZSPACE_OPENUSD_ROOT` in the host build. Toolsets has no USD option, IO
source, or `zspace::io` link dependency. The former `zTsIO` API is removed.
Alice's `userSrc/zspace/sketch_harmonic_slicer_test.cpp` contains the
application schema reader and writer and calls Core IO directly. Python exposes generic IO through
`zspace.core.USDDocument` and schema helpers through `zspace.io.slicing_usd`.
Deploy Core IO and the matching OpenUSD runtime/plugins only in hosts needing IO.

All hosts use `/AG_block`, `/Slicer_bracing`, `/Slicer_slot`, `/Slicer_setting`
and `/Slicer_data`. Geometry and authored dimensional settings become Z-up mm;
ratios remain unchanged. Optional settings use internal mm defaults. Exports
preserve input and replace `/TDP_gcode` (Scope), with custom `points`,
`printNormals`, `printWidths`, `printHeight` arrays and root `customLayerData`
statistics. Source units and up axis are restored on export. Applications own
this schema; Core owns only general USD operations.

## Browser boundary

The dedicated `wasm/slicing` target compiles the same engine and PrintSyn sources
into `wasm/out/zspace_toolset.js` and `.wasm`. It leaves the established core and
toolsets runtimes unchanged. Run `wasm/scripts/build-slicing-wasm.ps1` with the
Emscripten SDK installed. The default TetGen source directory is
`wasm/cache/tetgen`; provide `-TetgenDir` for another source checkout. The reviewed
build uses Emscripten 6.0.13 and libigl/tetgen commit
`e05aca7df74e3f531bc35733ed87d36d437266c5`. Without TetGen, capability bit 2 is
absent and harmonic mode is rejected. The runtime supports serial and OpenMP
execution with a hardware-sized Worker pool.

The Slicer ABI is version 2, with capability bits 1 (Topology) / 2 (Harmonic).
Create/destroy manage handles. Numeric array pointers/counts carry geometry,
settings, scalar fields and path attributes without JSON or a tagged packet.
See [Slicing array ABI](slicing-wasm-arrays.md) for slots and lifetime rules.
Core supplies allocation and its own unchanged bindings in the combined runtime.
Failed computation preserves the previous result.

Results contain owned slice/chart meshes, path/width fields, per-layer print
paths and the combined PrintSyn path with widths/heights/normals/track IDs.
The adapter copies native arrays into owned Float64Arrays and transfers their
buffers to the UI. Native point order is preserved without sorting.

Toolsets IO parses the Slicer USD schema through Core's generic USD API.
Both import and export run in an IO Worker using the same combined runtime.
JavaScript transfers file bytes and numeric arrays; it does not duplicate
the schema interpretation. Geometry and dimensional settings normalize to mm.
