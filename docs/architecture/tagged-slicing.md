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
- One `VertexTag` per **face**, in polygon order: `Bottom`, `Top`, or `Side`.
  Operators derive their vertex tags from these face tags.
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
`computePostProcess()` extracts paths and applies PrintSyn corner reconstruction,
2D resampling, seam alignment, projection and linking. FaceTag paths use
fixed-count sampling after seam alignment. Neither route resamples in 3D.
Failures retain diagnostic partial results
where the desktop implementation did so; check the return value and status
before consuming the combined path.

Per-path data includes the graph, widths, heights, normals and track IDs.
Topology intermediate outputs are available through `getTopologySlicer()`;
harmonic outputs through `getHarmonicSlicer()`. Changing mode replaces the
operator and invalidates references to its results.

### Feature targets and resampling

Production features are outer FaceTag junctions: vertices incident to two
distinct side edge tags from 3 through 6. Actual graph adjacency determines the
junction, independent of vertex storage order, angle or minimum per-point tag.
Tag 0, cap tags, slot endpoints and travel edges do not supply a side identity.
When track metadata exists, a junction must be outer (track 1); a single
perimeter is outer, while inner/bracing-controlled tracks remain 0.
`PrintPath::featureFlags` stores `OuterJunctionFeature` (bit 4) or zero at every
production target. Legacy `AngleFeature` (bit 1) and `IntersectionFeature`
(bit 2) remain for standalone compatibility but do not select features.
`featureVertexIds()` retains its angle argument for ABI compatibility and
ignores it. Reversal, seam operations, labeling, projection and linking refresh
the flags from surviving edge ownership. The added enum value does not change
the public `PrintPath` layout. Packaged runtimes require a separate deployment.

`FeatureAngleTol` uses the angle between the two incident edges: straight is
180 degrees, and an inner reconstruction/support candidate has **angle strictly
below the threshold**. Equality is excluded. Short nonzero edges participate;
coincident edges cannot define an angle. Open endpoints remain mandatory.
Outer runs do not reconstruct corners or insert targets at geometric turns.
Their dense contour supplies the reference arc for uniformly spaced targets;
chords between those targets naturally approximate that contour.

Each maximal tag0 run reconstructs independently, using supports only inside
that run. Its shared endpoints are fixed exactly, so no neighboring outer run
or junction can move. Inner corners remain ordinary targets, never yellow
junction features. Pure untagged generic paths use the same inner behavior.
`InnerResampleLength` is an independent nominal chart-arc distance, **20 mm by
default**, not a point count. Inner sampling steps by that distance between
mandatory corners/endpoints, with a terminal residual and 0.01 mm coincidence
cleanup. Counts can vary across layers. A closed inner contour always needs
at least three targets; a contour too small for three nominal intervals uses
bounded uniform subdivisions to preserve that topology.

The `ResampleLength` pass measures arc length in the **2D chart**, in model
units (normally millimetres). A shared count is computed for each outer FaceTag
across layers. Integer intervals are distributed length-proportionally across
its disconnected pieces, then each entire piece is sampled uniformly including
both endpoints. Ordinary outer angles never reset that phase; different piece
lengths and integer rounding can give different spacings between pieces.
A whole closed single-tag outer piece has a minimum of three intervals; its
topology minimum participates in the shared budget and allocation bounds.
This is chart-arc uniformity; changing surface stretch and bent chords can
produce different final3D chord lengths. Inner sampling phase restarts at its
mandatory reconstructed corners. A regular
sample within 0.01 mm of the following geometric corner or endpoint yields to
that exact point. Consecutive targets closer than 0.01 mm merge, including
across loop closure; a geometric corner wins over an ordinary seam start.
Tagged sampling first merges protected corner/category-transition spans whose
intervening chart or projected world arc is shorter than 0.01 mm. It retains
the exact surviving corner's coordinates and attributes, leaves dense ordinary
contour vertices intact, then computes the shared FaceTag budgets. A short
categorical run can disappear. Unchanged categorical ownership must permit
equal segment and unique point-label totals across layers; incompatible counts
fail explicitly. Per-span interval limits redistribute a feasible fixed count
without discarding generated targets after allocation. Final tagged targets
are checked for 0.01 mm clearance in both the chart and printable section.
Inner residuals and coincidence cleanup can differ from the nominal spacing.
USD readers convert both spacing settings into millimetres. Python exposes
`InnerResampleLength`; the slicing WASM JSON option is `innerResampleLength`.
The setting is appended to the native structure, so native consumers must be
rebuilt together even if previous field offsets are unchanged.

After sampling and seam work, mapping preserves the 2D target
counts, feature identities and order when producing 3D paths. The production
pipeline proceeds directly from projection to linking, without a second 3D
resampling pass. FaceTag budgets therefore remain intact; geometric chart
distortion can produce different spacings on the mapped 3D surface. The explicit
`resamplePrintPaths3D()` utility remains available for separate callers, but is
not part of `computePipeline()`.

Surface barycentrics update the matching 2D coordinates and interpolate existing
section heights and normals; widths and track IDs interpolate along the mapped
curve. Both path representations retain aligned attribute counts. This changes
ordinary target locations without introducing a new height computation method.
Spacing is a target on the pre-pass 3D polyline: projection and feature clearance
can change final intervals, so it is not a strict maximum segment length.

Explicit nonzero `SeamSize` still removes the cut arc, including features on
that arc, and logs the number of removed feature vertices. Surviving edge
ownership determines the refreshed masks; open endpoints cannot be junctions.
Zero-size seams remove the closing edge without deleting target coordinates.

### 2D-to-3D mapping contract and review

Each layer maps its targets barycentrically from its 2D chart to its matching
3D section. Each pair must have identical vertex order and face connectivity;
different harmonic layers can have independent topology. Mapping keeps graph
vertex count/order and connectivity, copies widths/track IDs/feature masks,
and computes normals and layer heights from the section. Source points outside
the chart or on incompatible overlapping chart triangles fail instead of
falling back to arbitrary nearest-triangle mapping. A mapping failure now
returns false even without an `outBadLayers` argument, and includes the first
failed layer/path/vertex and reason. Valid partial projections remain inspectable.

The focused native `zspace_toolsets.slicing.postprocessing` regression checks
an independent analytic piecewise-affine surface with anisotropic stretch,
a fold along the triangulation diagonal, rigid rotation and layer translation,
plus interpolated per-layer heights. It checks feature identity through both
seam-mapping strategies and the final linked targets. Full Topology/Harmonic
fixtures check serial/parallel preservation. This validates the existing
pointwise barycentric formulation. The separate `zspace_toolsets.slicing.sampling3d`
check covers the requested second pass on a folded, stretched surface, including
close mandatory features, half-spacing clearance, surface placement, aligned
attributes and open/closed curves.

Mapping remains subject to the existing containment/projection tolerances
(barycentric boundary expansion of 1e-4 and a projection-distance limit of
`max(0.1, 1e-4 * chart diagonal)` in model units). These tolerances and chart
distortion limit geometric accuracy. The second pass reduces the dependence of
ordinary spacing on 2D distortion. Edges crossing section triangles remain chords
between mapped target points; it does not insert every triangle crossing or
guarantee that entire edges lie on the section.

Both slicers retain all section/chart/height tables, including the exact top cap
as a fixed/reference mesh, but extract and print only the first `N - 1` layers.
The omitted cap contributes no printed path or print statistics. PrintSyn accepts
this printed prefix with the complete section tables; its combined path links
only printed layers.

With iso relaxation disabled, harmonic heights use the local normal-spacing estimate `deltaW / |grad(w)|`,
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

### Constrained harmonic iso relaxation

The harmonic pipeline extracts all layers, relaxes the intermediate meshes, then
unrolls them and builds the print fields. It uses the existing Eigen sparse solver
and zSpace geometry types. No libigl dependency is required. The scalar field is
piecewise linear in each tetrahedron; increasing its gradient magnitude does not
trigger geometric smoothing. Refining tetrahedra alone cannot remove a crease
already present in that field.

```cpp
auto &settings = session.getSettings();
settings._HIsoRelaxIterations = 8;       // 0 disables relaxation
settings._HIsoRelaxStrength = 0.25;      // finite, (0,1]
settings._HIsoRelaxMaxDeviation = 0.0;   // model units; 0 selects automatic limit
```

Each intermediate layer solves a screened biharmonic target with cotangent
Laplacian `L` and lumped vertex mass `M`:
`(M + strength * meanEdge^4 * L.transpose() * M.inverse() * L) X = M X_original`.
Projected, backtracked steps approach that target. Connectivity remains fixed.
Top and bottom cap coordinates are copied exactly, including genuine roof ridges.

Boundary edges retain their input side-width groups and connected planar support
patches. A boundary vertex can slide within its patch; at two patch junctions it
slides on the shared input edge chain; a junction of three or more patches stays
fixed. Entire boundary segments must remain covered by their input patch, so
sliding cannot cut a chord across a folded wall. This conservatively retains the
existing relationship to the input triangulation. A missing source width tag can
be recovered only through unique full-segment support, with a diagnostic count.
The surface/containment tolerance is `max(1e-9, 2e-7 * input diagonal)` model units,
with barycentric tolerances scaled by simplex altitude.

Total vertex displacement is limited to 20% of the initial minimum triangle-to-
triangle gap to either adjacent layer. A positive `MaxDeviation` further reduces
that limit. Each accepted step checks volume containment, full boundary support,
triangle reversal/collapse, intersections with the input surface, self-overlap
and adjacent-layer intersections. Float-scale sliver vertices are held fixed;
vertices without a valid transported UV remain at their previous position.
Rejected steps backtrack and eventually stop. `getIsoRelaxationReport()` exposes
accepted iterations, constrained/moved vertex counts, displacement, boundary error
and the last rejected-step reason; inspect it when constraints limit smoothing.
Local geometry violations activate fixed vertices rather than stopping unrelated
regions. Initially planar sections stay in their fitted plane, preventing float
sliver conditioning from introducing a spurious bend; the solve uses centred
coordinates. Triangle-intersection predicates use a tighter `0.001 * surface
tolerance`, separate from the expanded projection/containment slab.

Moved vertices receive fresh tetrahedral UV coordinates. Intermediate print normals
are recovered from area-weighted mesh normals and barycentrically interpolated
onto path targets. Unsupported float-scale normals are recovered from a reliable
neighbouring mesh ring. Fixed caps retain their outward winding and solver
normal convention. For fabrication only, the bottom cap supplies inward
area-weighted vertex normals to path projection, so its one-sided print height
extends toward increasing harmonic layers. This local direction follows the
cap geometry even when the solid is rotated; source meshes, tetrahedra and
height-solver normals are unchanged. The projection override is rebuilt after
relaxation and surface-height feedback, including zero relaxation iterations.
Heights are recomputed against neighbouring relaxed meshes along
the recovered normal, averaging the two full adjacent gaps for an intermediate
layer and using one full gap at each cap. If a trimmed side boundary prevents a
ray hit, a local tangent-plane estimate from the closest adjacent point and its
recovered normal is used and counted. These are geometric spacing estimates,
not a nozzle-clearance or extrusion optimisation.

Outer junction features follow the resulting path's edge ownership. Geometric
inner turns remain ordinary supports; outer turns do not add targets. PrintSyn carries aligned junction
metadata through projection and linking. Relaxation does not move already-generated path targets.
The optional per-section normal tables follow the same ordering contract as the
per-section height tables; an empty normal table uses the existing face-normal
mapping. Zero relaxation iterations retains the previous geometry, gradient-based
heights and face-normal path behavior. The requested strict previous-mesh-only
height rule was withdrawn; the existing height formulations described above
remain in use. The controls currently exist in C++. The local TASK-008 Alice2
build was rebuilt against this source; packaged Python/browser consumers still
need a coordinated rebuild and their runtimes have not been refreshed.

Relaxation caches triangle geometry and bounds, computes exact minimum layer
separation with a best-first pair-of-BVHs search, and runs independent per-layer
checks and updates in parallel when requested. Fixed-size traversal/contact
buffers remove allocations from repeated queries. These changes retain the same
constraints and targets. `zTsRelax::Settings::Parallel` follows the slicer's
`useParallel` argument. `Report::StageMilliseconds` records setup, gap bounds,
sparse targets, trial checks, normals/heights plus diagnostics, and total time.

On the shared `slicer_inputs/example.usda` (SHA-256
`9BE2169D3779878DC42CAD2EA95F81E5B232C1F26C2AA588F2761C73563D4E52`),
Release/MSVC, four OpenMP threads, one warmup per mode and five alternating
measured runs gave the following medians. Input parsing, rendering, export and
build time are excluded; all other functional settings match.

| Measurement | Relaxation disabled | Eight relaxation iterations |
| --- | ---: | ---: |
| Complete native pipeline | 6.340 s | 9.484 s |
| Contour stage, including relaxation | 0.148 s | 3.327 s |
| Relaxation alone | 0 s | 3.170 s |
| PrintSyn post-processing | 2.306 s | 2.296 s |

The pre-optimization profile measured relaxation alone at 27.905 s: 16.242 s
for gap bounds, 9.195 s for trial checks and 2.063 s for normal/height refresh.
The optimized median is about 8.8 times faster, with a remaining complete-pipeline
overhead of 3.144 s (49.6%) versus disabled relaxation. Optimized iso positions
and section height tables are bit-identical to the captured pre-optimization
result. All 2,339 retained angle/intersection targets on the example survive the
new sampling pass exactly. There are 64 harmonic reference sections/63 printed
layers and 66 topology sections/65 printed layers.

Actual Alice2 F1 captures confirm both pipelines run on this USD. The disabled
harmonic gradient-based height method produces extreme finite heights near the
top (maximum about `3.41e13` mm), making that unrelaxed print-mesh view unusable;
its path view remains inspectable. This pre-existing height behavior is retained
at the user's request and is not a valid physical extrusion result. Relaxed
example path heights range from 5.841 to 16.983 mm. Native checks passing do not
constitute printability approval.

`zTsUtils::buildPrintMesh()` lofts one rectangle per print-path vertex, with width
centered on the path (`-width/2` to `+width/2`) and height extending from `0` to
`height`. Sections are perpendicular to the neighbor-to-neighbor chord, with
the print normal projected into the section's height axis. Raw edge vectors
retain segment-length weighting at corners; endpoint tangents use the adjacent
edge and closed curves use wrapped neighbors. Chords stay forward throughout to
avoid flipping rectangle vertex order. A complete reversal uses the outgoing
edge; a print normal parallel to the section tangent is rejected.
Each graph component and
layer stays separate; open paths receive end caps. Zero widths/heights retain
collapsed sections, and invalid attributes or connectivity clear the mesh and
return a diagnostic with the layer/path location.

USD export preserves `/TDP_gcode` and also writes `/TDP_printmesh` (Mesh, display
name `print mesh`) and `/TDP_printpath` (linear BasisCurves, display name `print
path`). Pass the processed per-layer paths through the optional `layerPaths`
argument of `zTsIO::exportPrintPathUSD()` to avoid lofting travel connectors
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

Harmonic tetrahedral density is controlled in C++ by
`session.getSettings()._HTetEdgeLengthScale` (finite, in `(0,1]`). The default
`0.5` uses a target edge length of half the input surface's mean triangle
edge, hence a maximum volume target of `0.5^3 = 0.125` times the former
bound. `1.0` restores the original density rule. TetGen quality constraints
also influence actual counts; the block regression increases from 50 to 122
tetrahedra at `0.5`, whereas a `0.75` factor does not increase its density.
The input surface positions/connectivity and the other slicer settings stay
unchanged. TetGen may split its internal boundary facets while preserving the
same piecewise-planar domain and face markers. Recompute the harmonic field
and slices after changing the factor. This C++ control is not yet exposed as
a desktop/viewer slider or a Python property; downstream binaries must be
rebuilt to use the new default.

Cap extraction retains geometrically valid faces even when a vertex's UV
trace failed. UV validity controls the unfolding start/fallback, rather than
whether the 3D cap exists. This prevents added tetrahedral boundary vertices
from opening cap holes solely because their UV traces failed. The native
regression checks this explicitly, along with refinement counts/volumes,
source-mesh preservation, unchanged slicing settings and invalid factors.

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

For the optional desktop USD adapter, set
`ZSPACE_TOOLSETS_SLICING_WITH_OPENUSD=ON` and `ZSPACE_OPENUSD_ROOT` to an external
OpenUSD SDK. The adapter reuses core's OpenUSD CMake helper. Deploy the SDK's
runtime libraries and schema/plugin resources with the application. Its public
API is `zTsIO::readSlicerInputUSD()` / `exportPrintPathUSD()`; without the option,
these functions return a diagnostic explaining the missing capability.

The USD adapter preserves `/block`, `/bracing`, `/slotpoint`, face/width tags,
settings and print-path attributes. It is separate from the compute engine.

## Browser boundary

The dedicated `wasm/slicing` target compiles the same engine and PrintSyn sources
into `wasm/out/zspace_slicing.js` and `.wasm`. It leaves the established core and
toolsets runtimes unchanged. Run `wasm/scripts/build-slicing-wasm.ps1` with the
Emscripten SDK installed. The default TetGen source directory is
`wasm/cache/tetgen`; provide `-TetgenDir` for another source checkout. The reviewed
build uses Emscripten 6.0.13 and libigl/tetgen commit
`e05aca7df74e3f531bc35733ed87d36d437266c5`. Without TetGen, capability bit 2 is
absent and harmonic mode is rejected. OpenMP is disabled in this worker build.

The exported ABI has version 1 and capability bits 1 (topology) / 2 (harmonic).
Create/destroy manage handles; `zspace_slicing_compute_from_buffers` uploads
owned mesh arrays, face tags, optional width groups and bracing, anchor, mode and
settings using the existing ZST1 binary codec. Read results with
`zspace_slicing_result_buffers` and `zspace_toolsets_result_buffer_size`. Input
allocation uses `zspace_alloc`/`zspace_free`. Errors use
`zspace_toolsets_last_error`. Failed computation preserves the previous result.

Results contain owned slice/chart meshes, path/width fields, per-layer print
paths and the combined PrintSyn path with widths/heights/normals/track IDs. The
viewer computes in a module worker and copies the result before releasing its
native session. No live heap views reach viewer scene objects.

Browser USD IO is owned by the core OpenUSD runtime. Toolsets WASM excludes
`utils/IO.cpp` and accepts only owned geometry, tags, width groups and settings.
The viewer translates `/AG_block`, `/Slicer_bracing`, `/Slicer_slot`,
`/Slicer_data` and `/Slicer_setting` through generic core USD queries and converts
lengths to millimetres. See the viewer's `docs/development/slicing.md` for builds
and tests. Native desktop IO remains available independently.
