# zTsForcePolyhedra

Graph-based disjointed force polyhedra and a geometric 3D graphic-statics
reciprocation solver. The input graph is the **form diagram**: vertices are joints,
edges are axial members. Each joint owns an independent force-cell mesh. No
volumetric mesh or global topological dual is required, so chains, branching,
crossing members and arbitrary graph connectivity use the same API.

## References and procedure

Juney Lee, Tom Van Mele and Philippe Block, *Disjointed force polyhedra*,
Computer-Aided Design 99 (2018), 11–28, DOI 10.1016/j.cad.2018.02.004,
especially sections 3.2–3.4 and Figures 6–9. The implementation is independently
written C++ using the supplied paper as reference.

The geometric equilibrium procedure is reimplemented from
`legacy/pre-v3/cpp/source/zToolsets/statics/zTsPolytopal.cpp`:
`computeFormTargets`, `updateFormDiagram`, and `checkParallelity`. It extends
that form-only procedure to disjoint cells and movement of either diagram.
It does not use force-density form finding or an algebraic stiffness solve.

1. Snapshot the graph through `zFnGraph`; assign supports, loads and optional
   target member forces. Positive axial forces are tension; other members use
   compression unless explicitly assigned with `setTensionEdges`.
2. Build the EGI minor-arc arrangement. Split arcs at cross-adjacencies and
   insert zero-area EGI masses. Extract spherical regions and their dual mesh.
3. Add untargeted virtual caps for planar or collinear nodes. A temporarily
   unbalanced node receives a zero-area closing face, as in paper Figure 6c.
   Zero and virtual faces carry no physical force or graph member.
4. Build face goals from bounded normal rotation, planarity and area scaling.
   Area scale goals are limited to 0.8-1.25 per pass. Average incident goals
   into correction forces and integrate them with `zFnMeshDynamics`. Zero faces
   receive principal-line collapse goals. Polygon vertex order is established
   at seed creation and is never permuted during an iterative update.
5. Run 3DGS reciprocation. Read both corresponding outward normals per member,
   account for compression/tension, and blend with the graph direction using
   the movement weight. Neighbour-based form corrections are integrated by
   `zFnGraphDynamics`, retaining fixed supports and the legacy length clamp.
   Both corresponding normals and their mutual opposition are checked.

Both function sets own persistent particle systems on stable, separate backing
objects. Unit-mass Euler integration is overdamped: force, velocity and derivative
state are cleared after each step. Correction increments use rebased local
coordinates, then accumulate into double precision diagram positions, avoiding
float stagnation in shallow cells and distant world coordinates. No direct
least-squares vertex projection follows the dynamics step.

`timeStep` (0-1), `maxVertexStepRatio` (0-0.25) and `normalStepDegrees` (0-30)
control the update. Positive values are required. Default vertex motion is at
most 0.05 times the initial cell extent or longest form edge, per integration
step; a common force scale preserves relative correction directions when this
bound is active. Changing force directions does not relabel cell vertices or
reverse virtual cap identities. This bounded dynamics extension is an
implementation choice, not a verbatim reproduction of the paper prototype.

An adaptive time-step acceptance check measures relative area error, oriented
normal error, planarity, convexity and zero-face line-collapse error. Trial
increments are halved if they worsen this combined merit or rotate a face more
than the permitted angle. If the averaged face goals do not provide descent,
the same dynamics engine integrates the merit gradient instead. A failed
trial leaves the previous geometry intact; it never commits an inverted face
or relabels vertices. This is a numerical stability check, not a declaration
that the physical problem is infeasible.

## Diagram movement weight

`formWeight` is a **movement weight**:

| Value | Form graph | Force cells |
| --- | --- | --- |
| 0 | Fixed exactly | May move |
| 0.5 | Free vertices may move | May move |
| 1 | Free vertices may move | Fixed exactly |

Support vertices remain pinned. Applied load vectors, including self weight,
retain their prescribed directions and magnitudes. At weight 1, AP is not run
during reciprocation; to finish AP first, call `stepAreaPursuit`. Changing only
the movement weight preserves the current diagrams.

At weight 0 the form graph is fixed and only mesh particles move. Member area
targets redistribute gradually from the mean paired actual face contributions.
Applied load areas and optional prescribed member magnitudes remain fixed.
There is no force-density solve or fixed-form feasibility gate. Incompatible
constraints can leave errors nonzero; iteration alone is not equilibrium.

Reports identify the worst angle edge, worst area cell, largest moving cell,
maximum form/cell vertex steps, and per-cell dynamics step counts. The area,
normal/planarity, angle, contact-area and physical force tolerances must all
pass before `equilibrium` becomes true.

## C++ use

```cpp
#include <zspace/toolsets.h>
using namespace zSpace;

zTsForcePolyhedra gs;
gs.setInputGraph(graph);
gs.setSupports({0, 4});
gs.setVertexMasses({0, 1, 1, 1, 0});
gs.setGravity(zVector(0, 0, -9.81));
// Optional explicit loads are added to mass * gravity.
// gs.setVertexLoads(loads);
// Optional force-driven constraints, one positive magnitude per original edge:
// gs.setEdgeForceMagnitudes(targetMagnitudes);

gs.createForcePolyhedra(); // EGI + AP; form graph coordinates are unchanged.
gs.setFormWeight(0.5);
auto report = gs.solveEquilibrium(2000);
// Or gs.stepEquilibrium(1) for playback.
const auto& form = gs.resultGraph();
const auto& forceCell = gs.forceCell(2);
const auto& faceMetadata = gs.cellFaces(2);
```

Without specified member magnitudes, a load-scaled positive seed initializes
each member. During reciprocation the two measured contact-area contributions
are averaged into a common area target, permitting force redistribution along
loaded chains. Explicit member magnitudes remain fixed instead.

`setFormPositions` edits the form while retaining the force cells, useful for
force-fixed reciprocation. Load/support/mass/sign/target edits invalidate the
force cells. `reset` restores the input graph and clears computed diagrams while
retaining assignments. Input graph replacement resets the entire session.

Use consistent units. `P_i = explicitLoad_i + mass_i * gravity` includes loads
applied directly at supports. Mass defaults to zero and gravity to `(0,0,-1)`.
Support reactions include the support's own loads; unsupported nodes never
receive artificial reactions. Force-cell area is magnitude times `forceScale`.
Cells are centred locally at the origin; place them near their form nodes for
display. Virtual cap sizes are visual, with no physical interpretation.

## Mapping and convergence

`cellFaces(vertexId)` follows the corresponding mesh polygon order. Each face has
a target normal, target/measured areas, a type and physical-source records.
Members refer to original `zFnGraph::getEdgeData` edge IDs; loads and reactions
refer to original vertex IDs. Equal-direction forces share a convex face, with
all source contributions retained. In that case each source's share of the total
face area is used for contact matching; a merged face is not falsely assigned
in full to each member.

`equilibrium` requires all of:

- AP area, orientation, planarity and supporting-plane checks pass for every cell;
- every form edge is parallel to both corresponding signed face normals, and
  the two contact normals oppose each other within `angleTolerance` (degrees);
- the corresponding contact-area contributions differ by no more than
  `correspondingAreaTolerance` (absolute force-diagram area units);
- maximum physical free-node force residual is below `equilibriumTolerance`.

Reports separate `maxAngleDeviation`, `maxContactAreaMismatch`, `maxAreaError`,
`maxNormalError`, `maxPlanarityError`, and physical `maxResidual`/`residualNorm`.
`iterations` counts AP work and `equilibriumIterations` counts reciprocation
steps. An incomplete numerical solve is returned for inspection and continuation.
This does not assess structural rigidity, buckling or material response.

Limits: 1024 graph vertices, 2048 graph edges, 16 distinct physical directions
per cell and 2048 EGI directions after splitting. Self loops, duplicate members,
collapsed edges, invalid IDs and nonfinite inputs are rejected.

## WASM

The shared toolsets build includes `wasm/bridge/zspace_force_polyhedra_bridge.cpp`.
Force Polyhedra sessions have their own lifecycle and error getter:

- `zspace_forcepolyhedra_create` / `zspace_forcepolyhedra_destroy`
- `zspace_forcepolyhedra_compute_from_buffers(handle, pointer, bytes)`
- `zspace_forcepolyhedra_step(handle, iterations)` for AP
- `zspace_forcepolyhedra_equilibrium(handle, iterations)` for reciprocation
- `zspace_forcepolyhedra_set_form_weight(handle, weight)`
- `zspace_forcepolyhedra_set_form_from_buffers(handle, pointer, bytes)`
- `zspace_forcepolyhedra_result_buffers(handle)` and
  `zspace_forcepolyhedra_result_buffer_size(handle)`
- `zspace_forcepolyhedra_last_error`

Use the existing ZST1 typed-transfer contract and shared `zspace_alloc` /
`zspace_free`. Reacquire heap views after allocation; free uploads in `finally`;
copy result bytes before the next operation on the session. Input is an object
with `graph: {vertexPositions, edgeConnects}`, `supports`, optional `masses`,
`mass`, `gravity`, `loads`, `tensionEdges`, `edgeForceMagnitudes`, and `settings`.
Compute runs EGI/AP; optional `solveEquilibrium: true` also runs reciprocation,
with an `equilibriumIterations` budget. Form updates contain `vertexPositions`.

The result includes form graph geometry, each cell's polygon arrays and source
mapping, signed member forces, reactions, residuals and all reports. Failed
input computations retain the previous session. Legacy JSON file/string
compute/result exports are also available for parity checks.

## Verification and scope

```powershell
cmake --build build/native --target zspace_toolsets_force_polyhedra
ctest --test-dir build/native -R force_polyhedra --output-on-failure
powershell -NoProfile -ExecutionPolicy Bypass -File wasm/scripts/build-wasm.ps1
node tests/force_polyhedra_wasm.mjs
```

Native/WASM tests cover analytic force cells, loaded chains, branching, virtual
faces, irregular EGI cross-adjacencies, exact diagram freezing at the weight
endpoints, geometric reciprocation, mapping/export consistency, invalid input,
typed/legacy parity and independent sessions. The toolset and bridge require no
core changes. Alice provides the `Force Polyhedra` sketch with Rhino Directrix
import, equal arc-length sampling, support picking, area pursuit and weighted
3DGS controls. Its timeline advances bounded AP/3DGS steps, caches result frames,
and deterministically replays native operations when branching after rewind.
Reset recreates the initial cells on the imported graph. The sketch supplies
zero mass at supports, hides their boundary cells, mirrors displayed free-node
cells for the upward load-line convention, and displays geometric support
thrust as the negative mean outgoing edge vector. The paper's unified-diagram
visualization is not included.

Low-valence cells use virtual caps. Initial three-force cells are seeded as
closed force triangles; the subsequent reciprocity iterations update their
directions freely. Virtual caps follow the physical force plane and have no
physical load contribution.

`areaPursuitResolved` reports only whether every corresponding member-area pair
passes `abs(A0-A1) <= correspondingAreaTolerance` (default `1e-4`). It is separate
from `cellsConverged`, which includes prescribed target areas, normals and
planarity. `maxCorrespondingAreaDifference` and `worstAreaDifferenceEdge` identify
the limiting pair. `maxContactAreaMismatch` retains its relative diagnostic.
Merged faces use source-weighted contributions; display scale never affects
these areas. Equilibrium requires both criteria.

The WASM `zspace_forcepolyhedra_set_controls` operation applies corresponding-area,
angle, movement-weight and normal/planarity tolerances through `setSettings`
without reconstructing meshes, graphs or dynamics. The viewer retains geometry,
frame and history on these edits, reclassifies cached frames under active
tolerances, and records the controls used for each numerical step for replay.
