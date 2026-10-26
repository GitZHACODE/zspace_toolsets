# Parametric Stairs

`zTsParametricStairs` follows the modern C-Tubes header/source layout and uses
`zTsVault` for Best Fit FDM. Rhino SDK dependencies stay outside the toolset.
The Alice sketch imports straight line curves from the case-insensitive `graph`
layer using its existing Rhino kernel and line-joining helper. Shared endpoints
are welded with the join tolerance; endpoint-on-edge contacts split an edge.
Interior crossings are not automatically connected. Each component must have
at least two degree-one ends.

## Construction and defaults

`update(graph, settings)` first builds a coarse continuous quad ribbon: two
faces per graph edge, with only boundary edges, central graph edges, and the
cross-edges at graph nodes. There are no longitudinal sample rows at this stage.
Adjacent arms share centre and angular
miter vertices at junctions; there are no disconnected overlapping branch caps.
Widths may vary per graph vertex (`vertexWidths`, full widths in model units).
Junction offsets use XY
directions and preserve each graph node's elevation. This is an XY layout with
variable Z, including branched stair/bridge forms from the supplied workflow.

`numDivs` inserts transverse cross-sections along each graph arm, without adding
extra lengthwise lanes across the width. Each level doubles the longitudinal
segments: zero retains the coarse topology, one inserts a midpoint section,
and two gives four segments per arm. The original central graph remains the
only interior lengthwise lane before smoothing. Shared junction vertices and
the three vertices on each leaf cap remain intact at this stage.
`numSmoothDivs` then applies core Catmull-Clark subdivision, **before** FDM.
Both default to one. The smoothed mesh becomes the solver input and attraction
reference; FDM fixes that mesh's XY coordinates by default. Smoothing can move XY before
the solve. `lowPolyMesh()`, `subdividedMesh()` and `centerMesh()` expose the three
construction stages. Original cap supports stay pinned during smoothing; new
vertices on each cap inherit support membership. Explicit legacy `divisions`
inputs retain the previous longitudinal sampling behavior.

Every centre and side vertex across a valence-one graph end becomes a support.
All edges default to compression, with fixed XY and gravity (0,0,-1), uniform
positive vertex mass, and bounded force densities. Support elevations stay
fixed. The working mesh reference is used for attraction; FDM updates that mesh
and reports residuals, convergence, rank, and stopping status.
Initialisation fits force densities; `stepBestFit()` advances form finding.
An above-tolerance best fit is not reported as equilibrium.

Setters regenerate and initialise a complete candidate before replacing the
previous output. Sharp miters, collapsed/folded strips, invalid widths and
unsupported topology throw while retaining the previous valid result.
Mapping arrays expose graph-node centre vertices and the source graph vertex
for each retained graph-station vertex and refined cap support (-1 for other
refinement vertices). Centre vertices coincide with graph nodes in the coarse
mesh and can move during smoothing.

## Viewer workflow

Choose **Parametric Stairs** and load it.
The default input is `C:/Users/vishu.b/Desktop/Geometries/Stairs_Graph.3dm`.
When available, the sketch reads this document automatically. If the file cannot
be read, it uses the generic Y graph; **Reload Rhino file** retries the same path.
The input text remains editable. Invalid geometry in a readable file is reported
rather than silently replaced.
The Rhino file's complete layer tree is imported beneath a filename group,
including Landing, hidden and empty layers, colors, and object visibility.
Curves on Graph become the editable central graph; other supported objects are
tessellated as reference geometry in their original layers. They are reloaded
with the source file and do not enter the Stairs/FDM mesh. Unsupported reference
objects are reported while other layers continue loading.
Enter a local `.3dm` path to import its `graph` layer, or select an existing
zObjectGraph and use **Use selected graph** once. Press **Generate base mesh
from graph** to construct the mesh and enable live rebuilding. The generation
button also accepts the currently selected graph directly. Loading or assigning
a graph alone does not generate a mesh. Edit the published central
graph in vertex mode or transform it. Width, end-width and junction-width
controls regenerate the mesh automatically; no Compute action is needed.

The sketch checks the assigned graph in world coordinates by persistent object
ID every 100 ms, independent of selection. Parameters are read live. This
watcher updates geometry only. Use the existing bottom **Play/Pause** controls
to advance Best Fit FDM, until convergence, stall, or iteration limit. Step,
rewind and reset use the same transport; reset restores the smoothed input.
Loading this sketch selects its solver timeline automatically. Changing source
geometry or division settings stops playback and resets the solve to the new
reference. Failed
geometry updates preserve the last valid output and pause fitting until the
input changes. Removing the graph clears derived output. Switching sketches
suspends updates; disposal clears the timer and destroys the WASM handle.

The bottom transport dispatches to this sketch timeline before ordinary mesh
solver playback, including when a generated mesh with no solver nodes is selected.
**Snap FDM (fixed densities)** offers a direct sparse equilibrium solve, without
animated outer iterations. It holds the current fitted, bounded densities and supports fixed, respecting
the optional XY constraint. **Geometry attraction (Snap)** is a soft reference penalty:
0 solves equilibrium in the free axes; 1 restores them to the reference.
Between them it minimizes `||L z - b||² + lambda ||z - zReference||²`, with
`lambda = mean(diag(LᵀL)) * weight / (1-weight)` (a small scale floor handles zero
densities). This is an attraction weight, not a maximum-displacement constraint.
More attraction can leave more force residual. Full XYZ residuals, including
horizontal imbalance at fixed XY, remain the convergence criterion. This is
not a simultaneous global optimum over geometry and force densities. Repeated
snaps at unchanged density/weight solve the same system rather than animate.
Reset and rewind retain their history behavior. **Geometry update share
(iterative)** separately preserves the existing Best Fit FDM update weighting.

After Play or Snap, **Smooth after form finding** refines the current working
mesh in place. Existing supports remain pinned and new end-cap vertices inherit
supports. The refined topology is fitted again for force densities; subsequent
Play and Snap operate on that same mesh. Its refined geometry becomes the new
reference and reset point, with a fresh iteration/history sequence. Repeated
smoothing compounds refinement. Graph or construction-parameter changes rebuild
the working mesh from the graph.

Only one generated mesh is published, in **Center Mesh**, with a stable scene
object ID. Low-poly and subdivision stages remain internal API snapshots.
**Fix XY** defaults on; uncheck it to permit free vertices to move in XYZ.
Snap then solves all three equilibrium components using the same sparse
factorization; supports remain fixed. Geometry attraction displays two decimal
places and penalizes displacement from the current reference in the free axes.

**Show / hide thrust and compression** styles this mesh's edges: compression is
RGB (0,40,240), matching Vault, tension RGB (240,0,140), and support-to-support edges black.
Active edge width uses the Vault linear 1 to 5 mapping of absolute force-density
magnitude (constant density maps to 1), multiplied by **Edge weight multiplier**.
The multiplier affects display only. Green support reaction vectors and optional
gold residual vectors occupy the **Forces** layer. **Thrust vector scale** controls
arrow length. Thrust is the force applied to the support (the negative of the
solver support reaction), matching Best Fit FDM. Explicit density colors and
widths survive selection/component-mode styling. The HUD reports minimum and
maximum active density magnitude. Solver, smoothing and graph changes refresh
displayed forces.

The sketch publishes `graph`, `Center Mesh`, `Stair Assembly`, `Supports`, and `Forces` layers.
File path changes
import automatically; **Reload Rhino file**
rereads a file saved at the same path. Editing the source externally in Rhino
requires saving and reloading; live updates concern the editable viewer graph.

## Stair layout and clearance

The viewer enables **Generate stair assembly** by default. Native callers retain
ribbon-only behavior unless `generateSteps=true`. Inputs are tread depth, target
riser height, landing depth, structural thickness, and
step clearance, all in the model's existing units.

Flights follow maximal graph paths between endpoints, junctions, or nodes inside
Rhino landing outlines. Landing elevations remain fixed. The riser count is
`ceil(totalRise / targetRiser)`, and every riser in that flight uses
`actualRiser = totalRise / count`. There are `count-1` tread goings; the upper
landing supplies the final walking surface. Tread depth remains exact. Flat
approaches share unused run. Insufficient run reports the flight number,
available run, required run and riser count, preserving the previous result.
The editable walking graph retains the original node and edge topology, with
one straight edge per input connection. Stair stations stay internal; they do
not add editing handles. Intermediate input
node elevations are adapted to that profile.

Closed horizontal curves in Rhino's Landing/Landings layers are used as landing
footprints, including parent layers of that name. Nodes are associated by both
XY containment and elevation, so stacked levels remain separate. Other endpoints
and junctions receive generated convex landings. Imported outlines remain fixed;
solved nodes must remain inside their assigned outlines. Treads and landings inherit the centre mesh plan, face boundaries and subdivisions without subdividing the base structural mesh by step stations.
Step contours split the stair surface faces into horizontal treads and landings without vertical riser surfaces. Shared source
vertices and edge intersections on each horizontal level are welded into a zero-thickness
**Stair Assembly** mesh above the structural mesh. No riser, bottom or side
closure is generated. Imported outlines constrain landing reserves but do not
add separate slabs outside the centre mesh. The legacy `stepThickness` input is
accepted for compatibility and has no effect on this surface. The structural mesh remains the same **Center Mesh** object.

Structural face ownership is propagated through construction and Catmull-Clark
refinement. Each face stays below its own flight's step undersides, accounting for
clearance plus half the structural thickness as a vertical thickness allowance.
For other flights it may be above or below their surface/riser envelope, without
intersection. Projected triangle/polygon overlap checks include contained
vertices and crossing edges; vertex-only tests are insufficient. Initial mesh
creation and post-FDM smoothing lower free structural vertices to satisfy the
own-flight envelope, then check the other flights. Fixed supports never move.

The Vault candidate validator rejects unsafe Play trials, allowing its existing
line search to shorten them. Snap backtracks between current and proposed
geometry if needed, reporting `direct-constrained` for shortened moves or
`clearance-blocked` if no acceptable move is found. This is a feasible direct
move, not a constrained quadratic-program optimum. Residuals still determine
convergence; clearance does not imply equilibrium.

After each accepted solve, steps and the simple walking graph regenerate from the
solved centreline plan, keeping flight counts, landing levels, tread depth and
uniform risers. The same editable graph object receives the world-space result
through its retained transform. The sketch retains the design graph separately;
solver writeback does not trigger rebuilding. User graph edits establish a new
design reference. Reset and rewind restore the coupled layout for their mesh
state. Post-smoothing establishes the existing new reset reference.

Verification includes the local Stairs_Graph.3dm with all eight imported landing
outlines. Range optimisation at width 1.00, tread domain 0.25-0.35 and riser domain
0.15-0.18 keeps the bottom elevation fixed, updates necessary landing levels and
generates a layout with the requested 0.05 clearance. This checks layout and
clearance feasibility; force residuals remain the equilibrium diagnostic.

Loads continue to use the vertex-mass input; they are not derived automatically
from tread or landing volumes. Assembly intersections with itself, headroom,
railings and construction detailing are not solved by this layout stage.

## Optimise graph levels

**Optimise graph levels** works before mesh generation and on an existing live
layout. Tread minimum/maximum and riser minimum/maximum define allowed domains.
**Use tread / riser ranges** defaults on; switch it off to use the exact Tread
depth and Target riser height as single-value domains. **Minimum landing depth**
reserves at least half that depth at each end of a flight, or more where its
landing footprint extends further into the flight.

For each connected component, the lowest input node is anchored at its existing
elevation. A degree-two bottom node receives a landing footprint if needed.
Flights are traversed outwards from that anchor. The optimiser first searches
integer riser counts and dimensions inside the domains that retain a landing's
existing elevation. It prefers values close to the target tread/riser inputs.
Where none fits, it selects the closest feasible elevation; subsequent flights
are then fitted from the adjusted parent level. New landings cannot fall below
the component's anchored bottom. Originally flat links retain a shared level.
A loop with incompatible assigned landing levels is rejected atomically; this
is a bottom-outward feasibility pass, not a global mixed-integer optimum.

Selected riser counts, tread depths and landing-node identities are retained per
flight. Imported Rhino landing outlines move vertically with their associated
nodes; XY footprints and the Rhino file on disk are unchanged. Intermediate
graph levels adapt to the selected flight profiles. A graph-only optimisation
leaves mesh generation gated; use **Generate base mesh from graph** afterwards.
When a mesh already exists, the candidate complete layout is validated before
committing graph/reference changes, and fitting resets for the new reference.
The chosen plan persists through live division, mass and solver changes. Graph
edits or changes to dimension/domain controls invalidate it; run the button again
to refit the dimensions. The viewer reports whether levels were retained or how
many landing levels changed, and shows selected tread/riser diagnostics.

The typed `zspace_stairs_optimise_graph_from_buffers` and
`zspace_stairs_optimised_graph_buffers` entrypoints share the existing ownership
contract. Legacy equivalents accept a filesystem input and return JSON. Input
adds `ranges:{treadMin,treadMax,riserMin,riserMax}`. Output includes `graph`,
updated `landingPolygons`, `landingOffsets`, `bottomNodes`, `landingVertices`,
`adjustedLandings`, `flightRiserCounts`, `flightTreadDepths` and
`flightActualRisers`. The chosen counts, treads and landing vertices can be passed
back through ordinary Stairs settings to generate exactly that layout. Invalid
optimisation leaves the current fitted mesh untouched.

## WASM contract

Opaque Stairs sessions use `zspace_stairs_create`, `zspace_stairs_update`,
`zspace_stairs_step`, `zspace_stairs_reset`, `zspace_stairs_previous`,
`zspace_stairs_result`, and shared destroy/error exports.
Legacy update takes `{graph:{vertexPositions,edgeConnects}, settings:{...}}`
from a filesystem path. Additive `zspace_stairs_update_from_buffers` and
`zspace_stairs_result_buffers` use the existing version-1 typed transfer contract.
The viewer uses typed transfers; result arrays are copied before further calls.
Other toolset sessions cannot be used through Stairs entrypoints or vice versa.
`zspace_stairs_snap(handle, geometryWeight)` performs the weighted direct solve;
`snapFixedDensities(weight)` is available in the native Vault and Stairs APIs.
`zspace_stairs_smooth_after(handle, divisions)` refines the working solver mesh.
Stair-enabled packets additionally include `stairMesh`, `stairFaceTypes` (0 tread,
1 riser, 2 landing), `walkingGraph`, `flightRiserCounts`, `flightActualRisers`, and
`minimumClearance`. Input `landingPolygons` is a list of flattened XYZ outlines.
Result packets include solver `edgeConnects`,
`activeEdgeIds`, `axialForces`, `residualForces` and `supportReactions` for display.

## Limits and verification

At most 4096 graph vertices and 4096 edges. Both subdivision counts accept
0..5 levels. The previous 1,024 active-edge / 10,000 vertex refinement
restrictions are removed. Larger systems use sparse iterative least squares
instead of allocating a dense equilibrium matrix. The 2 face / 4 smooth
division regression generates 6,144 faces and runs an accepted FDM step.
The sparse fit reports `rank = -1` because it does not measure numerical rank;
residual and convergence diagnostics still refer to physical equilibrium.
Subdivision remains subject to representable indices and available memory;
the 100 ms timer is a scheduling interval, not a guarantee of responsiveness
for all meshes. Geometry and fitting currently execute on the UI thread.
Self-intersections between remote arms are not resolved. This produces the
structural centre surface plus tread, riser and landing assembly geometry.
Railings and construction-code checks remain outside the current scope.

Checks: native Stairs smoke tests, legacy WASM tests, typed/legacy parity,
actual Rhino line import, automatic source/parameter update, deletion/disposal,
and the viewer regression suite. Build using `scripts/build_toolsets.bat` and
`wasm/scripts/build-wasm.bat`; run `node tests/stairs_wasm.mjs`, then the viewer's
`npm run build` and `node --test tests/parametricStairs.test.mjs`.
