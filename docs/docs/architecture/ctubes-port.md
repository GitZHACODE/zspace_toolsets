# C-Tubes initial port

`zTsCTubes` lives in the fabrication domain of zspace_toolsets. Reference:
EPFL-LGG/C-tubes, commit `fa04ca1`; MIT attribution is retained in
`docs/third-party/C-tubes-LICENSE.txt`.

## Available operations

- `setDirectrix(zObjectGraph&)`: order an open chain or closed loop by connectivity.
- `setGeneratrix(zObjectGraph&)`: preserve a planar section graph and express it in
  a local 2D frame. By default its vertex centroid is the profile origin; the first traversed
  edge defines local X. In single-profile mode the centroid is placed at the path start.
  `preserveOrientation` retains the supplied world plane and rotation; the viewer enables
  this by default. False retains the original automatic frame alignment.
- `setSettings`: spline samples, initial twist (radians), profile scale, and a
  constant apex scale for generic construction.
- `computeTube`: generic projection or parallel transport; returns a quad mesh.
- `computeStrips`: one longitudinal mesh per generatrix edge.
- `unrollStrips`: call core `zFnMesh::unroll` on each strip, then lay out flat meshes
  with a configurable gap. Geometry remains owned by `zObjectMesh`.

Spline sampling uses uniform control-point parameters, natural open boundaries
and periodic closed boundaries. Polyline mode retains the input graph stations.
Closed curves retain a duplicate end station: end-profile closure is not imposed
or optimized. Tube ends are uncapped. The initial profile frame uses a fallback
axis when the first tangent is vertical. The generic construction uses bisecting
section planes. Projection singularities and reversing tangents report errors.
Spline mode supports at most 512 control points and 4096 output stations; the
output is limited to two million vertices.

Core unrolling preserves planar faces and unfolds across their shared edges.
Nonplanar faces are triangulated before unfolding. Incompatible cycles,
non-manifold edges, disconnected components, and degenerate faces are rejected.
This is isometric strip unfolding, not general UV parameterization or nesting.
Self-overlap of flat strips and tube self-intersection are not resolved.

Not yet ported: the full Python optimisation framework/autodiff, arbitrary objective
combinations, directrix optimisation, dihedral parameterization, symmetry, connectors/networks,
end caps, fabrication flaps, and strip splitting. The implementation has native
analytic geometry/metric tests; full numerical parity with upstream Python has
not yet been established.

## Toolsets WASM

`wasm/scripts/build-wasm.bat` compiles core, interface, IO, and modern toolset
sources into `wasm/out/zspace_toolsets.js` and `.wasm`. Core provides its existing
Emscripten shims. No desktop display or host SDK is linked. Generated output,
cache, and temporary files are ignored by git. The core checkout must include
the new `zFnMesh::unroll` declaration and implementation.

ABI version 1 has opaque integer session handles, an error accessor, input JSON
via the Emscripten filesystem, and bulk position/face arrays. The first session
type is C-Tubes. Copy arrays before the next mutation or memory growth. Destroy
handles when finished. Exceptions are caught at the exported operation boundary.
The module is independent from `zspace_core.wasm`; pointers are never shared.

Build and test:

```powershell
scripts/build_toolsets.bat
ctest --test-dir build/ninja-msvc -C Release --output-on-failure
wasm/scripts/build-wasm.bat
node tests/ctubes_wasm.mjs
```

Copy both runtime files from `wasm/out` into the sibling webviewer's ignored
`wasm` directory, then run its `npm run build`.

## Viewer

Select the C-Tubes sketch and show Sketch Controls (S) and Outliner (O).
Optional example graphs are supplied through an action button. Select the
directrix graph and assign it, then select/assign the generatrix. Choose Generic
or Parallel, press Compute C-Tube, then Compute strips and unroll.

Assignments are snapshots in world coordinates. Reassign after editing a source
graph. Parameter edits take effect on the next Compute; unrolling operates on
the last computed tube. Generated groups update in place. The current viewer
stores graph render buffers, so coincident graph vertex positions are rejected
when recovering connectivity. The sketch context exposes selected-graph snapshots,
named output groups, and the toolsets adapter; it never performs tube math in JS.

## Multiple target sections and optimisation

`setGeneratrices(zObjectGraphArray&)` supplies target cross-sections for ONE tube.
With two or more targets, `computeTube(Generic)` runs a parameter-based fit using
per-segment apex scales and, optionally, intermediate section-plane normals. It
always regenerates vertices using the generic C-Tube projection recurrence. It
does not morph or loft mesh vertices independently.

Target world positions and planes are retained. Target-plane intersections are
inserted into the sampled directrix; the closest intersection to each profile's
centroid is chosen if there are several. Targets are sorted along the path, with
the earliest section held exactly. Construction propagates both forwards and
backwards from that section. Remaining sections are least-squares targets, not
guaranteed hard constraints. Cyclic seam shifts and reversed winding are matched at each objective evaluation.

The fit uses analytic chain-rule gradients and LBFGS++ L-BFGS-B with
More-Thuente line search. Vendored MIT source is pinned in
`third_party/LBFGSpp/UPSTREAM.txt`. Settings follow the reference notebook's
absolute projected-gradient tolerance 1e-5, relative objective tolerance 1e-10,
and memory size 10; default budget is 200 iterations, maximum 2000.

Apex controls remain initial apexScale * exp(q), bounded q in [-3,3]. Free
normal tangent-plane coordinates remain bounded to [-3,3]; the seed normal,
directrix and seed profile are fixed. These are existing model restrictions,
not the paper's complete optional Q/theta degrees of freedom. LBFGS++ is a
C++ implementation of the same solver family, not SciPy's identical binary.

Production optimisation uses no finite differences. The optional native
checkFitGradient diagnostic compares every analytic derivative with central
differences at a perturbed state. Tests cover apex, normals, backwards seed
propagation, multiple targets, weights and smoothness. Singular trials return
infinite objective so line search can shorten its step; solver failure retains
the best valid evaluated state and reports the stopping reason.

`fitReport()` / WASM `zspace_ctubes_diagnostics` expose initial/final RMS, maximum
point error, iterations, convergence, target stations, and all fitted apex scales.
Converged means maximum paired point error <= fitTolerance in input coordinate
units. An incompatible target set may return a best fit above tolerance; the viewer
reports this explicitly. Directrix geometry is fixed. Single-profile scale/twist
controls do not change multi-profile world-space targets.

Initial fitting limits: open path, Generic mode, 2-16 planar targets of matching open/closed topology with
matching vertex counts, <=128 original path stations and <=128 vertices/profile.
Coincident target stations and non-intersecting target planes are rejected.
Parallel transport remains a single-profile alternative, with no apex optimisation.

## Frames, layers and strip metrics

C++ returns each actual section frame (origin, X, Y, normal). The viewer displays
selectable RGB triads in Frames; select one to show its station in the HUD and
coordinates in the channel box. Update frame HUD applies size/stride/visibility
controls without recomputing the geometry.

C Tube, Strips, Unrolled Strips, Frames and Sampled Inputs have separate layers.
Each unrolled strip identifies its 3D partner and displays source area, flat area,
signed difference, and absolute percentage deviation: 100*abs(Aflat-A3D)/A3D.
Areas are measured in C++ on triangulated copies of the computed meshes before
viewer transforms. Tiny floating-point differences are expected. Equal area does
not prove that a pattern is overlap-free or fabrication-ready.

## Weighted objectives and picked geometric origins

Global settings: `targetSectionWeight`, `targetPlaneWeight`,
`apexSmoothnessWeight`, `planeSmoothnessWeight`. Optional `sectionWeights` and
`planeWeights` arrays multiply the global weights per input generatrix (input
order, even if stations are sorted). Nonnegative finite weights up to 1e8 are
accepted; zero disables a term. Native/WASM smoothness defaults are zero for
compatibility; the sketch enables both at 0.01. Section and plane weights default
 to one. The seed section remains exact regardless of its weights.

Section energy is the p=10 norm (reference smooth_max) of per-target mean
paired squared 3D distances, multiplied by their section weights. There is no
seed-radius normalization; existing weight settings change numerical balance. Plane energy is mean squared cross product of fitted and target unit
normals (unoriented planes). Each plane passes through its fixed directrix
station; this objective adjusts orientation, not station position. Turning off
plane optimisation fixes all normals and makes plane weights ineffective.

Smoothness terms are half the mean squared second differences of actual apex
scales and unit plane normals, matching the reference discrete energy form.
They use station index, not arc-length spacing. Fewer than three values give
zero second-difference energy. Analytic-gradient L-BFGS-B minimises the combined energy with the existing
weak parameter regularisation (1e-10 squared control norm). It no longer
stops early merely because target fitting is within tolerance. Reports include
each weighted energy, total initial/final objective, per-input RMS and plane
angle errors. `converged` still means ALL target points are within fitTolerance,
including targets with zero weight; smoothness may deliberately trade this off.

`setGeneratrixOrigin(profileIndex, vertex)` picks an existing world-space vertex.
WASM accepts `originVertex` on each generatrix graph, using its original buffer
index (not traversal order). A single section is translated so the selected
vertex, instead of its centroid, is at the path start; scale/twist are applied
about that anchor. A target section is translated in its supplied plane so its
picked vertex lies at the plane/directrix intersection. Its shape and rotation
are retained. The anchor does not change the seam or vertex order. Non-seed
sections remain soft fitting targets. No picked origin preserves previous
single-profile centroid / multi-profile world-position behaviour.

Viewer workflow: assign a curve/graph; its sampled zObjectGraph is published in
Sampled Inputs. In vertex selection mode pick exactly one vertex on that graph,
then use `Use picked graph vertex as geometric origin`. The graph snapshot
includes viewer world transforms. Recompute after changing origins or weights.
Choose a 1-based target number and apply its section/plane multipliers using the
separate action; reset its origin to restore the default. Channel diagnostics
show objective terms and per-target fit errors. The original Rhino curve is kept.

Solver reports include `gradientNorm` (projected infinity norm), `evaluations`,
and `termination`, alongside geometric `converged`. Stationarity is distinct from
matching the targets. Directrix/seed constraints and incompatible section shapes
can still prevent an exact fit. Model units now affect the balance with smoothness.

Validation on the supplied equal-count (8/8) profiles, with both smoothing weights
0.01: previous solver at 100 iterations RMS 0.0951101, max 0.222430; L-BFGS-B
stopped at 101 iterations on objective-change tolerance, RMS 0.0951097, max
0.222429. A 2000-iteration limit produced the same result. This validates the
solver replacement, not exact matching or equivalence with all paper variables.

## Open profiles and separate Rhino lines

Generatrices may be planar connected graphs, including open chains, closed
polygons and branched T sections. Connected opposite collinear edges combine
into one fabrication strip, including through a branched junction. A T section
therefore unrolls as two strips: its continuous crossbar and its stem. Straight
subdivisions stay in the merged strip mesh, preserving its vertices and face areas.
Continuity is checked at every sweep station with a tolerance of 0.001 times
the shorter incident edge length plus coordinate float-precision tolerance;
noncollinear edges remain separate strips. Open chain
matching permits forward/reversed order, with no cyclic shifts. Branched targets
must share vertex counts and isomorphic connectivity; the fitter tests allowed
graph correspondences independently of source vertex IDs. Straight lines
have no intrinsic section plane: the fitter uses the nearest directrix tangent
projected perpendicular to the line to choose that plane. Single-section framing
uses the initial path tangent. Planarity and valid projection checks still apply.

The browser joins touching curve endpoints within each source layer using the
Advanced line join tolerance (default 0.001 model units), without joining
interior-to-interior crossings or curves across layers.
An endpoint touching an edge interior splits that edge into a shared junction.
Each connected component is one section graph; free ends and branches are
preserved. Separate Directrix lines similarly
must form one connected open path. Source curve IDs and layer metadata are
retained. Differently sampled open chains can share normalized-length subdivisions
without adding a closing edge; original corners and picked anchors are retained.

Section planes use noncollinear vertices from the whole graph, independent of
edge order. The section need not physically touch the directrix. Multiple target
planes still locate stations by intersection with the directrix. A selected
origin vertex is the geometric anchor; otherwise single-section generation uses
the vertex centroid. Nonplanar sections and branched directrices are rejected.
The internal sweep shares one position per section vertex. Because a T sweep
has three incident sheets along its junction, the manifold zObjectMesh output
uses coincident per-edge seam vertices there. Fabrication strips combine straight
edge chains and each resulting strip unrolls independently.

Target intersections replace nearby display stations within 1e-5 of the path
extent plus coordinate float-precision tolerance. The exact target intersection
is retained. This avoids near-endpoint sliver faces that can fall below the
unroller's global metric tolerance after optimisation. Endpoint replacement can
trim a correspondingly tiny portion of the input path; other stations remain.
## Branching boundary sketch and prescribed frames

The viewer's separate `Branching C-Tubes` sketch loads an editable absolute OBJ
path through a same-origin, uncached local `.obj` endpoint. It currently supports
Y quad ribbons with one interior degree-three junction and three end caps.
Removing the cap midpoint vertices from the boundary loop yields three open
boundary directrices. Interior longitudinal edges plus cap midpoints define the
central green graph. Invalid or ambiguous topology is rejected explicitly.

`Frame spacing` samples arc length on each boundary polyline, independently of
mesh transverse edges. The exact junction station is inserted in addition to
regular spacing stations, so that F1/F2/F3 share the central junction. Central
correspondence is interpolated along the original green graph edges; every
frame X axis runs from the boundary origin through that central point. Frame Y
comes from interpolated mesh surface normals, orthogonalized against X. Thus
section planes contain the central graph targets. Spacing and T lengths use
unscaled OBJ coordinate units.

Loading or reloading shows only the input mesh. `Generate C-Tubes` constructs
the initial three tubes with the supplied boundary directrices and preserves
their inward flange gaps. Optimisation is a separate action.
`Extend flanges to central graph` instead extends every initial section to its
central target while retaining boundary directrices and section planes. It is
a geometric preview: closure alone does not ensure developability.
The initial T flange/stem junction is on the boundary directrix; the stem
extends opposite the surface-normal axis.
`zTsCTubes::computeTubeWithFrames` performs generic directional projection
through these fixed section planes. The WASM input's optional `frames` array
selects this path; the existing unconstrained/optimised path stays unchanged.
The outward flange and stem dimensions describe the seed section; later projected
dimensions can vary. Optional WASM `flangeTargets` and `flangeVertex` constrain
the inward endpoint. Both flange edges stay collinear through their common anchor.
Adjacent boundary flanges share an identical central-edge sampling, including
central polyline corners and the union of the distance stations on both sides.
These extra seam stations prevent gaps between mismatched chords; F1/F2/F3
share the central junction. The inward flange meshes have coincident matching
seams but remain separate objects for fabrication.

`zTsCTubeBranching` uses generic conical directional projection with the green
central seam as a hard construction directrix. It optimises log apex scales and
section normal directions with L-BFGS-B, using analytic chain-rule gradients
checked against central differences. The generated flange/stem anchor fits the
original boundary: boundary directrices may move, while central seam endpoints
remain fixed. Plane priors and apex/normal smoothness regularise the fit. The
inward flange width varies to reach the shared seam. This is a constructive
extension of the paper's generic projection, rather than its full Q/theta model.

The browser runs optimisation in a WASM worker, then verifies exact seam
incidence and strict `zFnMesh::unroll` acceptance before publishing. Failed
results preserve the initial tubes. Diagnostics show boundary RMS/max movement,
iterations and gradient-check error. Advanced controls expose iteration count.
Collinear flange edges combine into one fabrication face column when coplanar
at the kernel's existing precision, avoiding redundant unfolding cycles through
the T anchor. Three tubes each provide a flange and stem strip; all six strips
show area comparisons. Closure and unroll checks do not certify absence of
geometric overlaps at the junction or global self-intersections.

`Optimise central graph (fixed directrices)` selects the coupled alternative.
Each boundary remains a hard construction directrix. Shared central coordinates
are single optimisation unknowns across the three tubes, so their tips coincide
exactly. Panel coplanarity uses scalar triple-product constraints, analytic
gradients, an augmented Lagrangian and L-BFGS-B inner solves. Central movement
and displacement smoothness regularise the result. The next stem ruling is
the intersection of its section plane and the plane through the previous stem
and fixed boundary edge. Its length stays at the input stem length; the angle
within the section may vary to maintain planarity. The outward flange remains on
the opposite side of its boundary anchor at the specified length; its panels
share the optimised inward flange's plane. Analytic Gauss-Newton projection
initialises and restores panel planarity around the L-BFGS-B solves.
Publishing requires panel/construction residual at most 1e-5 OBJ units and all
six strict unroll checks to pass. Diagnostics report central movement, seam gap
and construction residual. This mode extends the paper's planar-panel conditions
with central-vertex objectives and stem ruling transport; it does not use a
single common conical apex for the whole T section. The original boundary-
optimisation mode continues to use generic conical C-Tube projection.
The iteration parameter limits each inner solve; reports count total iterations.

`Fix XY (move only Z)` locks the imported world X/Y coordinates of the moving
graph using equal L-BFGS-B bounds. The Gauss-Newton projection also removes
its X/Y Jacobian columns, so it cannot bypass the lock. In central mode the
boundary rail remains fixed; in boundary mode the central rail remains fixed.
With this constraint enabled, both modes use the symmetric planar-panel solver
and developable stem rulings, rather than the boundary mode's conical apex
parameterisation. Section axes and stem directions may still rotate; Fix XY
applies to the graph vertices. Z-only constraints can be infeasible: failed
planarity or strict unfolding preserves the previous displayed geometry.
`tests/ctubes_branch_fixxy_wasm.mjs` verifies exact X/Y preservation in both
modes, small displacement for a local two-panel defect, and six strict unrolls.

The branching sketch also has `Generate C-section tubes`. The corrected profile
is an open chain: upper central endpoint, boundary anchor, lower outer corner,
lower central endpoint. It has a top flange, outer web and bottom flange; it
opens toward the green central graph. Adjacent Cs share both open endpoints,
so their cross sections form a closed loop without an internal central web.
Depth d is the shared world-Z offset from the upper to the lower central graph.
Theta is the outer web tilt from downward world Z, in degrees; zero is vertical
and positive theta moves its lower corner toward the centre. Section planes
contain the top flange and world Z. Both central graphs are displayed.
Initial generic sweeps use the existing C++ fixed-frame construction and may
have gaps away from the seed. Extension prescribes every upper/lower central
endpoint exactly, with the requested d and theta. Optimisation uses the existing
rail solver and builds the same paired C geometry on the solved rails, retaining
d and theta. Strict C++ prescribed-section and zFnMesh unroll checks still reject
invalid fabrication results before publication. Three strips are produced per
branch. No C++/WASM API changes are required for this profile. Viewer tests in
`tests/ctubeCSection.test.mjs` verify web angle, topology, and exact upper/lower
endpoint coincidence between two Cs at multiple stations.
