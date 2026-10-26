# Vault: Best Fit Force Density Method

`zTsVault` is the statics toolset entry point. This first version implements
mesh-edge Best Fit FDM, an optional fixed-XY constraint, and reconstruction of
TNA form/force diagrams. Regular FDM, interactive TNA editing, graph input, and
additional geometric constraints remain future extensions.

## Reference and interpretation

Lorenz Lachauer and Philippe Block, **Interactive Equilibrium Modelling**,
International Journal of Space Structures 29(1), 2014, pp. 25-37.
The implementation is independently written from sections 3.1-3.4. No Rhino or
NumPy code is bundled. The archived pre-v3 Vault implementation was consulted
for context; its hard-coded boundary constraints and bounds are not retained.

For branch incidence C and coordinates X, A stacks C_free^T diag(C X_c).
Positive q is tension; negative q is compression. External loads P are physical
vectors. The algebraic residual is A q - P; displayed free-node force is its
negative, P - A q. A support reaction is (A_all q - P) at a support, including
loads applied directly to that support. Axial force is q times branch length.

The initial fit solves min ||Aq-P|| with Eigen complete orthogonal decomposition,
then clamps q into its interval, matching the paper's two-stage initialization.
This is not a bounded least-squares optimum and is not a minimum-distance
geometry optimization. Iteration splits residuals using geometryWeight and
uses the equation 14 ratio sum(length*q)/sum(length*abs(q)); a zero denominator
gives zero geometric correction. Corrections have negative signs for the
documented algebraic residual, verified by tension and compression fixtures.

Numerical departures/clarifications: A is 3n by m (the printed 3n by 3m is a
typo); geometryWeight is in [0,1]; equation 14's ratio is used rather than the
prose's +/-1 claim. Bounds are projected at each integration stage. Adaptive
explicit Heun RK2 uses an Euler error estimate, bounded retries, and rejects
steps that do not decrease the residual when edge regularization is disabled.
Regularized runs instead use the combined merit described below. This is not a reproduction of the
paper prototype's unspecified adaptive Runge-Kutta implementation. Rank is
reported; normal equations and matrix inversion are avoided.

Convergence is ||r|| <= max(absoluteTolerance, relativeTolerance*initialResidual).
The initial residual is measured after clamping. Relative tolerance defaults to
0.01. A stalled or iteration-limit result does not claim equilibrium.

## C++ use

```cpp
#include <zspace/toolsets.h>
zSpace::zTsVault vault;
vault.setInputMesh(mesh); // Copies through zFnMesh, preserves vertex numbering.
vault.setSupports({0, 1, 2, 3});
vault.setTensionEdges(tensionEdgeIds); // zFnMesh::getEdgeData order.
vault.setVertexMass(1.0); // Or setVertexMasses: one nonnegative value per vertex.
vault.setGravity(zSpace::zVector(0, 0, -9.81));
// Optional setVertexLoads adds explicit force vectors to mass * gravity.
zSpace::zVaultBestFitSettings settings;
settings.qMin = 0.01;
settings.qMax = 1000;
vault.setBestFitSettings(settings);
auto report = vault.computeBestFitForceDensities(); // Geometry unchanged.
report = vault.stepBestFit(1); // Call repeatedly or solveBestFit(iterationBudget).
auto& result = vault.resultMesh(); // Treat as read-only while solving.
auto& q = vault.forceDensities();
vault.previousBestFit(); // Up to 256 accepted steps, including q and timestep.
vault.resetBestFit(); // Restores initial fitted state.
```

qMin/qMax are **magnitudes**, with 0 <= qMin <= qMax and qMax > 0.
Tension edges use [qMin,qMax], remaining edges use [-qMax,-qMin]. Setting
compressionOnly=false allows remaining edges in [-qMax,qMax]. Optional signed
per-edge `setForceDensityBounds(lower,upper)` intersects those intervals.
Conflicting intervals fail. Empty arrays clear the per-edge override.

Mass is explicit nodal mass, not computed tributary area or volumetric density.
P_i = explicitLoad_i + mass_i * gravity. Use consistent length/mass/force units.
The viewer defaults to gravity magnitude 1, so mass 1 generates a unit load;
use 9.81 for SI kg/m/s units. Support vertices may carry mass. Support-support
members are unidentifiable in free-node equilibrium, excluded from fitting,
and return q=0 with no active-edge flag, irrespective of their sign assignment.

Input changes invalidate the solver and require a new fit. `setInputMesh` resets
assignments and settings; set it first. Geometry is owned by zObjectMesh and
accessed using zFnMesh. Eigen state is private to the toolset. Existing core
operations suffice; this change does not modify zspace_core.

## Shared WASM

The existing zspace_toolsets.js/.wasm module includes Vault and C-Tubes.
ABI version 1 is extended additively. Opaque handles distinguish session type;
wrong-type operations fail. `zspace_vault_fit` reads JSON through the module FS
and replaces session state only after a successful fit. Input contains mesh
polygon arrays, supports, tension endpoint pairs, masses, gravity, settings,
and optionally loads and qLower/qUpper arrays. Endpoint pairs are resolved in
C++ against original mesh edges. qLower/qUpper use the returned kernel edge order.

`zspace_vault_step`, `reset`, `previous`, `select_mesh`, and `result` expose the
session. Geometry uses existing bulk mesh buffers; diagnostics use JSON with
double precision force densities. The adapter immediately copies buffer data.
Core and toolsets heaps never share pointers. Destroy sessions on sketch disposal.

## Build and verify

From a VS 2022 x64 Developer Command Prompt at the toolsets root:

```bat
cmake -S . -B build/native -G Ninja -DZSPACE_TOOLSETS_BUILD_TESTS=ON
cmake --build build/native --parallel 4
ctest --test-dir build/native --output-on-failure
powershell -NoProfile -ExecutionPolicy Bypass -File wasm/scripts/build-wasm.ps1
node tests/vault_wasm.mjs
node tests/ctubes_wasm.mjs
```

The WASM script accepts `-CoreDir` and `-EmsdkDir`; defaults are the sibling core
and `%USERPROFILE%/emsdk`. It discovers new toolset cpp files automatically.
Copy **both** generated runtime files from `wasm/out` to the sibling viewer's
ignored `wasm` directory. In the viewer run:

```powershell
node --test meshStatics.test.mjs ctubes.test.mjs solverPresentation.test.mjs
npm run build
npm run dev
```

## Viewer workflow

Select a mesh, then choose **Mesh Statics → Best Fit Force Density**. Mesh Statics
is located after Mesh Solver. The mesh gains a **Best Fit FDM** node under the
channel box's **Mesh Solver** history. No sketch controls open.

Use the Mesh Statics menu to set selected support vertices, boundary supports,
and selected tension edges, or clear either assignment. Edit mass, gravity,
qMin/qMax and numerical/display settings in the channel box. The mass field
assigns selected vertices when vertex components are selected; otherwise it
assigns all vertices. Other vertices retain their assigned masses.

**Fit Force Densities** fits the current mesh; **Resolve One Step** or the existing
bottom Step/Play controls updates that same mesh in place. Changed geometry,
constraints or numerical settings trigger a fresh fit on the next solve. Display
settings update overlays immediately. Previous/Reverse restores recorded states;
Reset returns to the last fitted input; frame seeking resets and advances.
Assignments and parameters participate in the existing object undo snapshots.

Black points indicate supports, magenta edges tension (RGB 240, 0, 140), blue edges
compression (RGB 0, 40, 240), and gold arrows residual forces. Edges connecting two
supports are black regardless of their tension assignment. Tension and compression
colors are editable in the channel box; color edits preserve the fitted session.
Signed min/max force densities are reported in the channel box and top status
strip, excluding inactive support-to-support edges. After fitting, active edge
widths map |q| linearly to 1–5 pixels, multiplied by Edge Factor. A constant range
uses 1 pixel before scaling. Before fitting and for inactive edges, existing
hard-edge thickness is retained and scaled. Edge Factor zero hides statics edges.
This follows the 1–5 display mapping in the C# `zTsVault.FDM.cs` implementation
called by `projects/zspace_components/zspace_rhino_components/src/zStatics_FDM.cs`,
but uses force density magnitude rather than its axial-force magnitude |q L|.
Display support thrust draws vectors in RGB 0, 240, 140, using the negative of
the reported support reaction (load exerted on the support), scaled by Force
display scale. Edge Factor and thrust visibility are display-only controls.
The shared top solver status strip reports assignments, iteration, residual
norm and maximum residual. Calculations use mesh-local coordinates and local −Z
gravity, preserving the object's transform. The standalone action sketch remains
as an unregistered API example; the normal UI uses native solver history.

## Limitations

- Small fits retain the dense rank-revealing solve. Above 1024 active edges or
  8 million matrix entries, initialization and gradient solves use sparse
  iterative least squares (up to 1000 iterations, relative tolerance 1e-8).
  These are dispatch thresholds, not input limits. The 10000-vertex cap is
  removed. Large rescue steps use sparse conjugate gradients on the damped
  normal system. Memory availability and representable indices still apply.
- Sparse iterative fits report rank=-1 (not measured). Their finite best iterate
  can be used even when the inner iteration budget is exhausted; outer residual
  and convergence diagnostics do not claim equilibrium unless their tolerance
  is met. No dense fallback is attempted on these larger systems.
- The paper assumes mechanisms m < 3n. Other counts are accepted using a
  least-squares fit, explicitly flagged in the viewer, without convergence claims.
- Fixed topology, axial members only. No material stiffness, bending, buckling,
  collision detection, stability assessment or structural design certification.
- Nonzero free-node loads required. No prestress-only nullspace selection.
- No line/plane constraints or user-fixed nonsupport nodes yet.
- Mixed-sign force cancellation can stall. Bounds may be infeasible; nonconvergence
  is reported, not classified as a proof of infeasibility.
- Numerical state is double; zSpace geometry and renderer buffers are float.
- History holds 256 accepted steps. Browser computation is on the main thread;
  one numerical step per animation frame does not guarantee a frame rate.
- Solver sessions are runtime-only; scene file persistence remains deferred.

## Verification record (2026-09-16)

Viewer integration revision: production build and 12 targeted regression checks
pass, including copied constraint history, session invalidation and force overlays.
The numerical kernel was unchanged for this revision. Chrome UI checks confirmed
menu placement, creation of a Best Fit FDM solver node, channel-box mass/q
controls, boundary support markers and HUD counts. The existing timeline reached
the toolsets solver and correctly rejected a fully supported mesh with no free
vertices. Shared WASM regression fixtures verify residual reduction and history.
Full interactive playback on a partially supported mesh was not reverified in
this UI revision.

The following browser measurements were made with the original action sketch:

Native Vault and C-Tubes tests pass. Shared-module Vault and C-Tubes WASM tests
pass. The selected 22 viewer regression tests pass; the production viewer builds.
In the browser, the example with qMax=0.2 reduced residual 0.68 to 0.0065609 in
100 accepted iterations. Step 10, Previous (10 to 9), Reset (back to residual
0.68), changed-parameter refit requirements and explicit mesh assignment were
verified. Assigning mass 2 changed the initial residual to 1.68 at the same bounds.
Native fixtures also cover mixed tension/compression and signed per-edge bounds.

Separate existing failure: `zspace_toolsets.3dp.smoke` returns failure without
diagnostic output; its code and the 3DP implementation were not changed. The
viewer retains its bundle-size warning. Initial shell builds also encountered
missing MSVC environment variables and sandbox process/temp-file permissions;
the developer-shell native build and permission-enabled WASM/viewer builds pass.
## Edge-length regularization

The Best Fit FDM channel box exposes **Min / max edge length ratio** (`0..1`)
and **Edge regularization strength** (nonnegative; default `0`, disabled).
For example, ratio `0.5` requests that the shortest active member be at least
half the longest. Set a positive strength, then Fit or Play. The current ratio
and normalized regularization error are shown separately from the force residual.
Lengths are full 3D lengths; support-to-support members are excluded.

For active lengths `Li`, define deficits `di = max(0, ratio * max(L) - Li)`.
The correction is the negative gradient of their mean squared deficit: short
edges are extended and longest edges shortened, sharing the latter contribution
when multiple edges tie for maximum length. Supports receive no displacement;
Fix XY projects the correction to Z. Geometry weight scales this correction,
so geometry weight zero still means q-only updates. These are geometric correction
forces, not additional physical loads or support reactions.

With regularization enabled, the accepted-step merit is
`||physicalResidual||² / ||loads||² + strength * geometryWeight * mean(di²) / Lref²`,
where `Lref` is the longest active edge at the initial fit. Physical residual may
temporarily rise as lengths improve. Convergence requires both the existing force
tolerance and normalized RMS length deficit at most `1e-4`. The ratio is a soft
objective: incompatible supports, q bounds, or Fix XY can lead to a compromise
or stalled solve. Strength zero preserves the previous residual-only behaviour.

C++ settings are `minMaxEdgeRatio` and `edgeRegularizationStrength`, passed with
the same names in shared WASM `settings`. Native/WASM tests cover an initially
balanced unequal mesh, ratio improvement, support preservation, Fix XY, and
parameter validation. Native tests also verify reset/rewind and disabled behaviour.

## Plan constraints and persistent reference

The Best Fit FDM channel box has **Plan constraint: Free / Penalize XY / Fix XY**
and **XY penalty weight** (nonnegative, default 1). Free is the default.
Penalize XY attracts free vertices towards the assigned mesh's XY coordinates;
Z is unpenalized. Higher weights resist horizontal displacement more strongly.
Fix XY retains the current fit input's XY exactly. Supports remain fixed in all
three coordinates in every mode. Existing `fixXY=true` assignments migrate to
the hard mode; the legacy native/WASM flag remains supported.

The viewer snapshots reference positions and topology when Best Fit FDM is
assigned. Refits, q-bound/mass/weight changes, playback, reset and rewind do not
rebase this snapshot. Legacy assignments without a reference capture one on
their next fit. **Mesh Statics → Set XY Reference from Current Mesh** explicitly
replaces the snapshot. Topology changes require that action before the next fit.
References are copied with solver undo/redo state, never shared by mutable alias.

Native callers use `zVaultBestFitSettings::planConstraint` with
`zVaultPlanConstraint::{Free,PenalizeXY,FixXY}` and `xyPenaltyWeight`.
`setInputMesh` initializes the reference; `setReferencePositions` supplies a
separate vertex-corresponding reference. Shared WASM uses numeric modes 0/1/2,
`settings.xyPenaltyWeight`, and optional flattened `referencePositions`.

The soft correction is proportional to `(reference.x-current.x,
reference.y-current.y,0)`, averaged over free vertices and scaled by geometry
weight and XY penalty weight. The accepted-step merit includes the weighted mean
squared XY displacement normalized by the squared longest reference edge.
This reference scale is independent of the moving mesh. Weight zero or geometry
weight zero disables soft attraction. RMS/max XY displacement are reported in
mesh units over free vertices, separately from physical force residuals.

This adapts the weighted target-vector idea in Bhooshan, Van Mele and Block,
*Discrete funicular structures with curve-folded moulds*, IASS 2015,
§3.4.1(4), pp. 9–10 (user-supplied `PAPBhooshan.pdf`). It integrates attraction
into the existing coupled geometry/q solver, rather than reproducing the paper's
complete q-update and FDM reconstruction loop. Attraction is a design correction,
not a physical load. A soft target can retain nonzero displacement at a compromise;
convergence checks physical residual and stationarity of combined geometry updates,
not zero XY displacement. Incompatible goals can produce a stalled solve.

The hard constraint removes X/Y geometry updates at both integration stages; it
does not remove horizontal equilibrium equations or residuals. Consequently,
an infeasible fixed plan may stall with a nonzero residual. Geometry weight
still applies: zero permits only q updates, and one permits only Z updates
after the initial q fit. Reset and rewind preserve the constraint.

## TNA diagram reconstruction (Stage 1)

After fitting, choose **Mesh Statics → Create / Update Form and Force Diagrams**.
The linked `_form` and `_force` objects inherit compression/tension colours
and Edge Factor from the source mesh. The source channel box controls visibility,
X offsets, force-diagram display scale, and automatic updates. Changed solver
inputs invalidate the diagrams; they are hidden until reconstructed.

`createFormFromResult()` projects the fitted result to XY.
`createForceFromForm(settings)` reconstructs both diagrams from the current
fit, and `updateForceFromDensities()` repeats the reconstruction. These methods
leave the fitted mesh and q unchanged. Geometry is available through `formMesh()`
and `forceMesh()`, with bidirectional halfedge maps and `diagramReport()`.
The shared WASM export `zspace_vault_diagrams(handle, tolerance)` returns copied
JSON geometry, signed densities, maps, and diagnostics.

The implementation follows the archived pre-v3 `zTsVault` boundary-region
construction, splitting exterior regions at supports. Eigen SparseQR reconstructs
dual vertex positions from signed parallel target vectors `q * (Xj-Xi)_XY`,
with one translation anchor per component. Force-edge lengths represent horizontal
force `|q| * Lxy`; display scaling does not change physical q. Support-support
members remain black in the form and are omitted from the force diagram.

Unsupported boundary loops, including internal openings, each form a single
exterior dual region. Loops containing supports are split at those supports;
supports are not required on every loop.

Stage 1 requires consistently oriented manifold faces,
at least three incident edges at free vertices, and no
collapsed projected active edges or horizontal loads on free vertices. Unsupported
or ambiguous dual mappings are rejected. Mixed signs are supported; force polygons
can overlap or self-intersect and are displayed as coloured edges.

The HUD reports horizontal residual norm, maximum target-polygon closure error,
and maximum reconstructed edge-vector mismatch. A diagram is labelled approximate
when closure or mismatch exceeds the absolute tolerance (default `1e-5`). This is
a horizontal reciprocity diagnostic, not a claim of full 3D equilibrium or stability.
Native and WASM tests exercise signed q round trips and closure; viewer tests cover
linked object reuse, colour propagation, display scale, and stale-output hiding.

### Stalled-step fallback

If the adaptive Heun update cannot accept a step, a sparse damped Gauss-Newton
fallback minimizes the existing normalized residual-plus-penalties merit.
Its analytic Jacobian includes free coordinates, force densities, XY attraction,
and active edge-length deficits (averaging derivatives at tied longest edges).
Coordinate and density columns are weighted by the square roots of geometry
weight and its complement. Supports and Fix XY are eliminated; fixed q intervals
are eliminated and other q updates are projected into their signed bounds.
Sparse LDLT solves damped normal equations, with increasing damping and backtracking.
Only finite, noncollapsed trials with strictly lower merit are accepted.
Fallback acceptance uses objective reduction rather than the Heun integration-error
tolerance, since this is an optimization step, not an integration step.
Accepted steps use the same history/reset path. Stalled now means both methods
failed. This is not a certificate of infeasibility or structural stability;
projected bounds, ill-conditioning, and nonsmooth penalties can still block progress.
