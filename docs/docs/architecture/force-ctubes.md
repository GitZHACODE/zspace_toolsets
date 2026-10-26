# Force C-Tubes

`zTsForceCTubes` composes `zTsForcePolyhedra` and `zTsCTubes` for frozen-force
section sizing of one open arch. The Alice sketch is **Force C-Tubes**, separate
from **Force Polyhedra**, and reuses its input, AP/3DGS and timeline implementation.

## Inputs and workflow

1. Import an open Directrix from Rhino, sampled by the existing nearest-spacing
   equal-arc-length method. Guide vertices must follow chain order (2–128 stations).
2. Create force polyhedra and advance AP/3DGS as needed. Unconverged snapshots may
   be inspected but are labelled provisional, independently of section sizing.
3. Supply a connected planar generatrix (3–128 vertices), either by selection or
   from the Rhino Generatrix layer. A generatrix number selects one connected
   profile; multiple profiles are not target constraints in this first version.
4. Choose a closed solid polygon, or thin-wall material with fixed thickness.
   Thin-wall properties integrate rectangular edge strips, including their local
   thickness inertia; overlaps at corners/junctions are neglected. Holes and
   exact thickened-profile Boolean unions are not supported.
5. **Generate C-Tube** constructs the baseline. **Optimise section areas** sizes it
   against the current force snapshot. Neither action advances/resets 3DGS.

The geometry guide remains fixed and separate from the solved thrust graph.
The imported profile plane and roll are preserved. When a profile vertex matches
the guide start within 1e-5 model units, it becomes the geometric anchor. Otherwise
the C-Tube vertex-centroid origin convention applies. Material centroids are always
calculated from area, independently of that geometric anchor.

The sketch defaults to solid material for its closed Rhino generatrix. Enable
Thin-wall profile explicitly for a wall centreline. In thin-wall mode, thickness
is fixed, so area grows linearly with profile scale. For a solid polygon, area
grows quadratically. The HUD reports the seed linear scale; the inspector also
reports input material area and first target area to make the conversion inspectable.

## Force conversion and objective

For each member, partition actual physical face area by source-magnitude shares.
Ignore virtual faces. Average the available two end contributions; a support end
without a cell uses the other end. Retain their absolute difference as a separate
force-area compatibility check. A shared node section uses the largest adjacent
member demand, not the total force-cell surface area.

`force = memberFaceArea / forceScale`

`targetSectionArea = force / referenceStress`

Units must be consistent. Default normalised loads/reference stress are exploratory
quantities until calibrated. This is an equality-fit objective, not minimum-weight
structural design: matching axial area does not guarantee an eccentric section passes.

Optimisation holds section planes and guide origins fixed. It changes the seed
scale and one positive apex scale per segment through
`zTsCTubes::computeTubeWithFrames`. Every trial uses the existing projection
construction rather than independent vertex/section scaling. A bounded, damped
least-squares iteration minimises log area ratios and second differences of log
area. Settings bound each control, not the cumulative tube envelope. Report the
maximum relative area error, termination and dimensionless strip-quad warp.
Bounds/smoothing/iteration limits can prevent meeting area tolerance; do not
classify such a result as resolved. No unroll/self-intersection certificate is implied.

## Eccentricity and section checks

Intersect finite segments of the thrust graph with each physical section plane.
No intersection or multiple distinct intersections leaves that section unresolved.
Coincident intersections at a load node use the larger adjacent normal-force demand.

The eccentricity vector connects the material area centroid to the intersection
within the plane. Compute the section's full second-moment matrix and the linear
normal-stress field for the normal force and its eccentric moment. Transverse
force, shear/torsion, buckling, crushing/joint models and cracking redistribution
are outside this check. Tension members fail the compression-only check.

`kernUtilisation = max(0, 1 - minStress / (normalForce / materialArea))`

Utilisation ≤ 1 means the linear normal-stress field has no tension; also require
`maxStress / referenceStress ≤ 1`. This recovers the middle-third boundary for a
rectangle under uniaxial eccentric compression without imposing `t/6` on C profiles.

## Display and lifecycle

The translucent tube uses the common green–yellow–red status domain. Centroids
are blue, with coloured segments to thrust intersections. The HUD adds section
area deviation, compression-kern and normal-stress utilisation to the force angle
diagnostic. Detailed input/model/warp information is in the tube inspector.

Changed force geometry/areas/signs or physical sizing settings mark the tube
**Update required**, with grey status colours. Regenerate/optimise explicitly.
Display/section-area tolerance edits update immediately without running sizing
or resetting AP/3DGS. Reimported guides clear derived tube results. Loads remain
frozen; tube-volume self-weight feedback and branching guide junctions are deferred.

## Verification

Native smoke test `zspace_toolsets_force_ctubes` covers analytical rectangular
area, middle-third boundary/stresses, thin-wall area, split branched-section seams,
nonuniform area fitting, projection-quad planarity and invalid/missing inputs.
Viewer `tests/forceCTubes.test.mjs` exercises the installed typed runtime and the
production sketch, including explicit generation and stale-result handling.

## Mesh Force C-Tubes: physical sizing (October 2026)

The Mesh Force C-Tubes sketch now separates force-diagram display scale from
physical sizing. Its calibration input kF has units N per raw diagram-area unit:

    designForce = rawMemberFaceArea * kF * designLoadFactor
    requiredMaterialArea = designForce / designStress

Design stress and elastic modulus use N/model-unit squared; all geometry and
explicit effective buckling lengths use the same model length unit. For example,
when model units are mm, N/mm² equals MPa. The default calibration, factor and
stress are 1 and must be set to the project's units/material. A known physical
load divided by its corresponding raw force-face area can establish kF; visual
cell scale is never a force calibration.

Proportional sizing uses one common scale for all source widths and the world-Z
section depth. Wall thickness stays fixed. For this centreline thin-wall model,
A scales with s and I with s³. The governing scale satisfies axial stress and,
when supplied, ideal Euler flexural buckling using pi² E Imin / Leff². Excess
material area is safe for this check and is not an area-matching error. Width
bounds may leave capacity unresolved. Independent-width mode retains its area
fit, then applies any additional common capacity scale. Loads remain frozen.

The paired half-C profile is integrated as six midline wall segments. Its
centroidal second-moment matrix is projected normal to the associated graph
member; its smaller principal inertia supplies the preliminary Euler check.
Ordinary shared cuts use the largest adjacent demand; each junction member has
its own paired cut, sharing the physical boundary variables. Compression uses
Euler; tension does not. After directrix fitting, inertia is rechecked and an
additional common scale is applied if required. Actual panel warp and seam gaps
are reported separately. The graph remains at section mid-depth.

An unset modulus or effective length leaves compression buckling pending.
Effective length is an explicit restraint-based input, never the graph's edge
sampling distance. A single supplied length is applied to all members. These
are preliminary axial and ideal Euler checks: they assume the two halves act
together and neglect junction overlaps, eccentricity effects, imperfections,
local/distortional/torsional buckling and code-specific resistance reductions.
They do not establish junction adequacy or fabrication developability.

Changing display scale does not stale or resize the tube. Changing force
calibration, material data or effective length marks sizing for an explicit
update without resetting the AP/3DGS history. The HUD shows axial utilisation,
Euler utilisation when available, and panel warp alongside force convergence.

The capacity exceedance tolerance permits utilisation up to 1 + tolerance; set it to zero for strict unity acceptance. The HUD uses the same threshold.
