# Wall panelisation

`zTsPanelisation` runs in the shared toolsets WASM module, following the C-Tubes
session/bridge architecture. The sketch packs Rhino polygons and publishes results;
all chart construction, stitching, clipping and panel adjacency run in C++.

## Geometry and dimensions

Input is triangle/quad meshes in local metres, with explicit crease endpoint pairs.
The viewer uses refined, limit-sampled SubD polygons and saved BREP render meshes,
retaining source object IDs, layers and the common origin.

1. Stitch geometrically corresponding complete boundary edges. Interior points are
   never position-welded. Ambiguous overlapping skins and collapsed-pole components
   stay separate; roll back joins that would create non-manifold edges.
2. Preserve Rhino creases and detect additional 35-degree folds. Grow charts through
   smooth face adjacency. Quad charts propagate opposite-edge row directions, so
   their coordinates follow walls, soffits and curved canopies without world-axis
   seeds or radial distance wavefronts. Extraordinary/periodic topology creates seams.
3. Integrate average physical edge length for each quad-row interval. X follows the
   dominant crease-rail family (the longer family when uncreased); Y crosses it.
   Mixed triangle/quad patches use least-squares conformal coordinates, area-scaled
   to metres and aligned with their boundary directions. Strongly turning mixed
   patches receive automatic seams; collapsed faces do not enter the linear solve.
4. Clip both intrinsic coordinates at X and Y intervals. Cache intersections by
   original edge and edge parameter. At chart seams, insert the union of both sides'
   cuts in both polygon boundaries, preserving shared vertex/edge topology without
   cracks or T-junctions. Chart transitions can terminate a grid line at a seam;
   this is not a globally seamless UV grid on arbitrary branching topology.
5. Group connected fragments by chart, X station, Y column and source object,
   stopping at creases. End remnants and small panels are retained, preserving area.

X/Y are nominal surface dimensions on the polygon approximation. Quad coordinates
use row-average lengths; conformal coordinates distribute curvature distortion.
Doubly curved surfaces therefore do not guarantee every panel edge equals X or Y.
This is not exact NURBS arc length or a fabrication-tolerance guarantee.

Output includes the mesh, face-to-panel/source maps, joints and panel records
(source, station, column, chart, area). Joint kind 1 is X/blue, 3 is Y/cyan,
2 is a crease/chart/source seam/red, and 0 is an outer boundary/grey.
Small areas are highlighted, never removed. Panels remain zero thickness;
physical gaps, substructure and service collision checking are subsequent stages.

## Viewer controls and lifetime

The sketch exposes only **Panel length X (m)**, **Panel width Y (m)** and
**Show services (0 / 1)**. Compute panels updates the layout on all WALLS objects.
Show source restores the layered surface inputs. Old saved axis, seed, crease,
representation and wall-index UI settings cannot override these automatic choices.

The shared WASM build reserves a 2 MB stack for Eigen sparse ordering and enables
stack overflow checks; the default 64 KB fails on the full-model chart solve.

The low-level kernel retains offset/reverse and tolerance options for callers;
legacy axis, seedBand and surfaceDistance fields no longer choose a distance field.
Failed computes leave the last successful session intact. The bridge uses shared
create/destroy/error/bulk-buffer functions and copies output arrays between heaps.
The sketch's `panelisationResult()` retains source IDs, settings and maps for the
next fabrication stage. Generated runtime assets remain ignored.

## Verification

```powershell
cmake --build build/ninja-msvc --config Release --target zspace_toolsets_panelisation
ctest --test-dir build/ninja-msvc -C Release -R "panelisation|ctubes|vault" --output-on-failure
./wasm/scripts/build-wasm.ps1
node tests/panelisation_wasm.mjs
node tests/ctubes_wasm.mjs
node tests/vault_wasm.mjs
```

Native checks cover independent X/Y dimensions, rotated surfaces, bent ribbons,
area conservation and seam continuity. WASM fixtures additionally cover a smooth
arch, mixed triangle/quad walls, shared source boundaries, atomic failure and all
27 WALLS objects in the local Rhino asset when present. The native full-model
result at X=Y=1 m has 2,929 panels and preserves approximately 795.651 m².
The optional second argument to the native smoke executable exports a JSON preview.

In the viewer run `node --test tests/panelisation.test.mjs tests/rhinoSubD.test.mjs`
and `npm run build`. Browser availability is independent of those integration tests.

## Previous planar arch-wall method (superseded in review)

The viewer assigns an optional orthonormal `panelFrame` only to the three Rhino
patch IDs listed in `orangeWallReview.wallIds`. U is signed horizontal distance
from the canopy crown-rail centreline; V is height above Rhino world Z=0. X stations
therefore repeat on both sides of the centreline, and Y rows share exact horizontal
levels through the wall, return and boundary fillets. Original geometry and crease
boundaries are retained. This intentionally uses horizontal/vertical cuts for this
patch, overriding surface-distance stationing only here. The canopy/crown rail
keep their original intrinsic charts (regression-tested for identical panel records).
The authored frame is for this supplied model/local origin and is not a generic
symmetry detector. New geometry patches should receive separately reviewed methods.

## Current topology-flow guide review

`topologyGuides` requests a non-cutting preview on the three orange wall/fillet
patches. Each regular quad chart selects its highest boundary side, seeds it at
equal arc-length intervals no greater than X, and transports the seed parameters
through opposite quad edges. At every transverse mesh row, measure surface-polyline
length between adjacent paths. If it exceeds X (with a small numerical tolerance),
insert enough equally spaced intermediate paths to bring that sampled gap below X.
The inserted paths follow the same topology thereafter; their birth is interpolated
between the previous and current rows at the width threshold. Blue denotes initial
paths; red denotes inserted paths. Crease/chart boundaries remain explicit, and
non-quad charts are counted as unsupported rather than traced through silently.

The guide graph is returned separately from the panel mesh. The orange review shows
it over the original wall geometry, without the rejected vertical/horizontal panel
grid. The canopy panel records are unchanged. Y is retained for later panel cuts,
not applied to these guides. These are review curves with free ends, not fabricated
panel boundaries; the width bound is row-sampled, not a continuous-surface guarantee.
The fixed-frame option remains available to low-level callers but is no longer
assigned by the viewer's orange review profile.

## Flow panels with Y cross-splits

Topology-enabled quad charts now use the traced path parameters for their U cuts.
V integrates mesh-edge distance down each topological column from the selected top boundary,
and cross-cuts occur every Y metres of that surface coordinate. Adaptive branches
are extended upstream to the preceding Y joint, creating closed panel regions
rather than interior slits. Cuts are clipped against the original triangles;
original-edge intersections and interior branch/Y junctions are inserted on both
sides. The widening-fan regression checks for unmatched internal edges, area loss
and the effect of reducing Y. Other charts, including the canopy, keep their
previous algorithm. Lengths remain mesh/row-average approximations.

The viewer publishes each orange panel as a separate mesh, with stable panel/source
mapping in `panelisationResult().panelObjects`. Material/source layers remain intact.
The review action now displays those mesh patches and actual X/Y panel joints;
it does not overlay the earlier free-ended red guide preview. The raw guide data
remains in the result for diagnostics. This supersedes the guide-only review above.

## Curved wall divider

The orange main wall profile selects a surveyed edge of the existing curved quad
rail outside the arch using Rhino surface vertex IDs 2120/3582 (refinement 2).
`panelSeamSeeds` carries mesh endpoint pairs to C++. The kernel extends each seed
through opposite incident edges at regular quad vertices, stopping at boundaries
or irregular topology. These dividers do not split/restart the intrinsic charts:
X paths and Y stationing remain unchanged. They prevent panel merging across the
rail after clipping. The existing union of edge cuts inserts each X/Y intersection
on both incident faces, producing shared panel corners instead of overlay curves.
Joint kind 4 identifies this divider; the sketch displays it in orange.
The profile rejects missing seed IDs or changed refinement for reselection.
WASM regressions check a synthetic fan and the supplied model for complete rail
length, two-ended continuity, two incident faces, separate panel IDs and shared Y
vertices. Canopy panel records remain unchanged.

### Y termination before the divider

The revised review removes Y joints that reach the curved divider. Adjacent panel
regions are merged across the complete affected Y boundary, so the joint stops at
the preceding flow edge and the panel beside the divider is uninterrupted. This
changes panel membership, not just line visibility; internal tessellation remains.
The divider and canopy remain intact. Regression checks now require no Y-joint
vertices on the divider at X=1 with Y=1 and Y=2, plus the synthetic fan case.

### Crease-first top stationing

Flow charts anchor top starts to crease edges entering the wall, including partial
creases that terminate inside a chart. Each span between anchors is divided into
equal surface lengths no greater than X (within numerical tolerance). Top boundary
creases do not turn every tessellation vertex into a seed. Uncreased spans retain
equal spacing. The partial-crease regression checks a rail at 0.95 m does not gain
an adjacent 1 m cut and its 0.05 m sliver. Curved dividers and Y termination are unchanged.

### Continuous rail colours and cyan cleanup

Crease-seeded flow paths retain their red classification beyond the authored
crease endpoint and through the curved divider. For the reviewed arch-wall source,
Y boundaries are removed across a whole neighbouring panel pair when they enter
the bending region: a sampled column has turned ten degrees from its top tangent.
Cyan remnants shorter than 10% of min(X,Y) are also merged on this source. Both
rules merge real panel regions, retain source geometry, and leave canopy charts
unchanged. The colour regression follows a partial crease through its continuation.

Cleanup joins connected Y-boundary spans across blue flow rails, but does not cross
red crease rails or source/chart boundaries. This places remaining cyan endpoints
on red rails. A supplied-model regression checks every cyan endpoint on the main
wall at X=1/Y=2 against the red-rail/source-boundary vertex set.

## Ceiling and staggered arch review

The source-only sketch profile now assigns `layoutMode=1` to the long ceiling and
its edge trim, with `alignmentSource` pointing to the reviewed main wall and a
surveyed horizontal `alignmentDirection`. Mode 2 applies to the arch canopy,
with `staggerFraction=0.5`. These optional fields travel through the shared bridge;
the geometry work remains in zTsPanelisation.

Both modes build continuous quad charts across crease rails, preserving those
creases as final panel boundaries. The ceiling uses a common metric centreline for
its skin and trim; projected stations turn continuously through the curved reach.
Only the wall's highest top starts anchor the aligned reach. Intervening and end
spans are equally subdivided using X, and Y subdivides between existing rails.
This is stationing on the sampled mesh, not exact NURBS arc length.

The arch retains the authored crease-bounded loops. Odd loops have half-X phase;
each target is snapped to its nearest intrinsic topology row and duplicate rows
are removed. Thus X is a target distance for the arch, and coarse topology can
produce longer or shorter panels. The union of seam subdivisions keeps staggered
neighbours conforming. Source/layer IDs are retained, and each reviewed ceiling,
arch and wall panel is separately selectable. The three approved wall layouts
remain unchanged. Synthetic regressions cover exact straight-ceiling wall station
alignment, alternating phases and non-integral X snapping; full-model tests cover
area, wall preservation and the prior curved-divider/cyan rules.

For the ceiling reference, only full-width topology rows determine the centreline;
partial rows around vertical end returns are excluded and the two straight
terminal tangents extend to cover the source extents. Cubic interpolation blends
tangent directions through the corner. Station distance is measured along this
horizontal ceiling reference, while every panel vertex remains on the original
sampled surface. This prevents end returns from skewing the ceiling's plan cuts.

## Clean section cuts (current review)

The ceiling restores its shared centreline and wall-aligned stationing, with
half-X stagger on alternate broad panel bands. Narrow crease returns do not
consume the alternating phase. Ceiling cuts use section planes normal to the
horizontal centreline tangent, creating straight plan lines through the fan.

Arch stations still snap to the nearest metric topology row, but their section
planes are perpendicular to the local arch tangent in the surveyed elevation
plane. This corrects the inclination carried by some original mesh rows.
Both modes clip source triangles against these exact planes; they no longer
contour an interpolated distance field or require whole-face-only partitions.
Vertices remain on the sampled surface, and shared edge intersections are reused.

Regression checks measure the actual transverse panel boundaries for straightness
in ceiling plan and arch elevation, alongside surface-area conservation, wall
layout preservation, stagger and wall-alignment checks.

## Equal ceiling rail divisions (latest review)

The ceiling now extracts ordered longitudinal edge rails from its quad chart.
Complete native rails include the attached vertical return down to ground. Every rail is divided by the common count
`n = ceil(shortest rail length / (2 * X))`, as written in the review sketch.
Adjacent broad strips use even and odd internal division indices respectively;
narrow crease returns retain the adjacent strip phase. Both endpoints of each
cross-cut use the same arc-length fraction on their respective rails. A section plane through those points follows the local surface direction, including
the ground return, while retaining
the sampled surface. This replaces ceiling wall-aligned stationing. Arch and
vertical-wall strategies remain unchanged. With alternate indices, the resulting
panel length is about twice a rail's division interval; X is the formula input,
not a ceiling maximum-width constraint in this trial.

Straight ceiling reaches now lock section normals to the surveyed longitudinal
or transverse axis (both rail tangents within 10 degrees); the ground return
likewise uses vertical stationing. Stations remain centred on the two rail
samples and retain their alternating phase. Rail-to-rail inclined sections apply
only through the bend. A skewed-endpoint regression checks that equal rail
fractions cannot tilt the straight-reach cuts.

Ceiling rows now share cuts in pairs: A,A,B,B. Each pair uses one set of section
planes derived from the pair's outer rails, guaranteeing the split meets at the
shared boundary. The next pair uses the alternate division indices. Narrow
crease returns do not advance the row count. Arch staggering is unchanged.

Latest review removes row pairing and restores A,B,A,B ceiling staggering.
Small charts at extraordinary vertices on the main ceiling inherit the cutting
planes of the closest main-chart band. They no longer apply the full division
count to their short local rails, avoiding concentrated tiny panels at the end.
Straight-reach parallelism and full-length rails to ground are retained.

The front arch now uses the same equal arc-length rail divisions as the ceiling.
Its shared count is ceil(shortest rail / (2X)), raised by one if even. Adjacent
crease-bounded strips use alternate division indices (A,B). Cuts are section
planes in the surveyed arch elevation; straight reaches retain parallel cuts.
Small charts inherit the main arch stations instead of restarting the division
count. The ceiling count and the approved vertical wall layout are unchanged.
