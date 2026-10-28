# Slicing WASM array ABI 2

The combined runtime is zspace_toolset.js / zspace_toolset.wasm.
Core owns its bindings, allocator and build configuration; Toolsets adds Slicer
bindings. The independent Core runtime remains available.

## Data flow

JavaScript TypedArrays ? Core allocator / WASM heap ? native numeric vectors ?
Slicer ? native result vectors ? owned JavaScript TypedArrays.

There is no JSON value tree, JSON text, tagged packet, encoding or decoding step
in the Slicer bridge. Worker messages contain numeric TypedArrays and small
structured-clone envelopes. ArrayBuffers are transferred between Workers and
the UI, not serialized to JSON. File bytes use the Emscripten filesystem for
Toolsets USD IO.

Arrays use contiguous float64 values at this boundary. Coordinates, fields and
attributes retain numeric precision. Indices are checked for finite integral
values and range before conversion to native integers; negative width-group
sentinels are preserved. This follows Core's pointer/count convention; it does
not require Core's float32 position storage format.

## Session and input

- zspace_slicing_abi_version() returns **2**.
- zspace_slicing_create() / zspace_slicing_destroy(handle) manage ownership.
- zspace_slicing_set_array(handle, slot, pointer, count) copies float64 data.
- zspace_slicing_compute(handle, mode, parallel, processed) computes;
  mode 0 is Topology, mode 1 is Harmonic.
- zspace_toolsets_last_error() supplies the last error string.
- zspace_alloc(bytes) / zspace_free(pointer) are provided by Core.

| Input slot | Values |
| --- | --- |
| 0 | Mesh XYZ positions |
| 1 | Polygon vertex counts |
| 2 | Polygon vertex indices |
| 3 | Per-face tags |
| 4 | Per-face width IDs, optional |
| 5 | Slot anchor XYZ |
| 6 | Bracing XYZ positions |
| 7 | Bracing edge index pairs |
| 8 | Seventeen scalar settings in the order below |
| 9 | Boundary width table |
| 10 | Per-curve bracing widths |

Settings order: layerSpacing, fieldCellSize, fieldScaleUp, bracingWidth,
overlapWidth, doubleBoundary, slotOffset, bracingSlot,
traceStepFraction, traceMaxSteps, equaliseIterations, arapIterations,
unrollStart, resampleLength, featureAngle, seamOffset, seamSize.
Boolean settings use 0/1. The JS adapter fills omitted settings with the existing
defaults. USD import supplies the resolved native settings.

## Results

Read zspace_slicing_record_count(handle) and each record's
zspace_slicing_record_kind(handle, record) / zspace_slicing_record_layer.
Then use zspace_slicing_array_ptr(handle, record, slot) and
zspace_slicing_array_count for the numeric data.

| Kind | Contents | Array slots |
| --- | --- | --- |
| 0 | Shared field mesh | 0 XYZ, 1 face counts, 2 connectivity |
| 1 | Combined print path | 0 XYZ, 1 edges, 2 normals, 3 widths, 4 heights, 5 track IDs, 6 original side FaceTags |
| 2 | Layer slice mesh | Same as kind 0 |
| 3 | Layer chart mesh | Same as kind 0 |
| 4 | Layer print path | Same as kind 1 |
| 5 | Layer scalar fields | 0 path field, 1 width field |

Global records use layer -1. Per-layer records use zero-based indices.
The graph and attribute order is preserved without sorting.

zspace_slicing_result_number(handle, field) returns: 0 compute seconds,
1 mode, 2 actual thread count, 3 layer count, 4 parallel flag,
5 postprocessed flag. zspace_slicing_status(handle) returns status text.
Algorithm timing excludes result extraction, transfers and display.
The viewer reports algorithm time, data exchange, and their sum; it excludes
runtime initialization and display.

Pointers belong to the native session and expire on mutation/destruction.
After any allocation, reacquire the heap view because memory may grow.
The adapter immediately copies results into owned Float64Arrays.
An invalid compute preserves the previous successful result.
Count queries return -1 on failure; mutating operations return 0 on failure.

## USD IO

zspace_slicing_import_usd(handle, file) calls Toolsets IO backed by Core USD.
Read its arrays using zspace_slicing_input_ptr/count(handle, slot).
Unknown-setting warnings use zspace_slicing_warning_count / warning.

For export:
1. zspace_slicing_export_begin(handle, layerCount).
2. Upload each path using input slots 0 XYZ, 1 edges, 2 normals, 3 widths,
   4 heights, optionally 6 FaceTags (one integer 0..6 per point).
3. zspace_slicing_export_add_path(handle, layer), using -1 for the combined
   path and zero-based layer indices for individual paths.
4. zspace_slicing_export_usd(handle, source, output, blockid).

The export uses the current Toolsets IO implementation, including Gcode,
PrintMesh and print-path prims. No USD schema parsing is duplicated in JS.

## Migration and validation

The old Slicer compute_from_buffers, result_buffers, transfer-version and
packet-size exports have been removed. Deploy the matching adapter and JS/WASM
pair together. Other Toolsets runtimes and existing Core bindings are unchanged.

Build with wasm/scripts/build-slicing-wasm.ps1.
Viewer coverage is in tests/slicingWasm.test.mjs,
tests/slicingNativeIO.test.mjs and tests/slicingUsdWasm.test.mjs.
Checks cover both algorithms, serial/parallel parity, USD roundtrips including
PrintMesh, invalid pointers/indices, result ownership and Worker cancellation.
