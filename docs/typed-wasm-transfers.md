# Toolsets WASM transfers

The current bridge provides transfer version 1, advertised by
`zspace_toolsets_transfer_version()`. Its legacy ABI version remains 1.

## Typed path

- `zspace_alloc(bytes)` / `zspace_free(pointer)` match the core allocation contract.
- C-Tubes compute, branching optimization, Vault fit and panelisation compute have
  additive `_from_buffers(handle, pointer, byteLength)` exports.
- Diagnostics, branching results, Vault results/diagrams and panelisation results
  have `_buffers` exports. Read `zspace_toolsets_result_buffer_size(handle)` and
  copy the returned bytes immediately, before the next result read/mutation.
- Mesh positions/counts/connectivity retain their existing typed pointer exports.

The ZST1 packet starts with little-endian `0x3154535a`. A tagged tree preserves
object/array structure, strings, booleans and numeric values; numerical leaves
are contiguous signed int32 or float64 blocks. Float64 avoids precision loss
for force densities, loads and input values before native coordinate validation.
Signed integers preserve negative sentinels and reject invalid indices instead
of wrapping them into uint32. Numeric array results remain public `number[]`
for compatibility; mesh payload results remain owned typed arrays.

Limits: 256 MiB per packet and nesting depth 64. Both decoders validate tags,
collection lengths, truncation and trailing bytes; the native entrypoint checks
that the input range is inside WASM memory. The existing toolset-specific size,
coordinate and index checks still apply. This change does not raise those limits.

The viewer rejects non-finite numerical inputs before allocation and uses the shared core allocation/copy/free helper, reading the heap
after allocations in case memory grows. The branch worker transfers input and
result ArrayBuffers, without JSON stringify/parse on the typed path.

## Compatibility and scope

The old file/string APIs remain callable. Older installed modules use an explicit
JSON fallback; current modules must use the typed path. Do not remove the fallback
until a separately approved ABI migration. The binary path does not use the
Emscripten filesystem for uploads or JSON text dump/parse for result transfers.

Internal native result construction and validation still use the existing Json
value model. Binary numeric values are decoded into it; zero-copy native numeric
views are a possible later optimization. This migration removes text conversion
and unifies transport, not every native allocation. No algorithm or topology
behavior was changed.

## Verification

Build with `wasm/scripts/build-wasm.bat`. Run the existing toolsets WASM scripts
for legacy API compatibility. After explicitly refreshing the matching viewer
JS/WASM pair, run `npm run test:focused -- tests/toolsetsTypedTransfers.test.mjs`,
`npm run build` and `npm test` in the viewer. Parity checks cover all toolsets,
mesh readback, result/diagram fields, rollback, malformed packets, allocation
failures and heap growth. Browser worker behavior also needs a browser check
when its execution path changes. Never rewrite historical reports during tests.
