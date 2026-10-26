# Toolsets WASM work

- Preserve dirty/concurrent changes and the existing C and toolset APIs.
- Use the typed buffer transfer contract in
  `../docs/architecture/typed-wasm-transfers.md` for bulk numeric data. Keep legacy
  JSON entrypoints compatible; do not route current viewer sessions through them.
- Share allocation/copy/free behavior with the viewer's core transfer helper.
  Reacquire heap views after allocation, free inputs in finally, and own readback.
- Keep numerical validation and toolset limits; do not wrap signed invalid indices.
- Build this runtime and run affected legacy and typed parity tests. Do not rebuild
  every SDK or run large benchmarks for unrelated changes.
- Generated outputs are ignored. Do not commit/push unless requested.
