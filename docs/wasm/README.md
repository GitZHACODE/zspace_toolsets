# Toolsets browser runtime

Build with `scripts/build-wasm.bat` from this folder, or
`wasm/scripts/build-wasm.bat` from the repository root. `scripts/build-wasm.ps1`
accepts `-CoreDir` and `-EmsdkDir`. Generated build/cache/tmp/out files stay ignored.

The toolsets bridge compiles against the current sibling core sources. Its API
and heap are independent of the core viewer runtime. Follow
[typed transfer contract](../docs/architecture/typed-wasm-transfers.md).

Preserve existing C exports; add transfer APIs without changing the legacy ABI.
Refresh the viewer's packaged JS/WASM pair together only when runtime integration
is in scope. Read applicable viewer instructions before editing adapters.
