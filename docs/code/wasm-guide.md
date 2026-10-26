# Toolsets WebAssembly and buffer ownership

## Toolsets WebAssembly and buffer ownership

The toolsets browser runtime is its own zspace\_toolsets.js/.wasm pair. Core and Toolsets do not share heaps or native pointers. Viewer adapters translate scene geometry into the receiving module’s input representation and copy the returned geometry into JavaScript-owned arrays.

The bridge reference covers the main toolsets bridge and the separate Force C-Tubes/Force Polyhedra bridge source files. Session handles, allocation/free functions and typed packet exports must be paired within the same module. Dispose a session when the owning UI/sketch is removed.

Typed geometry transfers use the ZST1 version 1 packet contract where supported. Keep older JSON paths only for runtimes that lack the typed capability. Allocate, refresh heap views after allocation, copy data, invoke the export, read copied results and free temporary storage in finally blocks. Treat diagnostics as part of the result, especially fit/convergence status.

[Toolsets bridge export reference](group-bridge/)

[Viewer toolsets adapter](https://app.gitbook.com/s/lxNYT4JxjtqJ22h7b4hp/code/group-kernels/module-src-kernels-zspacetoolsetswasmadapter)

[Core bridge ownership](https://app.gitbook.com/s/FKPJDjBzZtKzxuy72aOX/code/wasm-guide)

[Typed packet layout and compatibility](../typed-wasm-transfers.md)

[wasm/README.md](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/wasm/README.md#L1)
