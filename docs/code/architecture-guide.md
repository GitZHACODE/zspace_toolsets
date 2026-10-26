# Architecture and method ownership

zspace\_toolsets builds focused project methods on top of zspace\_core. Maintained declarations are under include/zspace/zToolsets/ and implementations under src/zToolsets/. legacy/pre-v3 is historical reference and is excluded from this API reference.

The native build locates its sibling Core through ZSPACE\_CORE\_DIR, defaulting to ../zspace\_core, and links the kernel modules it uses. Toolsets accept and return Core geometry objects; they own project-specific settings, staged computation and reports.

### Find the method family

Fabrication includes 3D-print slicing and synthesis, C-Tubes and branching tubes. Statics includes Vault Best Fit FDM, Force Polyhedra and Force C-Tubes. Architecture includes panelisation and parametric stairs. Each class page lists its current inputs, settings, compute methods, result accessors and report types.

[Fabrication APIs](group-fabrication/)

[Statics APIs](group-statics/)

[Architecture APIs](group-architecture.md)

[Core geometry ownership](https://app.gitbook.com/s/FKPJDjBzZtKzxuy72aOX/code/architecture-guide)

[Viewer integration](https://app.gitbook.com/s/lxNYT4JxjtqJ22h7b4hp/docs/code/kernel-guide)

[CMakeLists.txt](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/CMakeLists.txt#L1)

[README.md](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/README.md#L1)
