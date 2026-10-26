# Project Toolsets

zspace\_toolsets contains project-specific tools and algorithms built on zspace\_core. It owns the higher-level construction, optimisation and fabrication workflows, while core provides shared geometry types, function sets and IO.

### Responsibilities

Maintained C++ toolsets are grouped under include/zspace/zToolsets and src/zToolsets. Existing documentation covers 3DP slicing, C-Tubes, force polyhedra, Vault force-density fitting, panelisation and parametric stairs. Legacy pre-v3 code is reference material.

Native builds link the core, interface and IO libraries. The browser bridge builds a separate zspace\_toolsets JavaScript/WASM pair against the sibling core sources; it does not share a heap or pointers with the core module loaded by the viewer.

### Relationship to the other repositories

[Geometry Kernel](https://app.gitbook.com/o/0Lv7pgeoKI2MXprc3zOv/s/FKPJDjBzZtKzxuy72aOX/) defines the common geometry objects and operations used by each toolset. Keep reusable geometry primitives and general kernel changes there.

[Web Viewer](https://app.gitbook.com/o/0Lv7pgeoKI2MXprc3zOv/s/lxNYT4JxjtqJ22h7b4hp/) supplies project inputs, parameter controls and interactive presentation. Its adapters send snapshots and buffers into toolset sessions and copy results back for display.

### Follow a project workflow

Read [C-Tubes initial port](ctubes-port.md) for the C++ toolset boundary, [Toolsets WASM transfers](typed-wasm-transfers.md) for transport, and [Force C-Tubes viewer workflow](https://app.gitbook.com/s/lxNYT4JxjtqJ22h7b4hp/docs/force-ctubes) for interactive use.

The [kernel graph ribbon implementation](https://app.gitbook.com/s/FKPJDjBzZtKzxuy72aOX/docs/architecture/graph-mesh) is a shared geometry operation; toolset-specific fitting and project constraints remain in this repository.

[Source repository and review branch](https://github.com/GitZHACODE/zspace_toolsets/tree/docs/gitbook-site-sync)

### Code documentation

Read the architecture and ownership guides, then use the source-derived API reference. The source-version page records the code revision and generation scope.

[Code Documentation](code/)

[Architecture and API ownership](code/architecture-guide.md)

[Source versions and generation](code/generation.md)
