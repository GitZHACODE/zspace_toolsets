# zspace_toolsets v3 migration

The pre-v3 Visual Studio/header-only style tree has been moved to `legacy/pre-v3` and is kept as algorithm reference material.

New toolsets should follow `zspace_core` structure:

```text
include/zspace/zToolsets/<domain>/   Public declarations
src/zToolsets/<domain>/              Implementations
examples/<domain>/                   Small usage examples
tests/smoke/                         Compile/runtime smoke tests
scripts/                             Build helpers
```

Rules:

- Use modern zSpace core objects such as `zObjectMesh`, `zObjectGraph`, and `zObjectMeshScalarField`.
- Use `zFn*` function sets for create/edit/query methods.
- Use `zIO` for mesh and graph import/export.
- Keep legacy `zObj*` usage inside `legacy/pre-v3` only.
- Add one focused toolset at a time and keep old code out of the default build.
