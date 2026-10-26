# zspace\_toolsets\_wasm\_bridge.cpp declarations

## wasm/bridge/zspace\_toolsets\_wasm\_bridge.cpp

[Source file](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/wasm/bridge/zspace_toolsets_wasm_bridge.cpp#L1)

Free functions, aliases, enums and constants declared in this file. Class and struct members have separate pages.

#### zspace\_toolsets\_abi\_version

{% code collapsedlinecount="10" %}
```cpp
int zspace_toolsets_abi_version()
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/wasm/bridge/zspace_toolsets_wasm_bridge.cpp#L92)

#### zspace\_toolsets\_transfer\_version

{% code collapsedlinecount="10" %}
```cpp
int zspace_toolsets_transfer_version()
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/wasm/bridge/zspace_toolsets_wasm_bridge.cpp#L93)

#### zspace\_alloc

{% code collapsedlinecount="10" %}
```cpp
void * zspace_alloc(size_t size)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: size. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/wasm/bridge/zspace_toolsets_wasm_bridge.cpp#L94)

#### zspace\_free

{% code collapsedlinecount="10" %}
```cpp
void zspace_free(void *pointer)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: pointer. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/wasm/bridge/zspace_toolsets_wasm_bridge.cpp#L95)

#### zspace\_panelisation\_create

{% code collapsedlinecount="10" %}
```cpp
int zspace_panelisation_create()
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/wasm/bridge/zspace_toolsets_wasm_bridge.cpp#L96)

#### zspace\_panelisation\_compute

{% code collapsedlinecount="10" %}
```cpp
int zspace_panelisation_compute(int handle, const char *path)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: handle, path. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/wasm/bridge/zspace_toolsets_wasm_bridge.cpp#L99)

#### zspace\_panelisation\_result

{% code collapsedlinecount="10" %}
```cpp
const char * zspace_panelisation_result(int handle)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: handle. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/wasm/bridge/zspace_toolsets_wasm_bridge.cpp#L127)

#### zspace\_panelisation\_select\_mesh

{% code collapsedlinecount="10" %}
```cpp
int zspace_panelisation_select_mesh(int handle)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: handle. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/wasm/bridge/zspace_toolsets_wasm_bridge.cpp#L130)

#### zspace\_toolsets\_last\_error

{% code collapsedlinecount="10" %}
```cpp
const char * zspace_toolsets_last_error()
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/wasm/bridge/zspace_toolsets_wasm_bridge.cpp#L135)

#### zspace\_ctubes\_create

{% code collapsedlinecount="10" %}
```cpp
int zspace_ctubes_create()
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/wasm/bridge/zspace_toolsets_wasm_bridge.cpp#L136)

#### zspace\_stairs\_create

{% code collapsedlinecount="10" %}
```cpp
int zspace_stairs_create()
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/wasm/bridge/zspace_toolsets_wasm_bridge.cpp#L157)

#### zspace\_stairs\_update

{% code collapsedlinecount="10" %}
```cpp
int zspace_stairs_update(int h, const char *path)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: h, path. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/wasm/bridge/zspace_toolsets_wasm_bridge.cpp#L160)

#### zspace\_stairs\_optimise\_graph

{% code collapsedlinecount="10" %}
```cpp
int zspace_stairs_optimise_graph(int h, const char *path)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: h, path. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/wasm/bridge/zspace_toolsets_wasm_bridge.cpp#L167)

#### zspace\_stairs\_optimised\_graph

{% code collapsedlinecount="10" %}
```cpp
const char * zspace_stairs_optimised_graph(int h)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: h. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/wasm/bridge/zspace_toolsets_wasm_bridge.cpp#L175)

#### zspace\_stairs\_step

{% code collapsedlinecount="10" %}
```cpp
int zspace_stairs_step(int h, int iterations)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: h, iterations. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/wasm/bridge/zspace_toolsets_wasm_bridge.cpp#L178)

#### zspace\_stairs\_reset

{% code collapsedlinecount="10" %}
```cpp
int zspace_stairs_reset(int h)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: h. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/wasm/bridge/zspace_toolsets_wasm_bridge.cpp#L181)

#### zspace\_stairs\_snap

{% code collapsedlinecount="10" %}
```cpp
int zspace_stairs_snap(int h, double geometryWeight)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: geometryWeight, h. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/wasm/bridge/zspace_toolsets_wasm_bridge.cpp#L184)

#### zspace\_stairs\_previous

{% code collapsedlinecount="10" %}
```cpp
int zspace_stairs_previous(int h)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: h. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/wasm/bridge/zspace_toolsets_wasm_bridge.cpp#L187)

#### zspace\_stairs\_smooth\_after

{% code collapsedlinecount="10" %}
```cpp
int zspace_stairs_smooth_after(int h, int divisions)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: divisions, h. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/wasm/bridge/zspace_toolsets_wasm_bridge.cpp#L190)

#### zspace\_stairs\_result

{% code collapsedlinecount="10" %}
```cpp
const char * zspace_stairs_result(int h)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: h. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/wasm/bridge/zspace_toolsets_wasm_bridge.cpp#L196)

#### zspace\_toolsets\_destroy

{% code collapsedlinecount="10" %}
```cpp
int zspace_toolsets_destroy(int handle)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: handle. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/wasm/bridge/zspace_toolsets_wasm_bridge.cpp#L216)

#### zspace\_ctubes\_compute

{% code collapsedlinecount="10" %}
```cpp
int zspace_ctubes_compute(int handle, const char *path)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: handle, path. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/wasm/bridge/zspace_toolsets_wasm_bridge.cpp#L219)

#### zspace\_ctubes\_unroll

{% code collapsedlinecount="10" %}
```cpp
int zspace_ctubes_unroll(int handle, double gap)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: gap, handle. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/wasm/bridge/zspace_toolsets_wasm_bridge.cpp#L277)

#### zspace\_ctubes\_branch\_optimise

{% code collapsedlinecount="10" %}
```cpp
int zspace_ctubes_branch_optimise(int handle, const char *path)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: handle, path. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/wasm/bridge/zspace_toolsets_wasm_bridge.cpp#L279)

#### zspace\_ctubes\_branch\_result

{% code collapsedlinecount="10" %}
```cpp
const char * zspace_ctubes_branch_result(int handle)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: handle. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/wasm/bridge/zspace_toolsets_wasm_bridge.cpp#L289)

#### zspace\_ctubes\_diagnostics

{% code collapsedlinecount="10" %}
```cpp
const char * zspace_ctubes_diagnostics(int handle)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: handle. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/wasm/bridge/zspace_toolsets_wasm_bridge.cpp#L290)

#### zspace\_ctubes\_strip\_count

{% code collapsedlinecount="10" %}
```cpp
int zspace_ctubes_strip_count(int handle)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: handle. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/wasm/bridge/zspace_toolsets_wasm_bridge.cpp#L310)

#### zspace\_ctubes\_select\_mesh

{% code collapsedlinecount="10" %}
```cpp
int zspace_ctubes_select_mesh(int handle, int kind, int index)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: handle, index, kind. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/wasm/bridge/zspace_toolsets_wasm_bridge.cpp#L315)

#### zspace\_toolsets\_positions

{% code collapsedlinecount="10" %}
```cpp
const float * zspace_toolsets_positions(int h)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: h. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/wasm/bridge/zspace_toolsets_wasm_bridge.cpp#L329)

#### zspace\_toolsets\_positions\_count

{% code collapsedlinecount="10" %}
```cpp
int zspace_toolsets_positions_count(int h)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: h. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/wasm/bridge/zspace_toolsets_wasm_bridge.cpp#L330)

#### zspace\_toolsets\_face\_counts

{% code collapsedlinecount="10" %}
```cpp
const int * zspace_toolsets_face_counts(int h)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: h. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/wasm/bridge/zspace_toolsets_wasm_bridge.cpp#L331)

#### zspace\_toolsets\_face\_count

{% code collapsedlinecount="10" %}
```cpp
int zspace_toolsets_face_count(int h)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: h. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/wasm/bridge/zspace_toolsets_wasm_bridge.cpp#L332)

#### zspace\_toolsets\_face\_connects

{% code collapsedlinecount="10" %}
```cpp
const int * zspace_toolsets_face_connects(int h)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: h. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/wasm/bridge/zspace_toolsets_wasm_bridge.cpp#L333)

#### zspace\_toolsets\_face\_connects\_count

{% code collapsedlinecount="10" %}
```cpp
int zspace_toolsets_face_connects_count(int h)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: h. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/wasm/bridge/zspace_toolsets_wasm_bridge.cpp#L334)

#### zspace\_vault\_create

{% code collapsedlinecount="10" %}
```cpp
int zspace_vault_create()
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/wasm/bridge/zspace_toolsets_wasm_bridge.cpp#L336)

#### zspace\_vault\_fit

{% code collapsedlinecount="10" %}
```cpp
int zspace_vault_fit(int h, const char *path)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: h, path. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/wasm/bridge/zspace_toolsets_wasm_bridge.cpp#L339)

#### zspace\_vault\_step

{% code collapsedlinecount="10" %}
```cpp
int zspace_vault_step(int h, int iterations)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: h, iterations. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/wasm/bridge/zspace_toolsets_wasm_bridge.cpp#L387)

#### zspace\_vault\_reset

{% code collapsedlinecount="10" %}
```cpp
int zspace_vault_reset(int h)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: h. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/wasm/bridge/zspace_toolsets_wasm_bridge.cpp#L388)

#### zspace\_vault\_previous

{% code collapsedlinecount="10" %}
```cpp
int zspace_vault_previous(int h)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: h. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/wasm/bridge/zspace_toolsets_wasm_bridge.cpp#L389)

#### zspace\_vault\_diagrams

{% code collapsedlinecount="10" %}
```cpp
const char * zspace_vault_diagrams(int h, double tolerance)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: h, tolerance. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/wasm/bridge/zspace_toolsets_wasm_bridge.cpp#L390)

#### zspace\_vault\_select\_mesh

{% code collapsedlinecount="10" %}
```cpp
int zspace_vault_select_mesh(int h)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: h. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/wasm/bridge/zspace_toolsets_wasm_bridge.cpp#L419)

#### zspace\_vault\_result

{% code collapsedlinecount="10" %}
```cpp
const char * zspace_vault_result(int h)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: h. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/wasm/bridge/zspace_toolsets_wasm_bridge.cpp#L423)

#### zspace\_toolsets\_result\_buffer\_size

{% code collapsedlinecount="10" %}
```cpp
int zspace_toolsets_result_buffer_size(int h)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: h. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/wasm/bridge/zspace_toolsets_wasm_bridge.cpp#L437)

#### zspace\_stairs\_update\_from\_buffers

{% code collapsedlinecount="10" %}
```cpp
int zspace_stairs_update_from_buffers(int h, const void *bytes, size_t size)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: bytes, h, size. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/wasm/bridge/zspace_toolsets_wasm_bridge.cpp#L438)

#### zspace\_stairs\_result\_buffers

{% code collapsedlinecount="10" %}
```cpp
const void * zspace_stairs_result_buffers(int h)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: h. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/wasm/bridge/zspace_toolsets_wasm_bridge.cpp#L439)

#### zspace\_ctubes\_compute\_from\_buffers

{% code collapsedlinecount="10" %}
```cpp
int zspace_ctubes_compute_from_buffers(int h, const void *bytes, size_t size)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: bytes, h, size. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/wasm/bridge/zspace_toolsets_wasm_bridge.cpp#L440)

#### zspace\_ctubes\_branch\_optimise\_from\_buffers

{% code collapsedlinecount="10" %}
```cpp
int zspace_ctubes_branch_optimise_from_buffers(int h, const void *bytes, size_t size)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: bytes, h, size. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/wasm/bridge/zspace_toolsets_wasm_bridge.cpp#L441)

#### zspace\_vault\_fit\_from\_buffers

{% code collapsedlinecount="10" %}
```cpp
int zspace_vault_fit_from_buffers(int h, const void *bytes, size_t size)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: bytes, h, size. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/wasm/bridge/zspace_toolsets_wasm_bridge.cpp#L442)

#### zspace\_panelisation\_compute\_from\_buffers

{% code collapsedlinecount="10" %}
```cpp
int zspace_panelisation_compute_from_buffers(int h, const void *bytes, size_t size)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: bytes, h, size. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/wasm/bridge/zspace_toolsets_wasm_bridge.cpp#L443)

#### zspace\_ctubes\_branch\_result\_buffers

{% code collapsedlinecount="10" %}
```cpp
const void * zspace_ctubes_branch_result_buffers(int h)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: h. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/wasm/bridge/zspace_toolsets_wasm_bridge.cpp#L444)

#### zspace\_ctubes\_diagnostics\_buffers

{% code collapsedlinecount="10" %}
```cpp
const void * zspace_ctubes_diagnostics_buffers(int h)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: h. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/wasm/bridge/zspace_toolsets_wasm_bridge.cpp#L445)

#### zspace\_vault\_result\_buffers

{% code collapsedlinecount="10" %}
```cpp
const void * zspace_vault_result_buffers(int h)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: h. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/wasm/bridge/zspace_toolsets_wasm_bridge.cpp#L446)

#### zspace\_panelisation\_result\_buffers

{% code collapsedlinecount="10" %}
```cpp
const void * zspace_panelisation_result_buffers(int h)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: h. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/wasm/bridge/zspace_toolsets_wasm_bridge.cpp#L447)

#### zspace\_vault\_diagrams\_buffers

{% code collapsedlinecount="10" %}
```cpp
const void * zspace_vault_diagrams_buffers(int h, double tolerance)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: h, tolerance. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/wasm/bridge/zspace_toolsets_wasm_bridge.cpp#L448)

#### zspace\_stairs\_optimise\_graph\_from\_buffers

{% code collapsedlinecount="10" %}
```cpp
int zspace_stairs_optimise_graph_from_buffers(int h, const void *bytes, size_t size)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: bytes, h, size. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/wasm/bridge/zspace_toolsets_wasm_bridge.cpp#L452)

#### zspace\_stairs\_optimised\_graph\_buffers

{% code collapsedlinecount="10" %}
```cpp
const void * zspace_stairs_optimised_graph_buffers(int h)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: h. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/wasm/bridge/zspace_toolsets_wasm_bridge.cpp#L453)
