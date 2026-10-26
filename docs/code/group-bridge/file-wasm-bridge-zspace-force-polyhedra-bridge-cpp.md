# zspace\_force\_polyhedra\_bridge.cpp declarations

## wasm/bridge/zspace\_force\_polyhedra\_bridge.cpp

[Source file](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/wasm/bridge/zspace_force_polyhedra_bridge.cpp#L1)

Free functions, aliases, enums and constants declared in this file. Class and struct members have separate pages.

#### zspace\_forcepolyhedra\_create

{% code collapsedlinecount="10" %}
```cpp
int zspace_forcepolyhedra_create()
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/wasm/bridge/zspace_force_polyhedra_bridge.cpp#L103)

#### zspace\_forcepolyhedra\_destroy

{% code collapsedlinecount="10" %}
```cpp
int zspace_forcepolyhedra_destroy(int h)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: h. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/wasm/bridge/zspace_force_polyhedra_bridge.cpp#L104)

#### zspace\_forcepolyhedra\_last\_error

{% code collapsedlinecount="10" %}
```cpp
const char * zspace_forcepolyhedra_last_error()
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/wasm/bridge/zspace_force_polyhedra_bridge.cpp#L105)

#### zspace\_forcepolyhedra\_compute

{% code collapsedlinecount="10" %}
```cpp
int zspace_forcepolyhedra_compute(int h, const char *path)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: h, path. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/wasm/bridge/zspace_force_polyhedra_bridge.cpp#L106)

#### zspace\_forcepolyhedra\_compute\_from\_buffers

{% code collapsedlinecount="10" %}
```cpp
int zspace_forcepolyhedra_compute_from_buffers(int h, const void *bytes, size_t length)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: bytes, h, length. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/wasm/bridge/zspace_force_polyhedra_bridge.cpp#L110)

#### zspace\_forcepolyhedra\_step

{% code collapsedlinecount="10" %}
```cpp
int zspace_forcepolyhedra_step(int h, int iterations)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: h, iterations. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/wasm/bridge/zspace_force_polyhedra_bridge.cpp#L115)

#### zspace\_forcepolyhedra\_equilibrium

{% code collapsedlinecount="10" %}
```cpp
int zspace_forcepolyhedra_equilibrium(int h, int iterations)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: h, iterations. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/wasm/bridge/zspace_force_polyhedra_bridge.cpp#L116)

#### zspace\_forcepolyhedra\_set\_form\_weight

{% code collapsedlinecount="10" %}
```cpp
int zspace_forcepolyhedra_set_form_weight(int h, double weight)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: h, weight. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/wasm/bridge/zspace_force_polyhedra_bridge.cpp#L117)

#### zspace\_forcepolyhedra\_set\_controls

{% code collapsedlinecount="10" %}
```cpp
int zspace_forcepolyhedra_set_controls(int h, double areaTolerance, double angleTolerance, double weight, double normalTolerance)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: angleTolerance, areaTolerance, h, normalTolerance, weight. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/wasm/bridge/zspace_force_polyhedra_bridge.cpp#L119)

#### zspace\_forcepolyhedra\_set\_form\_from\_buffers

{% code collapsedlinecount="10" %}
```cpp
int zspace_forcepolyhedra_set_form_from_buffers(int h, const void *bytes, size_t length)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: bytes, h, length. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/wasm/bridge/zspace_force_polyhedra_bridge.cpp#L124)

#### zspace\_forcepolyhedra\_result

{% code collapsedlinecount="10" %}
```cpp
const char * zspace_forcepolyhedra_result(int h)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: h. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/wasm/bridge/zspace_force_polyhedra_bridge.cpp#L129)

#### zspace\_forcepolyhedra\_result\_buffers

{% code collapsedlinecount="10" %}
```cpp
const void * zspace_forcepolyhedra_result_buffers(int h)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: h. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/wasm/bridge/zspace_force_polyhedra_bridge.cpp#L130)

#### zspace\_forcepolyhedra\_result\_buffer\_size

{% code collapsedlinecount="10" %}
```cpp
size_t zspace_forcepolyhedra_result_buffer_size(int h)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: h. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/wasm/bridge/zspace_force_polyhedra_bridge.cpp#L131)
