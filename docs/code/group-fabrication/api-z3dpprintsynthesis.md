# z3DPPrintSynthesis

## zSpace::z3DPPrintSynthesis

Declared in include/zspace/zToolsets/fabrication/z3DPPrintSynthesis.h. Public class/struct declaration reference.

[Header source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/z3DPPrintSynthesis.h#L12)

### Public Func

#### setSettings

{% code collapsedlinecount="10" %}
```cpp
void zSpace::z3DPPrintSynthesis::setSettings(const z3DPSettings &settings)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: settings. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/z3DPPrintSynthesis.h#L15)

#### setInputs

{% code collapsedlinecount="10" %}
```cpp
void zSpace::z3DPPrintSynthesis::setInputs(zObjectGraphArray &contours, zObjectMeshArray &unrolledMeshes, zObjectMeshArray &sliceMeshes)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: contours, sliceMeshes, unrolledMeshes. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/z3DPPrintSynthesis.h#L16)

#### resampleContours

{% code collapsedlinecount="10" %}
```cpp
void zSpace::z3DPPrintSynthesis::resampleContours()
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/z3DPPrintSynthesis.h#L18)

#### identifyFeaturePoints

{% code collapsedlinecount="10" %}
```cpp
void zSpace::z3DPPrintSynthesis::identifyFeaturePoints()
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/z3DPPrintSynthesis.h#L19)

#### mapContoursToSlices

{% code collapsedlinecount="10" %}
```cpp
void zSpace::z3DPPrintSynthesis::mapContoursToSlices()
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/z3DPPrintSynthesis.h#L20)

#### sequenceToolpaths

{% code collapsedlinecount="10" %}
```cpp
void zSpace::z3DPPrintSynthesis::sequenceToolpaths()
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/z3DPPrintSynthesis.h#L21)

#### buildPrintMesh

{% code collapsedlinecount="10" %}
```cpp
void zSpace::z3DPPrintSynthesis::buildPrintMesh()
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/z3DPPrintSynthesis.h#L22)

#### printPaths

{% code collapsedlinecount="10" %}
```cpp
zObjectGraphArray & zSpace::z3DPPrintSynthesis::printPaths()
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/z3DPPrintSynthesis.h#L24)

#### featurePoints

{% code collapsedlinecount="10" %}
```cpp
zObjectPointCloudArray & zSpace::z3DPPrintSynthesis::featurePoints()
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/z3DPPrintSynthesis.h#L25)

#### printMesh

{% code collapsedlinecount="10" %}
```cpp
zObjectMesh & zSpace::z3DPPrintSynthesis::printMesh()
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/z3DPPrintSynthesis.h#L26)

#### printPaths

{% code collapsedlinecount="10" %}
```cpp
const zObjectGraphArray & zSpace::z3DPPrintSynthesis::printPaths() const
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/z3DPPrintSynthesis.h#L28)

#### printMesh

{% code collapsedlinecount="10" %}
```cpp
const zObjectMesh & zSpace::z3DPPrintSynthesis::printMesh() const
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/z3DPPrintSynthesis.h#L29)

### Related types

[z3DPSettings](api-z3dpsettings.md)
