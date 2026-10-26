# z3DPSlicer

## zSpace::z3DPSlicer

Declared in include/zspace/zToolsets/fabrication/z3DPSlicer.h. Public class/struct declaration reference.

[Header source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/z3DPSlicer.h#L14)

### Public Func

#### setInputMesh

{% code collapsedlinecount="10" %}
```cpp
void zSpace::z3DPSlicer::setInputMesh(zObjectMesh &mesh)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: mesh. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/z3DPSlicer.h#L17)

#### setSettings

{% code collapsedlinecount="10" %}
```cpp
void zSpace::z3DPSlicer::setSettings(const z3DPSettings &settings)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: settings. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/z3DPSlicer.h#L18)

#### identifySlicingFeatures

{% code collapsedlinecount="10" %}
```cpp
void zSpace::z3DPSlicer::identifySlicingFeatures()
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/z3DPSlicer.h#L20)

#### computeSliceCount

{% code collapsedlinecount="10" %}
```cpp
void zSpace::z3DPSlicer::computeSliceCount()
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/z3DPSlicer.h#L21)

#### computeSliceMeshes

{% code collapsedlinecount="10" %}
```cpp
void zSpace::z3DPSlicer::computeSliceMeshes()
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/z3DPSlicer.h#L22)

#### triangulateSlices

{% code collapsedlinecount="10" %}
```cpp
void zSpace::z3DPSlicer::triangulateSlices()
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/z3DPSlicer.h#L23)

#### unrollSlices

{% code collapsedlinecount="10" %}
```cpp
void zSpace::z3DPSlicer::unrollSlices()
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/z3DPSlicer.h#L24)

#### computeSDFFields

{% code collapsedlinecount="10" %}
```cpp
void zSpace::z3DPSlicer::computeSDFFields()
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/z3DPSlicer.h#L25)

#### extractSDFContours

{% code collapsedlinecount="10" %}
```cpp
void zSpace::z3DPSlicer::extractSDFContours()
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/z3DPSlicer.h#L26)

#### sliceMeshes

{% code collapsedlinecount="10" %}
```cpp
zObjectMeshArray & zSpace::z3DPSlicer::sliceMeshes()
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/z3DPSlicer.h#L28)

#### triangulatedMeshes

{% code collapsedlinecount="10" %}
```cpp
zObjectMeshArray & zSpace::z3DPSlicer::triangulatedMeshes()
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/z3DPSlicer.h#L29)

#### unrolledMeshes

{% code collapsedlinecount="10" %}
```cpp
zObjectMeshArray & zSpace::z3DPSlicer::unrolledMeshes()
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/z3DPSlicer.h#L30)

#### fields

{% code collapsedlinecount="10" %}
```cpp
zObjectMeshScalarFieldArray & zSpace::z3DPSlicer::fields()
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/z3DPSlicer.h#L31)

#### contours

{% code collapsedlinecount="10" %}
```cpp
zObjectGraphArray & zSpace::z3DPSlicer::contours()
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/z3DPSlicer.h#L32)

#### edgeLoops

{% code collapsedlinecount="10" %}
```cpp
zObjectGraphArray & zSpace::z3DPSlicer::edgeLoops()
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/z3DPSlicer.h#L33)

#### cornerPoints

{% code collapsedlinecount="10" %}
```cpp
zObjectPointCloud & zSpace::z3DPSlicer::cornerPoints()
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/z3DPSlicer.h#L34)

#### topMesh

{% code collapsedlinecount="10" %}
```cpp
zObjectMesh & zSpace::z3DPSlicer::topMesh()
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/z3DPSlicer.h#L35)

#### bottomMesh

{% code collapsedlinecount="10" %}
```cpp
zObjectMesh & zSpace::z3DPSlicer::bottomMesh()
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/z3DPSlicer.h#L36)

#### sliceMeshes

{% code collapsedlinecount="10" %}
```cpp
const zObjectMeshArray & zSpace::z3DPSlicer::sliceMeshes() const
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/z3DPSlicer.h#L38)

#### unrolledMeshes

{% code collapsedlinecount="10" %}
```cpp
const zObjectMeshArray & zSpace::z3DPSlicer::unrolledMeshes() const
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/z3DPSlicer.h#L39)

#### contours

{% code collapsedlinecount="10" %}
```cpp
const zObjectGraphArray & zSpace::z3DPSlicer::contours() const
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/z3DPSlicer.h#L40)

#### features

{% code collapsedlinecount="10" %}
```cpp
const z3DPSlicingFeatures & zSpace::z3DPSlicer::features() const
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/z3DPSlicer.h#L41)

### Related types

[z3DPSettings](api-z3dpsettings.md)

[z3DPSlicingFeatures](api-z3dpslicingfeatures.md)
