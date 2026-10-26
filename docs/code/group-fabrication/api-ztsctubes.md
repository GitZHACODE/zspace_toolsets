# zTsCTubes

## zSpace::zTsCTubes

Declared in include/zspace/zToolsets/fabrication/zTsCTubes.h. Public class/struct declaration reference.

[Header source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/zTsCTubes.h#L42)

### Public Func

#### setDirectrix

{% code collapsedlinecount="10" %}
```cpp
void zSpace::zTsCTubes::setDirectrix(zObjectGraph &graph)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: graph. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/zTsCTubes.h#L45)

#### setGeneratrix

{% code collapsedlinecount="10" %}
```cpp
void zSpace::zTsCTubes::setGeneratrix(zObjectGraph &graph)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: graph. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/zTsCTubes.h#L46)

#### setGeneratrices

{% code collapsedlinecount="10" %}
```cpp
void zSpace::zTsCTubes::setGeneratrices(zObjectGraphArray &graphs)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: graphs. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/zTsCTubes.h#L47)

#### setGeneratrixOrigin

{% code collapsedlinecount="10" %}
```cpp
void zSpace::zTsCTubes::setGeneratrixOrigin(int profileIndex, const zPoint &vertex)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: profileIndex, vertex. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/zTsCTubes.h#L48)

#### fitReport

{% code collapsedlinecount="10" %}
```cpp
const zCTubeFitReport & zSpace::zTsCTubes::fitReport() const
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/zTsCTubes.h#L49)

#### setSettings

{% code collapsedlinecount="10" %}
```cpp
void zSpace::zTsCTubes::setSettings(const zCTubeSettings &settings)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: settings. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/zTsCTubes.h#L50)

#### computeTube

{% code collapsedlinecount="10" %}
```cpp
void zSpace::zTsCTubes::computeTube(zCTubeMode mode=zCTubeMode::Generic)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: mode. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/zTsCTubes.h#L51)

#### computeTubeWithFrames

{% code collapsedlinecount="10" %}
```cpp
void zSpace::zTsCTubes::computeTubeWithFrames(const std::vector< zCTubeFrame > &frames, const zPointArray &flangeTargets={}, int flangeVertex=-1)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: flangeTargets, flangeVertex, frames. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/zTsCTubes.h#L54)

#### computeTubeWithFrames

{% code collapsedlinecount="10" %}
```cpp
void zSpace::zTsCTubes::computeTubeWithFrames(const std::vector< zCTubeFrame > &frames, const zPointArray &flangeTargets, int flangeVertex, const zDoubleArray &apexScales)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: apexScales, flangeTargets, flangeVertex, frames. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/zTsCTubes.h#L56)

#### computeTubeWithSections

{% code collapsedlinecount="10" %}
```cpp
void zSpace::zTsCTubes::computeTubeWithSections(const std::vector< zCTubeFrame > &frames, const zPointArray &sections)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: frames, sections. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/zTsCTubes.h#L59)

#### computeStrips

{% code collapsedlinecount="10" %}
```cpp
void zSpace::zTsCTubes::computeStrips()
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/zTsCTubes.h#L60)

#### unrollStrips

{% code collapsedlinecount="10" %}
```cpp
void zSpace::zTsCTubes::unrollStrips(double gap=0.1)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: gap. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/zTsCTubes.h#L61)

#### clearResults

{% code collapsedlinecount="10" %}
```cpp
void zSpace::zTsCTubes::clearResults()
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/zTsCTubes.h#L62)

#### frames

{% code collapsedlinecount="10" %}
```cpp
const std::vector< zCTubeFrame > & zSpace::zTsCTubes::frames() const
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/zTsCTubes.h#L63)

#### stripAreas

{% code collapsedlinecount="10" %}
```cpp
const std::vector< zCTubeStripArea > & zSpace::zTsCTubes::stripAreas() const
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/zTsCTubes.h#L64)

#### tubeMesh

{% code collapsedlinecount="10" %}
```cpp
zObjectMesh & zSpace::zTsCTubes::tubeMesh()
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/zTsCTubes.h#L65)

#### sampledDirectrix

{% code collapsedlinecount="10" %}
```cpp
zObjectGraph & zSpace::zTsCTubes::sampledDirectrix()
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/zTsCTubes.h#L66)

#### stripMeshes

{% code collapsedlinecount="10" %}
```cpp
zObjectMeshArray & zSpace::zTsCTubes::stripMeshes()
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/zTsCTubes.h#L67)

#### unrolledMeshes

{% code collapsedlinecount="10" %}
```cpp
zObjectMeshArray & zSpace::zTsCTubes::unrolledMeshes()
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/zTsCTubes.h#L68)

### Related types

[zCTubeFitReport](api-zctubefitreport.md)

[zCTubeFrame](api-zctubeframe.md)

[zCTubeSettings](api-zctubesettings.md)

[zCTubeStripArea](api-zctubestriparea.md)
