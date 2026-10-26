# zTsForceCTubes

## zSpace::zTsForceCTubes

Declared in include/zspace/zToolsets/statics/zTsForceCTubes.h. Public class/struct declaration reference.

[Header source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsForceCTubes.h#L37)

### Public Func

#### zTsForceCTubes

{% code collapsedlinecount="10" %}
```cpp
zSpace::zTsForceCTubes::zTsForceCTubes()
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsForceCTubes.h#L39)

#### \~zTsForceCTubes

{% code collapsedlinecount="10" %}
```cpp
zSpace::zTsForceCTubes::~zTsForceCTubes()
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsForceCTubes.h#L40)

#### zTsForceCTubes

{% code collapsedlinecount="10" %}
```cpp
zSpace::zTsForceCTubes::zTsForceCTubes(const zTsForceCTubes &)=delete
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsForceCTubes.h#L41)

#### operator=

{% code collapsedlinecount="10" %}
```cpp
zTsForceCTubes & zSpace::zTsForceCTubes::operator=(const zTsForceCTubes &)=delete
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsForceCTubes.h#L42)

#### setInput

{% code collapsedlinecount="10" %}
```cpp
void zSpace::zTsForceCTubes::setInput(zObjectGraph &directrix, zObjectGraph &generatrix)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: directrix, generatrix. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsForceCTubes.h#L43)

#### setForceDiagram

{% code collapsedlinecount="10" %}
```cpp
void zSpace::zTsForceCTubes::setForceDiagram(const zTsForcePolyhedra &force, double forceScale=1)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: force, forceScale. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsForceCTubes.h#L44)

#### setForceData

{% code collapsedlinecount="10" %}
```cpp
void zSpace::zTsForceCTubes::setForceData(zObjectGraph &form, const zDoubleArray &memberFaceAreas, const zDoubleArray &axialForces, const zDoubleArray &pairedAreaDifferences, bool sourceEquilibrium=false)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: axialForces, form, memberFaceAreas, pairedAreaDifferences, sourceEquilibrium. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsForceCTubes.h#L45)

#### setSettings

{% code collapsedlinecount="10" %}
```cpp
void zSpace::zTsForceCTubes::setSettings(const zForceCTubeSettings &settings)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: settings. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsForceCTubes.h#L48)

#### compute

{% code collapsedlinecount="10" %}
```cpp
zForceCTubeReport zSpace::zTsForceCTubes::compute(bool optimise=false)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: optimise. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsForceCTubes.h#L49)

#### tubeMesh

{% code collapsedlinecount="10" %}
```cpp
zObjectMesh & zSpace::zTsForceCTubes::tubeMesh()
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsForceCTubes.h#L50)

#### frames

{% code collapsedlinecount="10" %}
```cpp
const std::vector< zCTubeFrame > & zSpace::zTsForceCTubes::frames() const
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsForceCTubes.h#L51)

#### sections

{% code collapsedlinecount="10" %}
```cpp
const std::vector< zForceCTubeSection > & zSpace::zTsForceCTubes::sections() const
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsForceCTubes.h#L52)

#### apexScales

{% code collapsedlinecount="10" %}
```cpp
const zDoubleArray & zSpace::zTsForceCTubes::apexScales() const
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsForceCTubes.h#L53)

#### report

{% code collapsedlinecount="10" %}
```cpp
const zForceCTubeReport & zSpace::zTsForceCTubes::report() const
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsForceCTubes.h#L54)

### Related types

[zCTubeFrame](../group-fabrication/api-zctubeframe.md)

[zForceCTubeReport](api-zforcectubereport.md)

[zForceCTubeSection](api-zforcectubesection.md)

[zForceCTubeSettings](api-zforcectubesettings.md)

[zTsForcePolyhedra](api-ztsforcepolyhedra.md)
