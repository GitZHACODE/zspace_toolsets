# zTsForcePolyhedra

## zSpace::zTsForcePolyhedra

Declared in include/zspace/zToolsets/statics/zTsForcePolyhedra.h. Public class/struct declaration reference.

[Header source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsForcePolyhedra.h#L65)

### Public Func

#### zTsForcePolyhedra

{% code collapsedlinecount="10" %}
```cpp
zSpace::zTsForcePolyhedra::zTsForcePolyhedra()
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsForcePolyhedra.h#L67)

#### \~zTsForcePolyhedra

{% code collapsedlinecount="10" %}
```cpp
zSpace::zTsForcePolyhedra::~zTsForcePolyhedra()
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsForcePolyhedra.h#L68)

#### zTsForcePolyhedra

{% code collapsedlinecount="10" %}
```cpp
zSpace::zTsForcePolyhedra::zTsForcePolyhedra(const zTsForcePolyhedra &)=delete
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsForcePolyhedra.h#L69)

#### operator=

{% code collapsedlinecount="10" %}
```cpp
zTsForcePolyhedra & zSpace::zTsForcePolyhedra::operator=(const zTsForcePolyhedra &)=delete
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsForcePolyhedra.h#L70)

#### clear

{% code collapsedlinecount="10" %}
```cpp
void zSpace::zTsForcePolyhedra::clear()
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsForcePolyhedra.h#L71)

#### setInputGraph

{% code collapsedlinecount="10" %}
```cpp
void zSpace::zTsForcePolyhedra::setInputGraph(zObjectGraph &graph)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: graph. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsForcePolyhedra.h#L72)

#### setFormPositions

{% code collapsedlinecount="10" %}
```cpp
void zSpace::zTsForcePolyhedra::setFormPositions(const zPointArray &positions)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: positions. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsForcePolyhedra.h#L75)

#### setSupports

{% code collapsedlinecount="10" %}
```cpp
void zSpace::zTsForcePolyhedra::setSupports(const zIntArray &vertexIds)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: vertexIds. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsForcePolyhedra.h#L76)

#### setTensionEdges

{% code collapsedlinecount="10" %}
```cpp
void zSpace::zTsForcePolyhedra::setTensionEdges(const zIntArray &edgeIds)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: edgeIds. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsForcePolyhedra.h#L77)

#### setEdgeForceMagnitudes

{% code collapsedlinecount="10" %}
```cpp
void zSpace::zTsForcePolyhedra::setEdgeForceMagnitudes(const zDoubleArray &magnitudes)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: magnitudes. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsForcePolyhedra.h#L80)

#### setVertexMass

{% code collapsedlinecount="10" %}
```cpp
void zSpace::zTsForcePolyhedra::setVertexMass(double mass)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: mass. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsForcePolyhedra.h#L81)

#### setVertexMasses

{% code collapsedlinecount="10" %}
```cpp
void zSpace::zTsForcePolyhedra::setVertexMasses(const zDoubleArray &masses)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: masses. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsForcePolyhedra.h#L82)

#### setGravity

{% code collapsedlinecount="10" %}
```cpp
void zSpace::zTsForcePolyhedra::setGravity(const zVector &gravity)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: gravity. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsForcePolyhedra.h#L83)

#### setVertexLoads

{% code collapsedlinecount="10" %}
```cpp
void zSpace::zTsForcePolyhedra::setVertexLoads(const zVectorArray &loads)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: loads. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsForcePolyhedra.h#L84)

#### setSettings

{% code collapsedlinecount="10" %}
```cpp
void zSpace::zTsForcePolyhedra::setSettings(const zForcePolyhedraSettings &settings)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: settings. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsForcePolyhedra.h#L85)

#### setFormWeight

{% code collapsedlinecount="10" %}
```cpp
void zSpace::zTsForcePolyhedra::setFormWeight(double weight)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: weight. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsForcePolyhedra.h#L86)

#### createForcePolyhedra

{% code collapsedlinecount="10" %}
```cpp
zForcePolyhedraReport zSpace::zTsForcePolyhedra::createForcePolyhedra()
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsForcePolyhedra.h#L90)

#### stepAreaPursuit

{% code collapsedlinecount="10" %}
```cpp
zForcePolyhedraReport zSpace::zTsForcePolyhedra::stepAreaPursuit(int iterations=1)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: iterations. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsForcePolyhedra.h#L91)

#### stepEquilibrium

{% code collapsedlinecount="10" %}
```cpp
zForcePolyhedraReport zSpace::zTsForcePolyhedra::stepEquilibrium(int iterations=1)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: iterations. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsForcePolyhedra.h#L95)

#### solveEquilibrium

{% code collapsedlinecount="10" %}
```cpp
zForcePolyhedraReport zSpace::zTsForcePolyhedra::solveEquilibrium(int maxIterations=2000)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: maxIterations. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsForcePolyhedra.h#L96)

#### reset

{% code collapsedlinecount="10" %}
```cpp
void zSpace::zTsForcePolyhedra::reset()
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsForcePolyhedra.h#L97)

#### resultGraph

{% code collapsedlinecount="10" %}
```cpp
const zObjectGraph & zSpace::zTsForcePolyhedra::resultGraph() const
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsForcePolyhedra.h#L98)

#### edgeConnects

{% code collapsedlinecount="10" %}
```cpp
const zIntArray & zSpace::zTsForcePolyhedra::edgeConnects() const
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsForcePolyhedra.h#L99)

#### axialForces

{% code collapsedlinecount="10" %}
```cpp
const zDoubleArray & zSpace::zTsForcePolyhedra::axialForces() const
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsForcePolyhedra.h#L100)

#### residualForces

{% code collapsedlinecount="10" %}
```cpp
const zVectorArray & zSpace::zTsForcePolyhedra::residualForces() const
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsForcePolyhedra.h#L101)

#### supportReactions

{% code collapsedlinecount="10" %}
```cpp
const zVectorArray & zSpace::zTsForcePolyhedra::supportReactions() const
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsForcePolyhedra.h#L102)

#### forceCell

{% code collapsedlinecount="10" %}
```cpp
const zObjectMesh & zSpace::zTsForcePolyhedra::forceCell(int vertexId) const
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: vertexId. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsForcePolyhedra.h#L105)

#### cellFaces

{% code collapsedlinecount="10" %}
```cpp
const std::vector< zForcePolyhedraFace > & zSpace::zTsForcePolyhedra::cellFaces(int vertexId) const
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: vertexId. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsForcePolyhedra.h#L106)

#### cellReport

{% code collapsedlinecount="10" %}
```cpp
const zForcePolyhedraCellReport & zSpace::zTsForcePolyhedra::cellReport(int vertexId) const
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: vertexId. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsForcePolyhedra.h#L107)

#### report

{% code collapsedlinecount="10" %}
```cpp
const zForcePolyhedraReport & zSpace::zTsForcePolyhedra::report() const
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsForcePolyhedra.h#L108)

### Related types

[zForcePolyhedraCellReport](api-zforcepolyhedracellreport.md)

[zForcePolyhedraFace](api-zforcepolyhedraface.md)

[zForcePolyhedraReport](api-zforcepolyhedrareport.md)

[zForcePolyhedraSettings](api-zforcepolyhedrasettings.md)
