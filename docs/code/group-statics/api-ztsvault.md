# zTsVault

## zSpace::zTsVault

Declared in include/zspace/zToolsets/statics/zTsVault.h. Public class/struct declaration reference.

[Header source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsVault.h#L42)

### Public Func

#### zTsVault

{% code collapsedlinecount="10" %}
```cpp
zSpace::zTsVault::zTsVault()
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsVault.h#L44)

#### \~zTsVault

{% code collapsedlinecount="10" %}
```cpp
zSpace::zTsVault::~zTsVault()
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsVault.h#L45)

#### zTsVault

{% code collapsedlinecount="10" %}
```cpp
zSpace::zTsVault::zTsVault(const zTsVault &)=delete
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsVault.h#L46)

#### operator=

{% code collapsedlinecount="10" %}
```cpp
zTsVault & zSpace::zTsVault::operator=(const zTsVault &)=delete
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsVault.h#L47)

#### setInputMesh

{% code collapsedlinecount="10" %}
```cpp
void zSpace::zTsVault::setInputMesh(zObjectMesh &mesh)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: mesh. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsVault.h#L48)

#### setReferencePositions

{% code collapsedlinecount="10" %}
```cpp
void zSpace::zTsVault::setReferencePositions(const zPointArray &positions)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: positions. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsVault.h#L50)

#### setGeometryValidator

{% code collapsedlinecount="10" %}
```cpp
void zSpace::zTsVault::setGeometryValidator(std::function< bool(const zPointArray &)> validator)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: validator. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsVault.h#L52)

#### setSupports

{% code collapsedlinecount="10" %}
```cpp
void zSpace::zTsVault::setSupports(const zIntArray &vertexIds)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: vertexIds. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsVault.h#L53)

#### setTensionEdges

{% code collapsedlinecount="10" %}
```cpp
void zSpace::zTsVault::setTensionEdges(const zIntArray &edgeIds)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: edgeIds. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsVault.h#L54)

#### setVertexMass

{% code collapsedlinecount="10" %}
```cpp
void zSpace::zTsVault::setVertexMass(double mass)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: mass. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsVault.h#L55)

#### setVertexMasses

{% code collapsedlinecount="10" %}
```cpp
void zSpace::zTsVault::setVertexMasses(const zDoubleArray &masses)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: masses. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsVault.h#L56)

#### setGravity

{% code collapsedlinecount="10" %}
```cpp
void zSpace::zTsVault::setGravity(const zVector &gravity)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: gravity. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsVault.h#L57)

#### setVertexLoads

{% code collapsedlinecount="10" %}
```cpp
void zSpace::zTsVault::setVertexLoads(const zVectorArray &loads)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: loads. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsVault.h#L58)

#### setForceDensityBounds

{% code collapsedlinecount="10" %}
```cpp
void zSpace::zTsVault::setForceDensityBounds(const zDoubleArray &lower, const zDoubleArray &upper)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: lower, upper. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsVault.h#L61)

#### setBestFitSettings

{% code collapsedlinecount="10" %}
```cpp
void zSpace::zTsVault::setBestFitSettings(const zVaultBestFitSettings &settings)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: settings. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsVault.h#L62)

#### computeBestFitForceDensities

{% code collapsedlinecount="10" %}
```cpp
zVaultSolveReport zSpace::zTsVault::computeBestFitForceDensities()
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsVault.h#L63)

#### stepBestFit

{% code collapsedlinecount="10" %}
```cpp
zVaultSolveReport zSpace::zTsVault::stepBestFit(int iterations=1)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: iterations. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsVault.h#L64)

#### solveBestFit

{% code collapsedlinecount="10" %}
```cpp
zVaultSolveReport zSpace::zTsVault::solveBestFit(int maxIterations)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: maxIterations. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsVault.h#L65)

#### snapFixedDensities

{% code collapsedlinecount="10" %}
```cpp
zVaultSolveReport zSpace::zTsVault::snapFixedDensities(double geometryWeight=0.0)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: geometryWeight. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsVault.h#L68)

#### resetBestFit

{% code collapsedlinecount="10" %}
```cpp
void zSpace::zTsVault::resetBestFit()
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsVault.h#L69)

#### previousBestFit

{% code collapsedlinecount="10" %}
```cpp
bool zSpace::zTsVault::previousBestFit()
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsVault.h#L70)

#### clear

{% code collapsedlinecount="10" %}
```cpp
void zSpace::zTsVault::clear()
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsVault.h#L71)

#### resultMesh

{% code collapsedlinecount="10" %}
```cpp
zObjectMesh & zSpace::zTsVault::resultMesh()
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsVault.h#L72)

#### forceDensities

{% code collapsedlinecount="10" %}
```cpp
const zDoubleArray & zSpace::zTsVault::forceDensities() const
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsVault.h#L73)

#### axialForces

{% code collapsedlinecount="10" %}
```cpp
const zDoubleArray & zSpace::zTsVault::axialForces() const
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsVault.h#L74)

#### residualForces

{% code collapsedlinecount="10" %}
```cpp
const zVectorArray & zSpace::zTsVault::residualForces() const
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsVault.h#L75)

#### supportReactions

{% code collapsedlinecount="10" %}
```cpp
const zVectorArray & zSpace::zTsVault::supportReactions() const
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsVault.h#L76)

#### activeEdgeIds

{% code collapsedlinecount="10" %}
```cpp
const zIntArray & zSpace::zTsVault::activeEdgeIds() const
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsVault.h#L77)

#### edgeConnects

{% code collapsedlinecount="10" %}
```cpp
const zIntArray & zSpace::zTsVault::edgeConnects() const
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsVault.h#L78)

#### report

{% code collapsedlinecount="10" %}
```cpp
const zVaultSolveReport & zSpace::zTsVault::report() const
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsVault.h#L79)

#### createFormFromResult

{% code collapsedlinecount="10" %}
```cpp
zVaultDiagramReport zSpace::zTsVault::createFormFromResult()
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsVault.h#L81)

#### createForceFromForm

{% code collapsedlinecount="10" %}
```cpp
zVaultDiagramReport zSpace::zTsVault::createForceFromForm(const zVaultDiagramSettings &settings={})
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: settings. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsVault.h#L82)

#### updateForceFromDensities

{% code collapsedlinecount="10" %}
```cpp
zVaultDiagramReport zSpace::zTsVault::updateForceFromDensities()
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsVault.h#L83)

#### formMesh

{% code collapsedlinecount="10" %}
```cpp
const zObjectMesh & zSpace::zTsVault::formMesh() const
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsVault.h#L84)

#### forceMesh

{% code collapsedlinecount="10" %}
```cpp
const zObjectMesh & zSpace::zTsVault::forceMesh() const
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsVault.h#L85)

#### formHalfEdgeToForceHalfEdge

{% code collapsedlinecount="10" %}
```cpp
const zIntArray & zSpace::zTsVault::formHalfEdgeToForceHalfEdge() const
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsVault.h#L86)

#### forceHalfEdgeToFormHalfEdge

{% code collapsedlinecount="10" %}
```cpp
const zIntArray & zSpace::zTsVault::forceHalfEdgeToFormHalfEdge() const
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsVault.h#L87)

#### diagramReport

{% code collapsedlinecount="10" %}
```cpp
const zVaultDiagramReport & zSpace::zTsVault::diagramReport() const
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/statics/zTsVault.h#L88)

### Related types

[zVaultBestFitSettings](api-zvaultbestfitsettings.md)

[zVaultDiagramReport](api-zvaultdiagramreport.md)

[zVaultDiagramSettings](api-zvaultdiagramsettings.md)

[zVaultSolveReport](api-zvaultsolvereport.md)
