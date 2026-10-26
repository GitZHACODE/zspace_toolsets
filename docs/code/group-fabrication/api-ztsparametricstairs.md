# zTsParametricStairs

## zSpace::zTsParametricStairs

Declared in include/zspace/zToolsets/fabrication/zTsParametricStairs.h. Public class/struct declaration reference.

[Header source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/zTsParametricStairs.h#L36)

### Public Static Func

#### optimiseGraph

{% code collapsedlinecount="10" %}
```cpp
static zStairsGraphOptimisation zSpace::zTsParametricStairs::optimiseGraph(zObjectGraph &graph, const zParametricStairsSettings &settings, double treadMin, double treadMax, double riserMin, double riserMax)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: graph, riserMax, riserMin, settings, treadMax, treadMin. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/zTsParametricStairs.h#L38)

### Public Func

#### update

{% code collapsedlinecount="10" %}
```cpp
void zSpace::zTsParametricStairs::update(zObjectGraph &graph, const zParametricStairsSettings &settings={})
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: graph, settings. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/zTsParametricStairs.h#L41)

#### stepBestFit

{% code collapsedlinecount="10" %}
```cpp
zVaultSolveReport zSpace::zTsParametricStairs::stepBestFit(int iterations=1)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: iterations. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/zTsParametricStairs.h#L42)

#### snapFixedDensities

{% code collapsedlinecount="10" %}
```cpp
zVaultSolveReport zSpace::zTsParametricStairs::snapFixedDensities(double geometryWeight=0.0)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: geometryWeight. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/zTsParametricStairs.h#L43)

#### smoothResult

{% code collapsedlinecount="10" %}
```cpp
void zSpace::zTsParametricStairs::smoothResult(int divisions=1)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: divisions. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/zTsParametricStairs.h#L44)

#### resetBestFit

{% code collapsedlinecount="10" %}
```cpp
void zSpace::zTsParametricStairs::resetBestFit()
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/zTsParametricStairs.h#L45)

#### previousBestFit

{% code collapsedlinecount="10" %}
```cpp
bool zSpace::zTsParametricStairs::previousBestFit()
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/zTsParametricStairs.h#L46)

#### centerMesh

{% code collapsedlinecount="10" %}
```cpp
zObjectMesh & zSpace::zTsParametricStairs::centerMesh()
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/zTsParametricStairs.h#L47)

#### lowPolyMesh

{% code collapsedlinecount="10" %}
```cpp
zObjectMesh & zSpace::zTsParametricStairs::lowPolyMesh()
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/zTsParametricStairs.h#L48)

#### subdividedMesh

{% code collapsedlinecount="10" %}
```cpp
zObjectMesh & zSpace::zTsParametricStairs::subdividedMesh()
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/zTsParametricStairs.h#L49)

#### resultMesh

{% code collapsedlinecount="10" %}
```cpp
zObjectMesh & zSpace::zTsParametricStairs::resultMesh()
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/zTsParametricStairs.h#L50)

#### hasStairLayout

{% code collapsedlinecount="10" %}
```cpp
bool zSpace::zTsParametricStairs::hasStairLayout() const
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/zTsParametricStairs.h#L51)

#### stairMesh

{% code collapsedlinecount="10" %}
```cpp
zObjectMesh & zSpace::zTsParametricStairs::stairMesh()
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/zTsParametricStairs.h#L52)

#### walkingGraph

{% code collapsedlinecount="10" %}
```cpp
zObjectGraph & zSpace::zTsParametricStairs::walkingGraph()
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/zTsParametricStairs.h#L53)

#### stairFaceTypes

{% code collapsedlinecount="10" %}
```cpp
const zIntArray & zSpace::zTsParametricStairs::stairFaceTypes() const
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/zTsParametricStairs.h#L54)

#### flightRiserCounts

{% code collapsedlinecount="10" %}
```cpp
const zIntArray & zSpace::zTsParametricStairs::flightRiserCounts() const
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/zTsParametricStairs.h#L55)

#### flightActualRisers

{% code collapsedlinecount="10" %}
```cpp
const zDoubleArray & zSpace::zTsParametricStairs::flightActualRisers() const
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/zTsParametricStairs.h#L56)

#### minimumClearance

{% code collapsedlinecount="10" %}
```cpp
double zSpace::zTsParametricStairs::minimumClearance() const
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/zTsParametricStairs.h#L57)

#### refreshLayout

{% code collapsedlinecount="10" %}
```cpp
void zSpace::zTsParametricStairs::refreshLayout()
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/zTsParametricStairs.h#L58)

#### supports

{% code collapsedlinecount="10" %}
```cpp
const zIntArray & zSpace::zTsParametricStairs::supports() const
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/zTsParametricStairs.h#L59)

#### meshVertexGraphVertices

{% code collapsedlinecount="10" %}
```cpp
const zIntArray & zSpace::zTsParametricStairs::meshVertexGraphVertices() const
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/zTsParametricStairs.h#L60)

#### centerGraphMeshVertices

{% code collapsedlinecount="10" %}
```cpp
const zIntArray & zSpace::zTsParametricStairs::centerGraphMeshVertices() const
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/zTsParametricStairs.h#L61)

#### forceDensities

{% code collapsedlinecount="10" %}
```cpp
const zDoubleArray & zSpace::zTsParametricStairs::forceDensities() const
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/zTsParametricStairs.h#L62)

#### edgeConnects

{% code collapsedlinecount="10" %}
```cpp
const zIntArray & zSpace::zTsParametricStairs::edgeConnects() const
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/zTsParametricStairs.h#L63)

#### activeEdgeIds

{% code collapsedlinecount="10" %}
```cpp
const zIntArray & zSpace::zTsParametricStairs::activeEdgeIds() const
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/zTsParametricStairs.h#L64)

#### axialForces

{% code collapsedlinecount="10" %}
```cpp
const zDoubleArray & zSpace::zTsParametricStairs::axialForces() const
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/zTsParametricStairs.h#L65)

#### residualForces

{% code collapsedlinecount="10" %}
```cpp
const zVectorArray & zSpace::zTsParametricStairs::residualForces() const
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/zTsParametricStairs.h#L66)

#### supportReactions

{% code collapsedlinecount="10" %}
```cpp
const zVectorArray & zSpace::zTsParametricStairs::supportReactions() const
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/zTsParametricStairs.h#L67)

#### report

{% code collapsedlinecount="10" %}
```cpp
const zVaultSolveReport & zSpace::zTsParametricStairs::report() const
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/zTsParametricStairs.h#L68)

### Related types

[zParametricStairsSettings](api-zparametricstairssettings.md)

[zStairsGraphOptimisation](api-zstairsgraphoptimisation.md)

[zVaultSolveReport](../group-statics/api-zvaultsolvereport.md)
