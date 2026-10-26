# zTsPanelisation

## zSpace::zTsPanelisation

Declared in include/zspace/zToolsets/fabrication/zTsPanelisation.h. Public class/struct declaration reference.

[Header source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/zTsPanelisation.h#L38)

### Public Func

#### compute

{% code collapsedlinecount="10" %}
```cpp
void zSpace::zTsPanelisation::compute(std::vector< zPanelisationSource > &sources, const zPanelisationSettings &settings)
```
{% endcode %}

Parameter descriptions are absent in the current declaration comments for: settings, sources. Use the declaration and implementation source when interpreting these inputs.

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/zTsPanelisation.h#L40)

#### mesh

{% code collapsedlinecount="10" %}
```cpp
zObjectMesh & zSpace::zTsPanelisation::mesh()
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/zTsPanelisation.h#L41)

#### facePanels

{% code collapsedlinecount="10" %}
```cpp
const zIntArray & zSpace::zTsPanelisation::facePanels() const
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/zTsPanelisation.h#L42)

#### faceSources

{% code collapsedlinecount="10" %}
```cpp
const zIntArray & zSpace::zTsPanelisation::faceSources() const
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/zTsPanelisation.h#L43)

#### jointEdges

{% code collapsedlinecount="10" %}
```cpp
const zIntArray & zSpace::zTsPanelisation::jointEdges() const
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/zTsPanelisation.h#L44)

#### jointKinds

{% code collapsedlinecount="10" %}
```cpp
const zIntArray & zSpace::zTsPanelisation::jointKinds() const
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/zTsPanelisation.h#L45)

#### panels

{% code collapsedlinecount="10" %}
```cpp
const std::vector< zPanelInfo > & zSpace::zTsPanelisation::panels() const
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/zTsPanelisation.h#L46)

#### guides

{% code collapsedlinecount="10" %}
```cpp
const zPanelisationGuides & zSpace::zTsPanelisation::guides() const
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/zTsPanelisation.h#L47)

#### components

{% code collapsedlinecount="10" %}
```cpp
int zSpace::zTsPanelisation::components() const
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/zTsPanelisation.h#L48)

#### weldedEdges

{% code collapsedlinecount="10" %}
```cpp
int zSpace::zTsPanelisation::weldedEdges() const
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/zTsPanelisation.h#L49)

#### smallPanels

{% code collapsedlinecount="10" %}
```cpp
int zSpace::zTsPanelisation::smallPanels() const
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/zTsPanelisation.h#L50)

#### resolvedAxis

{% code collapsedlinecount="10" %}
```cpp
int zSpace::zTsPanelisation::resolvedAxis() const
```
{% endcode %}

[Declaration source](https://github.com/GitZHACODE/zspace_toolsets/blob/90f88442ee45b1760b7333a6f4aedb1f1e6ee80a/include/zspace/zToolsets/fabrication/zTsPanelisation.h#L51)

### Related types

[zPanelInfo](api-zpanelinfo.md)

[zPanelisationGuides](api-zpanelisationguides.md)

[zPanelisationSettings](api-zpanelisationsettings.md)

[zPanelisationSource](api-zpanelisationsource.md)
