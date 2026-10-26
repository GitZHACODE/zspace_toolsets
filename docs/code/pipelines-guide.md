# Working with toolset pipelines

## Working with toolset pipelines

### 3D printing: zTs3DP

Provide an input through fromMesh or readMesh, set the layer height, field resolution, print width/spacing and other slicing parameters, then computeAll or call the stages explicitly. The declared stages are computeSlicingFeatures, computeSlices, computeUnrolledSlices, computeUnrolledSDFs, computePrintPaths and computePrintMesh. Results include slice meshes, unrolled slices, contours, loops, corner points and print geometry. clearResults resets generated output; consult settings() and the stage APIs when building a partial pipeline.

[zTs3DP API](group-fabrication/api-zts3dp.md)

[Slicer internals exposed by the public API](group-fabrication/api-z3dpslicer.md)

[Print synthesis API](group-fabrication/api-z3dpprintsynthesis.md)

### C-Tubes: zTsCTubes

Assign a directrix and a generatrix, or multiple target generatrices, using Core graphs. Configure zCTubeSettings, run computeTube or a frames/sections variant, then computeStrips and unrollStrips as needed. Read fitReport and stripAreas with the output geometry. Multiple generatrices constrain one tube; they do not request separate sweeps. A non-converged target fit must be reported as a best fit.

[zTsCTubes API](group-fabrication/api-ztsctubes.md)

[C-Tube settings](group-fabrication/api-zctubesettings.md)

### Best Fit FDM: zTsVault

Set the input mesh, supports, tension edges, reference positions, vertex masses/loads, gravity and force-density bounds before solving. Configuration setters invalidate the current solve. computeBestFitForceDensities, stepBestFit and solveBestFit expose different solver entry points; previousBestFit and resetBestFit manage snapshots. Query the report, force densities, residuals and reactions rather than interpreting geometry alone as convergence.

Positive force density denotes tension. Explicit vertex loads combine with mass times gravity. Keep input assignments in the same mesh coordinate/index system. The resultMesh accessor is read-only during solving. Diagram creation and update methods have their own mappings and report; they are not interchangeable with the numerical solve.

[zTsVault API](group-statics/api-ztsvault.md)

[Vault settings and tolerances](group-statics/api-zvaultbestfitsettings.md)

### Read algorithm notes with the API

The existing method pages describe numerical assumptions, limits and verification for Force Polyhedra, Force C-Tubes, panelisation, stairs and the other migrated toolsets. Use those notes alongside the declarations, particularly when diagnosing convergence or comparing against an upstream algorithm.

[Vault method and verification](../vault-best-fit-fdm.md)

[Force Polyhedra method](../force-polyhedra.md)

[C-Tubes method](../ctubes-port.md)
