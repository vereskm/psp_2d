# Runtime Data
This file documents structures returned by functions.
Users should look towards the README.md and inputs.md to determine how to use the code.
|Symbol|Definition|Meaning|
|---|---|---|
|nGeometry|training.geometriesPerBin|Number of unique sampled geometry|
|nFrequency|training.freqenciesPerBin|Number of sampled frequencies in each bin|
|nBins|training.frequencyBins|Number of frequency bins|
|nFeat|numel(training.parameters)|Number of varied training parameters|
|nTrain|nFeat $\times$ nFrequency| Number of snapshots for each bin|
|nTotal|nTrain $\times$ nBin| Number of snapshots for each bin|
|nDoF| $N_x\times N_y$| Total number of unknowns|
|nModes| ||



## `snapshots`

Returned from `generateSnapshots(problem,training)`.
Hold snapshot data metadata regarding its construction.

| Field | Type or options | Units | Meaning|
|---|---|---|---|
|snapshots.wBinEdges|"TE" or "TM" |- |Selects the electromagnetic polarization used|
|snapshots.wSamples|struct|-|Defines simulation domain|
|problem.mesh|struct|-|Defines the discretization scheme|
|problem.parameterSamples|struct|-|Defines PML properties|
|problem.geometry|struct|-|Defines the geometry|
|problem.sources|struct|-|Defines sources|
|problem.probes|struct|-|Defines probes|
