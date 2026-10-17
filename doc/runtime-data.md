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
|snapshots.wBinEdges|array (double) |Edges of frequency bins|
|snapshots.wSamples|array: (nTrain,nBins)|-|Frequencies trained per bin|
|snapshots.parameterSamples|array: (nTrain,nFeat,nBins)|-|holds the features trained for each snapshot|
|snapshots.fields|cell array: (nTrain,nBins)|-|Holds the snapshots for each frequency bin|

