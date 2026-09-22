# Data Structures


###  Configuration Structures 

### `problem`

The problem structure contains physical and numerical definitions of the full order problem.

| Field | Type or options | Units | Meaning|
|---|---|---|---|
|problem.polarization|"TE" or "TM" |- |Selects the electromagnetic polarization used|
|problem.domain|struct|-|Defines simulation domain|
|problem.mesh|struct|-|Defines the discretization scheme|
|problem.pml|struct|-|Defines PML properties|
|problem.geometry|struct|-|Defines the geometry|
|problem.sources|struct|-|Defines sources|
|problem.probes|struct|-|Defines probes|

#### `problem.domain`

Contains information regarding the simulation domain

| Field | Type or options | Units | Meaning|
|---|---|---|---|
|problem.domain.Lx|positive scalar|m|Length of domain in x-direction |
|problem.domain.Ly|positive scalar|m|Length of domain in y-direction|


#### `problem.mesh`

Defines the discretization scheme used

| Field | Type or options | Units | Meaning|
|---|---|---|---|
|problem.mesh.Nx|positive integer|-|points in x-direction |
|problem.mesh.Ny|positive integer|-|points in y-direction|
|problem.mesh.subpixelSamples|positive integer|-|points used in subpixel sampling|

#### `problem.pml`

Contains parameters regarding the PML. See PML documentation for equations used.

| Field | Type or options | Units | Meaning|
|---|---|---|---|
|problem.pml.thickness|positive scalar|m|thickness of the PML|
|problem.pml|order|positive integer|PML order|
|problem.pml.reflection|positive scalar|-|reflection coefficient|

#### `problem.geometry`, `problem.sources`,`problem.probes`

These define the geometry, sources and probes used in the current simulation.
See geometry.md, sources.md and probes.md.


| Field | Type or options | Units | Meaning|
|---|---|---|---|
|problem.geometry|struct array|-|holds geometry elements|
|problem.sources|struct array|-|holds source elements|
|problem.probes|struct array|-|holds probe elements|


### `training`

This struct defines parameters used to collect snapshots


| Field | Type or options | Units | Meaning|
|---|---|---|---|
|training.outputDir|"string"|-|Directory name to store snapshots|
|training.deletePrevious|true or false|-|If set to true this will delete all snapshots in the training.outputDir|
|training.wavelengthRange|double (2)|1/m|Range of frequencies to train|
|training.frequenciesPerBin|integer|-|Number of frequency models to create|
|training.randomSeed|integer|-|Seed for rng|
|training.parameters|struct array|-|Holds the fields to train|


#### `training.parameters`

Array of structs specifying the DoF to collect snapshots for.
Example:

```
training.parameters(1).device=2;
training.parameters(1).field="er";
training.parameters(1).range=[8.1,8.2];
```

| Field | Type or options | Units | Meaning|
|---|---|---|---|
|training.parameters(i).device|positive integer|-|Device ID|
|training.parameters(i).field|"string"|-|Name of Field|
|training.parameters(i).range|double (2)|-|Training range|


