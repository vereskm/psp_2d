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
|problem.sources|struct|-|Defines sources|
|problem.probes|struct|-|Defines probes|

#### 'problem.domain'

Contains information regarding the simulation domain

| Field | Type or options | Units | Meaning|
|---|---|---|---|
|problem.domain.Lx|positive scalar|m|Length of domain in x-direction |
|problem.domain.Ly|positive scalar|m|Length of domain in y-direction|


#### 'problem.mesh'

Defines the discretization scheme used

| Field | Type or options | Units | Meaning|
|---|---|---|---|
|problem.mesh.Nx|positive integer|-|points in x-direction |
|problem.mesh.Ny|positive integer|-|points in y-direction|
|problem.mesh.subpixelSamples|positive integer|-|points used in subpixel sampling|


