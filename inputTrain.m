% geometry used for training
clear; clc; close all;

% Problem definition
problem.polarization="TE";
% Size of simulation domain
problem.domain.Lx=24e-6;
problem.domain.Ly=20e-6;
problem.mesh.Nx=700; % number of grid points in x
problem.mesh.Ny=900; % number of grid points in y
nworkers=30;
problem.mesh.subpixelSamples=10; % subsampling
problem.material.backgroundEr=8.1;
problem.pml.thickness
problem.pml.order=3;
problem.pml.reflection=1e-7;
problem.geometry=G;
problem.sources=Sources;
problem.probes=Probes;


% Offline training settings
training.outputDir="./TrainingData";
training.deletePrevious=true;
training.wavelengthRange=[]; 
training.frequencyBins=3; 
training.geometriesPerBin=90;
training.frequenciesPerBin=40;
training.randomSeed=1;


training.parameters(1).device=2;
training.parameters(1).field="er";
training.parameters(1).range=[];

% reduced order model settings
rom.pod.maxModes=40;
rom.pod.energyTolerance=1e-8;

rom.mdeim.maxModes=100;
rom.mdeim.energyTolerance=1e-10;


rom.localModels.parameters(1).device=2;
rom.localModels.parameters(1).field="innerRadius";
rom.localModels.parameters(1).edges = linspace(450e-9,550e-9,4);



offline.generateSnapshots()