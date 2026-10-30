% geometry used for training
clear; clc; close all;
projectRoot = fileparts(mfilename('fullpath'));
 addpath(fullfile(projectRoot,'src'));
device;
% Problem definition
problem.polarization="TE";
% Size of simulation domain
problem.domain.Lx=Lx;
problem.domain.Ly=Ly;
problem.mesh.Nx=700; % number of grid points in x
problem.mesh.Ny=900; % number of grid points in y
nworkers=30;
pool = gcp('nocreate');
if isempty(pool)
    parpool('local',nworkers);
elseif pool.NumWorkers ~= nworkers
    delete(pool);
    parpool('local',nworkers);
end


problem.mesh.subpixelSamples=10; % subsampling
problem.material.backgroundEr=8.1;
problem.pml.thickness=dpml;
problem.pml.order=3;
problem.pml.reflection=1e-7;
problem.geometry=G;
problem.sources=Sources;
problem.probes=Probes;


% Offline training settings
training.outputDir = fullfile(projectRoot,'TrainingData');
training.deletePrevious=true;
training.wavelengthRange=[lambda1,lambda2]; 
training.frequencyBins=1; 
training.geometriesPerBin=30;
training.frequenciesPerBin=20;
training.randomSeed=1;


training.parameters(1).device=2;
training.parameters(1).field="er";
training.parameters(1).range=[esilicon.*(1-.01),esilicon.*(1+.01)];

% reduced order model settings
rom.pod.maxModes=40;
rom.outputDir = fullfile(projectRoot,'Model');

rom.mdeim.geometriesPerBin=100;
rom.mdeim.frequenciesPerBin=100;
rom.mdeim.randomSeed=1;
rom.mdeim.maxModes=100;
rom.mdeim.energyTolerance=1e-10;


%rom.localModels.parameters(1).device=2;
%rom.localModels.parameters(1).field="innerRadius";
%rom.localModels.parameters(1).edges = linspace(450e-9,550e-9,4);
%rom.localModels.parameters(1).unit="nm";

rom.rbf.kernel="cubic";
rom.rbf.polynomialDegree=1;
rom.rbf.normalize=true;
rom.rbf.maxCenters=300;
rom.rbf.regularization=1e-11;

[snapshots,gridData]=offline.generateSnapshots(problem,training);
POD = offline.buildLibrary(problem,training,rom,snapshots,gridData);
save(fullfile(projectRoot,'Problem.mat'),'problem');

