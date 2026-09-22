function []=buildLibrary(problem,training,rom,snapshots)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lx = settings.Lx;
Ly = settings.Ly;
Nx = settings.Nx;
Ny = settings.Ny;
pol    = settings.pol;
numBin = settings.numBin;
neig   = settings.neig;
nfeat  = numel(TrainingParams);

x = linspace(0,Lx,Nx);
y = linspace(0,Ly,Ny);
[X,Y] = meshgrid(x,y);

outputDirectory = "./Lib";
if ~isfolder(outputDirectory)
    mkdir(outputDirectory);
end
pool = gcp("nocreate");
if isempty(pool)
    pool = parpool("Processes",nworkers);
elseif pool.NumWorkers ~= nworkers
    delete(pool);
    pool = parpool("Processes",nworkers);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get frequencies
wmin = min(wRange);
wmax = max(wRange);
wBinEdges=linspace(wmin,wmax,numBin+1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Both names must exist before the parfor below because MATLAB classifies
% every variable referenced in the loop, including the inactive
% polarization branch.
etaTE = cell(0,1);
etaTM = cell(0,1);

if pol == "TE" 
    [etaTE,lambdaTE,PODDiagnosticsTE] = buildPOD( ...
        trainDir,"TE",settings.numBin,settings.neig, ...
        settings,TrainingParams,Sources,wBinEdges);
end

if pol == "TM" 
    [etaTM,lambdaTM,PODDiagnosticsTM] = buildPOD( ...
        trainDir,"TM",settings.numBin,settings.neig, ...
        settings,TrainingParams,Sources,wBinEdges);
end


nfeat=numel(TrainingParams);
rangesM=cell(nfeat,1);
gridMat=cell(nfeat,1);
for i=1:nfeat
    sR=TrainingParams(i).range(1);
    eR=TrainingParams(i).range(2);
    nMatI=TrainingParams(i).nMatI;
    rangesM{i}=linspace(sR,eR,nMatI);
end

% Get the mesh needed to train the matrices
[gridMat{:}]=ndgrid(rangesM{:});
nGTMat=numel(gridMat{1});
% This compute matrix points to project to store a precomputed grid
% Allocate only reduced quantities
if pol == "TE" 
    ArTE = cell(numBin,nGTMat,numel(wTp));
    brTE = cell(numBin,nGTMat,numel(wTp));
    ProbeTE=cell(numBin,1);
end

if pol == "TM"
    ArTM = cell(numBin,nGTMat,numel(wTp));
    brTM = cell(numBin,nGTMat,numel(wTp));
    ProbeTM=cell(numBin,1);
end


parfor gi = 1:nGTMat
    Gref = G;

    fprintf( ...
        "Building and projecting matrices: geometry %d/%d\n", ...
        gi,nGTMat);

    for i = 1:nfeat
        d = TrainingParams(i).device;
        fieldname = TrainingParams(i).name;
        Gref(d).(fieldname) = gridMat{i}(gi);
    end

    e = gen_e(settings,X,Y,Gref);
    [sgmx,sgmy] = PML(settings,X,Y,e);

    localArTE = cell(numBin,numel(wTp));
    localBrTE = cell(numBin,numel(wTp));
    localArTM = cell(numBin,numel(wTp));
    localBrTM = cell(numBin,numel(wTp));

    for wi = 1:numel(wTp)
        wc = wTp(wi);

        % A shared bin-edge frequency belongs to both neighboring
        % interpolants, so project it onto both POD bases.
        matchingBins = find( ...
            wc >= wBinEdges(1:end-1) & ...
            wc <= wBinEdges(2:end));

        assert(~isempty(matchingBins), ...
            "wTp(%d) lies outside wRange.",wi);

        if pol == "TE" 
            [phi,beta] = MakeSource( ...
                e,wc,X,Y,Sources,settings,"TE",sgmx,sgmy);

            A = TE_FDFD(settings,wc,e,sgmx,sgmy);
            bv = QAAQ(X,Y,A,phi,beta,Sources);

            % Project before the next frequency replaces A.  Usually this
            % loop has one entry; at an internal edge it has two.
            for b = matchingBins
                eta = etaTE{b};
                localArTE{b,wi} = eta' * (A * eta);
                localBrTE{b,wi} = eta' * bv;
            end
        end

        if pol == "TM" 
            [phi,beta] = MakeSource( ...
                e,wc,X,Y,Sources,settings,"TM",sgmx,sgmy);

            A = TM_FDFD(settings,wc,e,sgmx,sgmy);
            bv = QAAQ(X,Y,A,phi,beta,Sources);

            for b = matchingBins
                eta = etaTM{b};
                localArTM{b,wi} = eta' * (A * eta);
                localBrTM{b,wi} = eta' * bv;
            end
        end
    end

    if pol == "TE"
        ArTE(:,gi,:) = reshape( ...
            localArTE,numBin,1,numel(wTp));

        brTE(:,gi,:) = reshape( ...
            localBrTE,numBin,1,numel(wTp));
    end

    if pol == "TM" 
        ArTM(:,gi,:) = reshape( ...
            localArTM,numBin,1,numel(wTp));

        brTM(:,gi,:) = reshape( ...
            localBrTM,numBin,1,numel(wTp));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build Prob library
GDefault=G; % get default permativities for the probes 
% (note the probes should not be in areas effected by parametric changes)
eD=gen_e(settings,X,Y,GDefault);
[sgmxD,sgmyD] = PML(settings,X,Y,eD);

  if pol == "TE"
      for b=1:numBin
          [CP,CN]=ProbLibPOD(eD,sgmxD,sgmyD,wTp,X,Y,etaTE{b},Probes,settings,pol);
          ProbeTE{b}.CP=CP;
          ProbeTE{b}.CN=CN;
      end
  elseif pol=="TM"
    for b=1:numBin
          [CP,CN]=ProbLibPOD(eD,sgmxD,sgmyD,wTp,X,Y,etaTM{b},Probes,settings,pol);
          ProbeTM{b}.CP=CP;
          ProbeTM{b}.CN=CN;
    end
  end


fprintf("Done building and projecting interpolation matrices.\n");

  fprintf( ...
        "Done Building Matrices\n");

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Save Meta data
if pol=="TE"
    save('./Lib/PODModes.mat', ...
        'etaTE','lambdaTE','PODDiagnosticsTE','wBinEdges');
elseif pol=="TM"
    save('./Lib/PODModes.mat', ...
        'etaTM','lambdaTM','PODDiagnosticsTM','wBinEdges');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot eigenvalues
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot the eigenvalues
        EigValPlot=figure;
        hold on;
        box on;
        ax=gca;
        set(ax,'Linewidth',2);
        if pol=="TE"
            for i=1:numBin
               % Bin limits in THz
               fL = wBinEdges(i)   /(2*pi)/1e12;
               fH = wBinEdges(i+1)/(2*pi)/1e12;
               semilogy(lambdaTE(i,:)/lambdaTE(i,1),'LineWidth',2,"DisplayName",sprintf('%.1f-%.1f THz',fL,fH));
            end
            
        elseif pol=="TM"
            for i=1:numBin
               % Bin limits in THz
               fL = wBinEdges(i)   /(2*pi)/1e12;
               fH = wBinEdges(i+1)/(2*pi)/1e12;
               semilogy(lambdaTM(i,:)/lambdaTM(i,1),'LineWidth',2,"DisplayName",sprintf('%.1f-%.1f THz',fL,fH));
            end

      
        end
        xlabel("Index","FontSize",15);
        ylabel("POD Eigenvalue","FontSize",15);
        box on;
        set(ax,'FontSize',15);
        set(gca,"YScale","log");
        legend("Location","Best","FontSize",15);
        savefig(EigValPlot,'./Lib/Eigenvalues.fig');
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % BUILD INTERPOLATOR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf("Printing Results\n");
    if  pol=="TE"
        InterpTE = struct([]);
    
        for b=1:numBin
             [InterpTE(b).A,InterpTE(b).b]=BuildInterp(ArTE,brTE,b,rangesM,wTp,wBinEdges);
              InterpTE(b).wRange = wBinEdges(b:b+1);
        end
        save('./Lib/InterpTE.mat', ...
        'InterpTE', ...
        'wBinEdges', ...
        'rangesM', ...
        'TrainingParams', ...
        'wTp','ProbeTE');
    elseif pol=="TM"
        InterpTM = struct([]);
        for b=1:numBin
             [InterpTM(b).A,InterpTM(b).b]=BuildInterp(ArTM,brTM,b,rangesM,wTp,wBinEdges);
              InterpTM(b).wRange = wBinEdges(b:b+1);
        end
         save('./Lib/InterpTM.mat', ...
        'InterpTM', ...
        'wBinEdges', ...
        'rangesM', ...
        'TrainingParams', ...
        'wTp','ProbeTM');
    end


end



function [eta,lambda,diagnostics] = buildPOD(trainDir,pol,numBin,neig, ...
    expectedSettings,expectedTrainingParams,expectedSources,expectedBinEdges)

pol = upper(string(pol));
files = dir(fullfile(trainDir,pol + "_*.mat"));
fieldName = "Output" + pol;

assert(~isempty(files),"No %s training files found.",pol);

eta   = cell(numBin,1);
lambda = nan(numBin,neig);
diagnostics = repmat(struct( ...
    "Polarization",pol, ...
    "Bin",NaN, ...
    "NumSnapshots",NaN, ...
    "NumRetainedModes",NaN, ...
    "RetainedEnergyFraction",NaN, ...
    "LastToFirstEigenvalueRatio",NaN),numBin,1);

for b = 1:numBin
    snapshots = {};

    for f = 1:numel(files)
        filePath = fullfile(files(f).folder,files(f).name);
        data = load(filePath,fieldName,"settings","TrainingParams", ...
            "Sources","wBinEdges","wSamples");

        requiredMetadata = [fieldName,"settings","TrainingParams", ...
            "Sources","wBinEdges","wSamples"];
        for metadataName = requiredMetadata
            assert(isfield(data,metadataName), ...
                "Training file %s lacks metadata '%s'. Retrain that file before combining it.", ...
                filePath,metadataName);
        end

        % Seeds, sample counts, sampled values, and parameter ranges may
        % differ between files.  Only the discretization and model layout
        % must agree for their field snapshots to share one POD basis.
        structuralFields = ["Lx","Ly","Nx","Ny","eback", ...
            "m","dPML","R","Nsub","numBin"];
        for settingName = structuralFields
            assert(isfield(data.settings,settingName) && ...
                isequaln(data.settings.(settingName), ...
                         expectedSettings.(settingName)), ...
                "Training file %s has incompatible setting '%s'.", ...
                filePath,settingName);
        end
        assert(data.settings.pol == pol, ...
            "Training file %s has incompatible polarization metadata.",filePath);
        assert(isequaln(data.Sources,expectedSources), ...
            "Training file %s uses a different source definition.",filePath);
        assert(isequaln(data.wBinEdges,expectedBinEdges), ...
            "Training file %s uses different frequency-bin edges.",filePath);

        assert(numel(data.TrainingParams) == numel(expectedTrainingParams), ...
            "Training file %s uses a different number of trained features.",filePath);
        for featureIndex = 1:numel(expectedTrainingParams)
            assert(data.TrainingParams(featureIndex).device == ...
                    expectedTrainingParams(featureIndex).device && ...
                   string(data.TrainingParams(featureIndex).name) == ...
                    string(expectedTrainingParams(featureIndex).name), ...
                "Training file %s uses a different definition for feature %d.", ...
                filePath,featureIndex);
        end

        output = data.(fieldName);
        assert(size(output,2) == numBin && ...
               isequal(size(data.wSamples),size(output)), ...
            "Training file %s has inconsistent snapshot/bin dimensions.",filePath);
        expectedFieldLength = (expectedSettings.Nx-2)*(expectedSettings.Ny-2);
        assert(all(cellfun(@(u) isvector(u) && ...
            numel(u) == expectedFieldLength,output(:))), ...
            "Training file %s contains snapshots with the wrong grid dimension.",filePath);
        assert(all(cellfun(@(u) all(isfinite(u),"all"),output(:))), ...
            "Training file %s contains NaN or Inf snapshot values.",filePath);
        snapshots = [snapshots; output(:,b)]; %#ok<AGROW>
    end

    snapshots = cellfun(@(u) u(:),snapshots, ...
        "UniformOutput",false);

    snapshotMatrix = [snapshots{:}];
    k = min([neig,size(snapshotMatrix)]);

    [eta{b},S,~] = svds(snapshotMatrix,k);
    lambda(b,1:k) = diag(S).^2;

    retainedEnergy = sum(lambda(b,1:k))/norm(snapshotMatrix,"fro")^2;
    diagnostics(b).Bin = b;
    diagnostics(b).NumSnapshots = size(snapshotMatrix,2);
    diagnostics(b).NumRetainedModes = k;
    diagnostics(b).RetainedEnergyFraction = retainedEnergy;
    diagnostics(b).LastToFirstEigenvalueRatio = lambda(b,k)/lambda(b,1);

    fprintf("POD %s bin %d: %d snapshots, %d modes, " + ...
        "retained energy %.8f, lambda(last)/lambda(1) %.3g.\n", ...
        pol,b,diagnostics(b).NumSnapshots,k,retainedEnergy, ...
        diagnostics(b).LastToFirstEigenvalueRatio);
end

end
