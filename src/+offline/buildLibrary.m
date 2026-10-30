
function [POD] = buildLibrary(problem,training,rom,snapshots,gridData)
    libraryTimer = tic;
    fprintf('[Library] Starting POD and MDEIM library construction.\n');
    %BUILDLIBRARY
    % for this code you can either pass in the data of the file of the data
    if ischar(snapshots)|| isstring(snapshots)
        data=load(snapshots,'snapshots');
        snapshots=data.snapshots;
    end
    if nargin<5
        % if the user did not pass in the gridData (that is OK) we do it
        % here. Its OK!
        gridData=mesh.create(problem);
    end
    if ~isfolder(rom.outputDir)
        mkdir(rom.outputDir);
    end

    % TODO: now it only support 1 group
    %AHHHHHHHH, this is just so complicated, my head is spinning. Please
    %get me some help. Maybe it wasn't the best to refactor everything. I
    %am in some sort of hell.

    LocalDim=[]; % number of local dimensions (can only be 1 for now)
    nBins = training.frequencyBins;
    if isfield(rom,'localModels')
        assert(numel(rom.localModels.parameters)==1,"For now the code only supports 1 local parameter");
        nLoc=numel(rom.localModels.parameters);
        LocalDim=zeros(nLoc,1);
        for i=1:nLoc
            LocalDim(i)=numel(rom.localModels.parameters(i).edges)-1;
        end
    end
    if isempty(LocalDim)
        POD.modes=cell(nBins,1);
    else
        POD.modes = cell([nBins,LocalDim]);
    end
    % match the same size
    POD.lambda=cell(size(POD.modes));
    POD.wBinEdges=snapshots.wBinEdges;
    POD.trainingParameters=training.parameters;
    POD.localFeatureIndex=[];
    POD.localEdges=[];
    % now collect snapshots and generate POD modes
    fprintf('[POD] Generating bases for %d frequency bins.\n',nBins);
    maxModes=rom.pod.maxModes;
    for b=1:nBins
        if isfield(rom,'localModels')
            % find the index of the training
            device=rom.localModels.parameters(1).device;
            field=rom.localModels.parameters(1).field;
            Index=[];
            for p=1:numel(training.parameters)
                if device==training.parameters(p).device && field==training.parameters(p).field
                    Index=p;
                end
            end
            assert(~isempty(Index),"Local parameter is not in training.parameters.");
            edges=rom.localModels.parameters(1).edges;
            POD.localFeatureIndex=Index;
            POD.localEdges=edges;
            % get the modes
            localValues=snapshots.parameterSamples(:,Index,b);
            localIndicies=discretize(localValues,edges);
            for localBin=1:LocalDim(1)
                mask=(localIndicies==localBin);
                localSnapshots=snapshots.fields(mask,b);
                fprintf('[POD] Bin %d/%d, region %d: SVD of %d snapshots.\n',b,nBins,localBin,numel(localSnapshots));
                [POD.modes{b,localBin},...
                    POD.lambda{b,localBin}]=generatePODModes(localSnapshots,maxModes);
            end
        else
            fprintf('[POD] Bin %d/%d: SVD of %d snapshots.\n',b,nBins,size(snapshots.fields,1));
            [POD.modes{b,1},POD.lambda{b,1}]=generatePODModes(snapshots.fields(:,b),maxModes);
        end
    end

    for b=1:nBins
        for d=1:size(POD.modes,2)
            fprintf('[POD] Bin %d, region %d: retained %d modes.\n',b,d,size(POD.modes{b,d},2));
        end
    end

    % now plot the eigenvalues onto a single plot
    eigenmodeFig=figure;
    hold on;
    for b=1:nBins

        if isfield(rom,'localModels')
            for d=1:LocalDim
                startr=rom.localModels.parameters(1).edges(d);
                endr=rom.localModels.parameters(1).edges(d+1);
                if(isfield(rom.localModels.parameters(1),"unit"));
                    label=string(b)+":"+string(startr)+"-"+string(endr)+...
                        rom.localModels.parameters(1).unit;
                else
                    label=string(b)+":"+string(startr)+"-"+string(endr);
                end
                semilogy(POD.lambda{b,d}/max(max(POD.lambda{b,d}),realmin),'-',"LineWidth",2,'DisplayName',label);
            end
        else
            label=string(b);
            axis([1,inf,-inf,inf]);
            semilogy(POD.lambda{b,1}/max(max(POD.lambda{b,1}),realmin),'-',"LineWidth",2,'DisplayName',label);
        end
    end
    ax=gca;
    set(ax,'Yscale','log');
    set(ax,"Linewidth",2);
    set(ax,"Fontsize",14);
    xlabel('Mode index',"FontSize",18);
    ylabel('Normalized eigenvalue (\lambda/\lambda_{max})',"FontSize",18);
    box on;
    legend("show","FontSize",18);
    POD.nBins=nBins;
    % save the eigenfigure
    savefig(eigenmodeFig,fullfile(rom.outputDir,'eigenValues.fig'));

    % now we will create the reduced matrix to store!!!

    % now create the A matrices. See write up how I did it
    pol=problem.polarization;
    nFeat=numel(training.parameters);
    % get the geometries to train A and Q and DL and interpolartions
    nGeometry=rom.mdeim.geometriesPerBin;
    nFrequency=rom.mdeim.frequenciesPerBin;
    Adata=cell(size(POD.modes));
    POD.Abasis=cell(size(POD.modes));
    POD.Aproj=cell(size(POD.modes));
    POD.Api=cell(size(POD.modes));
    Qdata=cell(size(POD.modes));
    POD.Qbasis=cell(size(POD.modes));
    POD.Qproj=cell(size(POD.modes));
    POD.Qpi=cell(size(POD.modes));
    DLdata=cell(nBins,1);
    POD.DLbasis=cell(size(POD.modes));
    POD.DLproj=cell(size(POD.modes));
    POD.DLpi=cell(size(POD.modes));
    POD.interpolationSpace='deimCoefficients';
    POD.Arbf=cell(size(POD.modes));
    POD.Qrbf=cell(size(POD.modes));
    POD.DLcoefficients=cell(size(POD.modes));
    POD.geometrySamples=cell(size(POD.modes));
    POD.wpoints=cell(nBins,1);
    POD.referenceFrequency=zeros(nBins,1);

    lower=zeros(1,nFeat);
    upper=zeros(1,nFeat);
    for i=1:nFeat
        lower(i)=training.parameters(i).range(1);
        upper(i)=training.parameters(i).range(2);
    end
    assert(nFrequency>=2,'At least two MDEIM frequency samples are required.');
    if isfield(rom.mdeim,'randomSeed')
        rng(rom.mdeim.randomSeed);
    end

    for b=1:nBins
        fprintf('[MDEIM] Bin %d/%d: %d geometry samples, %d frequency samples.\n',b,nBins,nGeometry,nFrequency);
        A=cell(nGeometry,1); 
        Q=cell(nGeometry,1);
        DL=cell(nFrequency,1);
        DLPOD=cell(nFrequency,1);
        wmin=snapshots.wBinEdges(b);
        wmax=snapshots.wBinEdges(b+1);
        wpoints=linspace(wmin,wmax,nFrequency);
        l=(wmin+wmax)/2; % get the center frequency
        POD.wpoints{b}=wpoints;
        POD.referenceFrequency(b)=l;
        geometrySamples=zeros(nGeometry,nFeat);
        % see if there is a local mode
        geometryUnit=lhsdesign(nGeometry,nFeat);
        % we are assuming that the geometry changes does not impact the
        % PML
        e = geometry.generatePermittivity(problem,gridData);
        [sgmx,sgmy] = fdfd.buildPML(problem,gridData,e);

        fprintf('[MDEIM] Building full frequency corrections.\n');
        parfor gi=1:nFrequency
            w=wpoints(gi);
            Gref=problem.geometry;
            % build the correction factors
            DL{gi} = buildFrequencyCorrection(gridData,w,l,e,sgmx,sgmy,problem);
            fprintf('[MDEIM] Bin %d: frequency correction %d/%d finished.\n',b,gi,nFrequency);
        end

        if isfield(rom,'localModels')
            % find the index of the training
            device=rom.localModels.parameters(1).device;
            field=rom.localModels.parameters(1).field;
            Index=[];
            edges=rom.localModels.parameters(1).edges;
            for p=1:numel(training.parameters)
                if device==training.parameters(p).device && field==training.parameters(p).field
                    Index=p;
                end
            end
            for i = 1:nFeat
                if i==Index
                    continue;
                else
                    low = training.parameters(i).range(1);
                    high = training.parameters(i).range(2);
                    geometrySamples(:,i) = low + (high-low)*geometryUnit(:,i);
                end
            end

            for d=1:LocalDim
                % assign the local dimension
                low=edges(d);
                high=edges(d+1); 
                geometrySamples(:,Index) = low + (high-low)*geometryUnit(:,Index);
                V=POD.modes{b,d}; % get the needed modes 
                parfor gi=1:nFrequency
                    DLPOD{gi}=V'*DL{gi}*V;
                end

                fprintf('[MDEIM] Sampling and projecting A and Q.\n');

                parfor gi=1:nGeometry
                    problemRef=problem;
                    Gref=problem.geometry;
                    for i=1:nFeat
                        device= training.parameters(i).device;
                        fieldname = training.parameters(i).field;
                        Gref(device).(fieldname) = geometrySamples(gi,i);
                    end
                    problemRef.geometry=Gref;
                    e = geometry.generatePermittivity(problemRef,gridData);

                    A{gi}=V'*buildReferenceOperator(gridData,l,e,sgmx,sgmy,problem)*V;
                    Q{gi}=V'*buildMassMatrix(gridData,e,problem)*V;
                    fprintf('[MDEIM] Bin %d: geometry sample %d/%d finished.\n',b,gi,nGeometry);
                end
                Adata{b,d}=A;
                fprintf('[DEIM] Bin %d: compressing A.\n',b);
                [POD.Abasis{b,d},POD.Aproj{b,d},POD.Api{b,d}] = greedyDEIM(A,rom.mdeim.maxModes,rom.mdeim.energyTolerance);
                fprintf('[DEIM] A: retained %d matrix modes.\n',size(POD.Abasis{b,d},2));
                Qdata{b,d}=Q;
                fprintf('[DEIM] Bin %d: compressing Q.\n',b);
                [POD.Qbasis{b,d},POD.Qproj{b,d},POD.Qpi{b,d}] = greedyDEIM(Q,rom.mdeim.maxModes,rom.mdeim.energyTolerance);
                fprintf('[DEIM] Q: retained %d matrix modes.\n',size(POD.Qbasis{b,d},2));
                DLdata{b,d}=DLPOD;
                fprintf('[DEIM] Bin %d: compressing DL.\n',b);
                [POD.DLbasis{b,d},POD.DLproj{b,d},POD.DLpi{b,d}] = greedyDEIM(DLPOD,rom.mdeim.maxModes,rom.mdeim.energyTolerance);
                fprintf('[DEIM] DL: retained %d matrix modes.\n',size(POD.DLbasis{b,d},2));

                % Fit geometry-to-DEIM-coefficient maps for A and Q.
                fprintf('[Interpolation] Fitting A/Q coefficient RBFs and storing frequency coefficients.\n');
                localLower=lower;
                localUpper=upper;
                localLower(Index)=edges(d);
                localUpper(Index)=edges(d+1);
                POD.geometrySamples{b,d}=geometrySamples;
                POD.Arbf{b,d}=fitOperatorRBF(A,POD.Api{b,d},POD.Aproj{b,d}, ...
                    geometrySamples,localLower,localUpper,rom);
                POD.Qrbf{b,d}=fitOperatorRBF(Q,POD.Qpi{b,d},POD.Qproj{b,d}, ...
                    geometrySamples,localLower,localUpper,rom);
                % Interpolate correction coefficients in frequency; convert offline.
                DLcolumns=cellfun(@(M) M(:),DLPOD,'UniformOutput',false);
                DLmatrix=[DLcolumns{:}];
                POD.DLcoefficients{b,d}=(POD.DLproj{b,d}*DLmatrix(POD.DLpi{b,d},:)).';
            end
        else
            for i = 1:nFeat

                low = training.parameters(i).range(1);
                high = training.parameters(i).range(2);
                geometrySamples(:,i) = low + (high-low)*geometryUnit(:,i);
            end
            V=POD.modes{b}; % get the needed modes 
            parfor gi=1:nFrequency
                DLPOD{gi}=V'*DL{gi}*V;
            end

            fprintf('[MDEIM] Sampling and projecting A and Q.\n');

            parfor gi=1:nGeometry
                problemRef=problem;
                Gref=problem.geometry;
                for i=1:nFeat
                    device= training.parameters(i).device;
                    fieldname = training.parameters(i).field;
                    Gref(device).(fieldname) = geometrySamples(gi,i);
                end
                problemRef.geometry=Gref;
                e = geometry.generatePermittivity(problemRef,gridData);

                A{gi}=V'*buildReferenceOperator(gridData,l,e,sgmx,sgmy,problem)*V;
                Q{gi}=V'*buildMassMatrix(gridData,e,problem)*V;
                fprintf('[MDEIM] Bin %d: geometry sample %d/%d finished.\n',b,gi,nGeometry);
            end
            Adata{b}=A;
            fprintf('[DEIM] Bin %d: compressing A.\n',b);
            [POD.Abasis{b},POD.Aproj{b},POD.Api{b}] = greedyDEIM(A,rom.mdeim.maxModes,rom.mdeim.energyTolerance);
            fprintf('[DEIM] A: retained %d matrix modes.\n',size(POD.Abasis{b},2));
            Qdata{b}=Q;
            fprintf('[DEIM] Bin %d: compressing Q.\n',b);
            [POD.Qbasis{b},POD.Qproj{b},POD.Qpi{b}] = greedyDEIM(Q,rom.mdeim.maxModes,rom.mdeim.energyTolerance);
            fprintf('[DEIM] Q: retained %d matrix modes.\n',size(POD.Qbasis{b},2));
            DLdata{b}=DLPOD;
            fprintf('[DEIM] Bin %d: compressing DL.\n',b);
            [POD.DLbasis{b},POD.DLproj{b},POD.DLpi{b}] = greedyDEIM(DLPOD,rom.mdeim.maxModes,rom.mdeim.energyTolerance);
            fprintf('[DEIM] DL: retained %d matrix modes.\n',size(POD.DLbasis{b},2));

            % Fit geometry-to-DEIM-coefficient maps for A and Q.
            fprintf('[Interpolation] Fitting A/Q coefficient RBFs and storing frequency coefficients.\n');
            POD.geometrySamples{b,1}=geometrySamples;
            POD.Arbf{b,1}=fitOperatorRBF(A,POD.Api{b},POD.Aproj{b}, ...
                geometrySamples,lower,upper,rom);
            POD.Qrbf{b,1}=fitOperatorRBF(Q,POD.Qpi{b},POD.Qproj{b}, ...
                geometrySamples,lower,upper,rom);
            % Interpolate correction coefficients in frequency; convert offline.
            DLcolumns=cellfun(@(M) M(:),DLPOD,'UniformOutput',false);
            DLmatrix=[DLcolumns{:}];
            POD.DLcoefficients{b,1}=(POD.DLproj{b}*DLmatrix(POD.DLpi{b},:)).';
        end
    end

    fprintf('[Ports] Precomputing reduced port maps.\n');
    POD.ports=cell(size(POD.modes)); % define the ports!
    e = geometry.generatePermittivity(problem,gridData);
    [sgmx,sgmy] = fdfd.buildPML(problem,gridData,e);
    for b=1:nBins
        wpoints = POD.wpoints{b};
        if isfield(rom,'localModels')
            for d=1:LocalDim
                V=POD.modes{b,d};
                [CP,CN]=ports.buildReducedProbe(e,sgmx,sgmy,wpoints,gridData,V,problem);
                POD.ports{b,d}.CP=CP;
                POD.ports{b,d}.CN=CN;
            end
        else
            V=POD.modes{b};
            [CP,CN]=ports.buildReducedProbe(e,sgmx,sgmy,wpoints,gridData,V,problem);
            POD.ports{b}.CP=CP;
            POD.ports{b}.CN=CN;
        end
    end

    % Precompute reduced source vectors after the port coefficients.
   
    % This assumes source-region geometry is fixed across parameter samples.
    fprintf('[Source] Precomputing reduced source vectors.\n');
    POD.br=cell(size(POD.modes));
    for b=1:nBins
        wpoints=POD.wpoints{b};
        for d=1:size(POD.modes,2)
            POD.br{b,d}=complex(zeros(size(POD.modes{b,d},2),numel(wpoints)));
        end

        for wi=1:numel(wpoints)
            w=wpoints(wi);
            [phi,beta]=fdfd.buildSource(e,gridData,w,problem,sgmx,sgmy);
            Aw=buildReferenceOperator(gridData,w,e,sgmx,sgmy,problem);
            bFull=fdfd.applyQAAQ(gridData,Aw,phi,beta,problem);

            % Reuse the full source for every local basis in this bin.
            for d=1:size(POD.modes,2)
                V=POD.modes{b,d};
                POD.br{b,d}(:,wi)=V'*bFull;
            end
            fprintf('[Source] Bin %d: frequency %d/%d finished.\n',b,wi,numel(wpoints));
        end
    end

    fprintf('[Library] Saving %s ...\n',fullfile(rom.outputDir,'POD.mat'));
    save(fullfile(rom.outputDir,'POD.mat'),'POD','-v7.3'); % save the POD model
    fprintf('[Library] Saved successfully. Total elapsed: %.1f s.\n',toc(libraryTimer));
end
