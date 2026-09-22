function [snapshots]=generateSnapshots(problem,training)
    % GENERATESNAPSHOTS: Generate training snapshots
    %output: snapshots
    % snapshots.wBinEdges: edges of the frequency bins
    % snapshots.wSamples: Frequencies used during training
    %snapshots.parameterSamples: each row corisponds to a training sample
    % each column is the feature values, each depth (third index) is the
    % frequency bin
    % snapshots.fields: outputed fields (nTrainxbins) THIS IS A CELL ARRAY
    C = physics.constants();
    rng(training.randomSeed); % set the random seed
    pol=problem.polarization; % get the polarization 
    % Create the mesh
    gridData = mesh.create(problem);
    % get the number of bins and frequencies we need:
    
    % get the number of parameters using during training
    nFeat=numel(training.parameters);
    
    % get the number of geometries to get snapshots for
    nGeometry=training.geometriesPerBin; 
    
    % get the frequencies and bins we are using
    numBin=training.frequencyBins;
    nFrequency=training.frequenciesPerBin;
    nTrain=nGeometry*nFrequency;

    % Convert the wavelength limits to angular frequency. The reversed
    % wavelength limits keep the angular-frequency edges increasing.
    lambdaMin = min(training.wavelengthRange);
    lambdaMax = max(training.wavelengthRange);
    wmin = 2*pi*C.c0/lambdaMax;
    wmax = 2*pi*C.c0/lambdaMin;
    snapshots.wBinEdges=linspace(wmin,wmax,numBin+1);
    snapshots.wSamples = zeros(nTrain,numBin);
     

   
    % get the features
    snapshots.parameterSamples = zeros(nTrain, ...
        nFeat,numBin);

     snapshots.fields=cell(nTrain,numBin);
  
   
    
    % Latin-hypercube sampling: every feature uses each of nGeometry
    % equal-width strata exactly once.
    geometryUnit=lhsdesign(nGeometry,nFeat); 
    
    % each row coresponds to a geometry: columns are features
    geometrySamples = zeros(nGeometry,nFeat);
    for i = 1:nFeat
        low = training.parameters(i).range(1);
        high = training.parameters(i).range(2);
        geometrySamples(:,i) = low + (high-low)*geometryUnit(:,i);
    end

    for b = 1:numBin
        wl = snapshots.wBinEdges(b);
        wh = snapshots.wBinEdges(b+1);
        frequencyValues = linspace(wl,wh,nFrequency).';
        
        % each row corisponds to a particular geometry at a particular
        % frequency
        snapshots.parameterSamples(:,:,b) = repelem( ...
            geometrySamples,nFrequency,1);
        snapshots.wSamples(:,b) = repmat(frequencyValues,nGeometry,1);
    end
   
    
% Run training solves
for b = 1:numBin
    fieldsBins=cell(nTrain,1)
    % Copy bin data into sliced/local variables
    wSamplesB = snapshots.wSamples(:,b);
    parameterSamplesB = snapshots.parameterSamples(:,:,b);

    parfor gi = 1:nTrain
        problemRef = problem;
        Gref = problem.geometry;

        for i = 1:nFeat
            d = training.parameters(i).device;
            fieldname = training.parameters(i).field;

            Gref(d).(fieldname) = parameterSamplesB(gi,i);
        end
        problemRef.geometry = Gref;

        wSample = wSamplesB(gi);

        e = geometry.generatePermittivity(problemRef,gridData);
        [sgmx,sgmy] = fdfd.buildPML(problemRef,gridData,e);

        fprintf( ...
            "Solving training structure %d/%d for bin %d/%d\n", ...
            gi,nTrain,b,numBin);

        if pol == "TE"
            [phiTE,betaTE] = fdfd.buildSource( ...
                e,gridData,wSample,problemRef,sgmx,sgmy);

            ATE = fdfd.assembleTE(gridData,wSample,e,sgmx,sgmy);

            bTE = fdfd.applyQAAQ( ...
                gridData,ATE,phiTE,betaTE,problemRef);

            fieldsBins{gi} = ATE\bTE;
        end

        if pol == "TM"
            [phiTM,betaTM] = fdfd.buildSource( ...
                e,gridData,wSample,problemRef,sgmx,sgmy);

            ATM = fdfd.assembleTM(gridData,wSample,e,sgmx,sgmy);

            bTM = fdfd.applyQAAQ( ...
                gridData,ATM,phiTM,betaTM,problemRef);

            fieldsBins{gi} = ATM\bTM;
        end
    end

    % Store completed bin after parfor
    snapshots.fields(:,b) = fieldsBins;
    
   
   
end
    
    % save the result
    saveSnapshots(problem,training,snapshots);


end

