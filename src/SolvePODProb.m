function [cPlus,cMinus] = SolvePODProb(Data,bin,wq,a,pol)
%SOLVEPODPROB Evaluate precomputed probe overlaps for a POD solution.
% Each stored column is the probe response of one POD basis vector. Since
% the overlap operation is linear in the field, the response of V*a is the
% same linear combination of the stored modal responses.

    arguments
        Data
        bin (1,1) double {mustBeInteger,mustBePositive}
        wq (1,1) double
        a (:,1) double
        pol (1,1) string
    end

    assert(pol == "TE" || pol == "TM", ...
        'pol must be "TE" or "TM".');

    if pol == "TE"
        assert(isfield(Data,'ProbeTE'), ...
            "InterpTE.mat does not contain ProbeTE. Rebuild the library.");
        probe = Data.ProbeTE{bin};
    else
        assert(isfield(Data,'ProbeTM'), ...
            "InterpTM.mat does not contain ProbeTM. Rebuild the library.");
        probe = Data.ProbeTM{bin};
    end

    assert(isfield(Data,'wTp') && numel(Data.wTp) >= 2, ...
        "The interpolation library does not contain a valid wTp grid.");

    wGrid = Data.wTp(:);
    nModesAvailable = size(probe.CP,1);
    nProbes = size(probe.CP,3);
    nModesUsed = numel(a);
    assert(nModesUsed <= nModesAvailable, ...
        "Requested %d probe modes, but only %d are precomputed.", ...
        nModesUsed,nModesAvailable);
    assert(size(probe.CP,2) == numel(wGrid) && ...
        isequal(size(probe.CP),size(probe.CN)), ...
        "The saved probe-library dimensions are inconsistent with wTp.");
    assert(wq >= wGrid(1) && wq <= wGrid(end), ...
        "Probe frequency is outside the precomputed frequency range.");

    cPlus = complex(zeros(nProbes,1));
    cMinus = complex(zeros(nProbes,1));
    for pi = 1:nProbes
        cpGrid = reshape(probe.CP(:,:,pi),nModesAvailable,[]).';
        cnGrid = reshape(probe.CN(:,:,pi),nModesAvailable,[]).';
        cpAtFrequency = interp1(wGrid,cpGrid,wq,'linear');
        cnAtFrequency = interp1(wGrid,cnGrid,wq,'linear');
        cPlus(pi) = cpAtFrequency(1:nModesUsed)*a;
        cMinus(pi) = cnAtFrequency(1:nModesUsed)*a;
    end
end
