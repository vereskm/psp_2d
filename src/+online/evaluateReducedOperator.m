function [Ar,info] = evaluateReducedOperator(POD,w,parameters)
    %EVALUATEREDUCEDOPERATOR Reconstruct A + delta(k^2)*Q + deltaL.
    % Parameters follow POD.trainingParameters order. This evaluates the LHS
    % only; source interpolation and the final field solve are separate.
    % Interpolators return DEIM coefficients, so no projector is applied here.
    if iscell(parameters), parameters=cell2mat(parameters); end
    validateattributes(parameters,{'numeric'},{'real','finite'});
    validateattributes(w,{'numeric'},{'scalar','real','positive','finite'});
    parameters=reshape(parameters,1,[]);
    assert(numel(parameters)==numel(POD.trainingParameters), ...
        'online:evaluateReducedOperator:Features','Parameter count does not match the library.');
    assert(isfield(POD,'interpolationSpace') && ...
        strcmp(POD.interpolationSpace,'deimCoefficients'), ...
        'online:evaluateReducedOperator:CoefficientLibrary', ...
        'Rebuild this library to interpolate DEIM coefficients directly.');
    edges=POD.wBinEdges;
    b=find(w>=edges(1:end-1) & w<=edges(2:end),1,'first');
    assert(~isempty(b),'online:evaluateReducedOperator:Frequency', ...
        'Frequency lies outside the trained range.');
    d=1;
    if ~isempty(POD.localFeatureIndex)
        d=discretize(parameters(POD.localFeatureIndex),POD.localEdges);
        assert(~isnan(d),'online:evaluateReducedOperator:Region', ...
            'The local parameter lies outside the trained regions.');
    end
    assert(~isempty(POD.Arbf{b,d}),'online:evaluateReducedOperator:EmptyModel', ...
        'The selected local model has not been trained.');
    r=size(POD.modes{b,d},2);
    cA=online.evaluateRBF(POD.Arbf{b,d},parameters);
    cQ=online.evaluateRBF(POD.Qrbf{b,d},parameters);
    cL=interp1(POD.wpoints{b}(:),POD.DLcoefficients{b,d},w,'linear',NaN).';
    assert(all(isfinite(cL)),'online:evaluateReducedOperator:Correction', ...
        'Frequency correction interpolation returned nonfinite values.');
    A=reshape(POD.Abasis{b,d}*cA,r,r);
    Q=reshape(POD.Qbasis{b,d}*cQ,r,r);
    deltaL=reshape(POD.DLbasis{b,d}*cL,r,r);
    C=physics.constants();
    deltaK2=(w/C.c0)^2-(POD.referenceFrequency(b)/C.c0)^2;
    Ar=A+deltaK2*Q+deltaL;
    info.frequencyBin=b;
    info.localRegion=d;
    info.numModes=r;
end
