function samples = evaluateRBF(rbf,parameters)
    %EVALUATERBF Predict trained outputs for one parameter row.
    if iscell(parameters), parameters=cell2mat(parameters); end
    validateattributes(parameters,{'numeric'},{'real','finite'});
    parameters=reshape(parameters,1,[]);
    assert(numel(parameters)==rbf.numFeatures,'online:evaluateRBF:Features', ...
        'Parameters must follow the training feature order.');
    lower=rbf.bounds(1,:); upper=rbf.bounds(2,:);
    tolerance=32*eps(max(max(abs(rbf.bounds),[],1),realmin));
    assert(all(parameters>=lower-tolerance & parameters<=upper+tolerance), ...
        'online:evaluateRBF:OutsideBounds','Query lies outside this local model range.');
    parameters=min(max(parameters,lower),upper);
    if strcmp(rbf.kind,'constant')
        samples=rbf.value(:);
        return
    end
    active=rbf.activeFeatures;
    X=(parameters(active)-rbf.lower(active))./rbf.scale(active);
    distanceSquared=sum((rbf.centers-X).^2,2).';
    polynomial=1;
    if rbf.polynomialDegree==1, polynomial=[1 X]; end
    samples=(distanceSquared.^1.5*rbf.weights+polynomial*rbf.polynomial).';
    assert(all(isfinite(samples)),'online:evaluateRBF:Nonfinite', ...
        'RBF evaluation returned nonfinite output values.');
end
