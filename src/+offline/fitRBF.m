function rbf = fitRBF(samplePoints,values,lower,upper,options)
    %FITRBF Fit output vectors using a cubic RBF and polynomial tail.
    % samplePoints: nSamples-by-nFeatures; values: nSamples-by-nOutputs.
    % Real parameter coordinates and complex output values are supported.
    if nargin<5, options=struct(); end
    validateattributes(samplePoints,{'numeric'},{'2d','real','finite'});
    validateattributes(values,{'numeric'},{'2d','nonempty','finite'});
    lower=reshape(lower,1,[]);
    upper=reshape(upper,1,[]);
    assert(size(values,1)==size(samplePoints,1) && ...
        size(samplePoints,2)==numel(lower) && numel(lower)==numel(upper), ...
        'offline:fitRBF:Dimensions','Sample rows and parameter bounds must match.');
    assert(all(isfinite([lower upper])) && all(upper>=lower), ...
        'offline:fitRBF:Bounds','Bounds must be finite and ordered.');
    assert(all(all(samplePoints>=lower & samplePoints<=upper)), ...
        'offline:fitRBF:OutsideBounds','Training points lie outside the declared bounds.');
    kernel=string(setting(options,'kernel','cubic'));
    degree=setting(options,'polynomialDegree',1);
    ridge=setting(options,'regularization',1e-11);
    maxCenters=setting(options,'maxCenters',size(samplePoints,1));
    normalize=setting(options,'normalize',true);
    assert(kernel=="cubic" && ismember(degree,[0 1]), ...
        'offline:fitRBF:Options','Use a cubic kernel with polynomial degree 0 or 1.');
    validateattributes(ridge,{'numeric'},{'scalar','real','finite','nonnegative'});
    validateattributes(maxCenters,{'numeric'},{'scalar','integer','positive','finite'});
    validateattributes(normalize,{'logical','numeric'},{'scalar','binary'});

    rbf.bounds=[lower;upper];
    rbf.numFeatures=size(samplePoints,2);
    rbf.activeFeatures=find(upper>lower);
    rbf.lower=lower;
    rbf.scale=upper-lower;
    rbf.scale(rbf.scale==0)=1;
    if ~normalize
        rbf.lower=zeros(size(lower));
        rbf.scale=ones(size(lower));
    end
    rbf.samplePoints=samplePoints;
    rbf.trainingValues=values;
    variation=norm(values-values(1,:),'fro');
    if variation<=32*eps*max(norm(values,'fro'),realmin)
        rbf.kind='constant';
        rbf.value=values(1,:);
        rbf.numCenters=0;
        return
    end

    [points,first,group]=unique(samplePoints,'rows');
    uniqueValues=values(first,:);
    assert(norm(values-uniqueValues(group,:),'fro')<= ...
        100*eps*max(norm(values,'fro'),realmin), ...
        'offline:fitRBF:DuplicateValues','Repeated parameter points have different values.');
    active=rbf.activeFeatures;
    X=(points(:,active)-rbf.lower(active))./rbf.scale(active);
    assert(~isempty(X),'offline:fitRBF:Coverage', ...
        'Varying outputs require varying parameter samples.');
    P=ones(size(X,1),1);
    if degree==1, P=[P X]; end
    assert(rank(P)==size(P,2),'offline:fitRBF:Coverage', ...
        'Samples must span the polynomial tail. Add more independent geometries.');
    count=min(maxCenters,size(X,1));
    assert(count>=size(P,2),'offline:fitRBF:Centers', ...
        'maxCenters is too small for the polynomial tail.');
    [~,~,pivot]=qr(P.','vector');
    selected=pivot(1:size(P,2));
    distance=inf(size(X,1),1);
    for k=selected
        distance=min(distance,sum((X-X(k,:)).^2,2));
    end
    distance(selected)=-inf;
    while numel(selected)<count
        [~,next]=max(distance);
        selected(end+1)=next; %#ok<AGROW>
        distance=min(distance,sum((X-X(next,:)).^2,2));
        distance(selected)=-inf;
    end
    centers=X(selected,:);
    P=P(selected,:);
    distanceSquared=max(0,sum(centers.^2,2)+sum(centers.^2,2).'-2*(centers*centers.'));
    Phi=distanceSquared.^1.5;
    ridge=ridge*max(1,norm(Phi,'fro')/count);
    system=[Phi+ridge*eye(count),P; P.',zeros(size(P,2))];
    assert(rcond(system)>eps,'offline:fitRBF:Condition', ...
        'The RBF system is numerically singular. Check sample coverage and scaling.');
    solution=system\[uniqueValues(selected,:);zeros(size(P,2),size(values,2))];
    rbf.kind='cubic';
    rbf.centers=centers;
    rbf.weights=solution(1:count,:);
    rbf.polynomial=solution(count+1:end,:);
    rbf.polynomialDegree=degree;
    rbf.numCenters=count;
end

function value=setting(options,name,default)
    if isfield(options,name), value=options.(name); else, value=default; end
end
