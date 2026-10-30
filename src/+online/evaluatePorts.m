function [cPlus,cMinus] = evaluatePorts(POD,bin,localRegion,w,a)
    %EVALUATEPORTS Interpolate saved probe maps and apply reduced coefficients.
    % The operation is linear in a; use no conjugation and only its retained
    % modes. No full field or online port eigensolve is required.
    validateattributes(w,{'numeric'},{'scalar','real','finite'});
    validateattributes(a,{'numeric'},{'column','nonempty','finite'});
    assert(isfield(POD,'ports') && ~isempty(POD.ports{bin,localRegion}), ...
        'online:evaluatePorts:Library','Precompute and save the port maps before evaluating them.');
    probe=POD.ports{bin,localRegion};
    grid=POD.wpoints{bin}(:);
    assert(numel(grid)>=2 && all(isfinite(grid)) && all(diff(grid)>0), ...
        'online:evaluatePorts:Grid','Port frequencies must be increasing.');
    assert(w>=grid(1) && w<=grid(end),'online:evaluatePorts:Frequency', ...
        'Frequency is outside the precomputed port range.');
    assert(isequal(size(probe.CP),size(probe.CN)) && size(probe.CP,2)==numel(grid), ...
        'online:evaluatePorts:Dimensions','Port arrays do not match their frequency grid.');
    nModes=numel(a);
    assert(nModes<=size(probe.CP,1),'online:evaluatePorts:Modes', ...
        'More modes were requested than are present in the saved port maps.');
    nProbes=size(probe.CP,3);
    cPlus=complex(zeros(nProbes,1));
    cMinus=complex(zeros(nProbes,1));
    for port=1:nProbes
        cpGrid=reshape(probe.CP(1:nModes,:,port),nModes,[]).';
        cnGrid=reshape(probe.CN(1:nModes,:,port),nModes,[]).';
        cp=interp1(grid,cpGrid,w,'linear',NaN);
        cn=interp1(grid,cnGrid,w,'linear',NaN);
        cPlus(port)=cp*a;
        cMinus(port)=cn*a;
    end
    assert(all(isfinite(cPlus)) && all(isfinite(cMinus)), ...
        'online:evaluatePorts:Nonfinite','Port interpolation returned nonfinite amplitudes.');
end
