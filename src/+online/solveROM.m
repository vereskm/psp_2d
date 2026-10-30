function [a,info,u] = solveROM(POD,w,parameters,nModes)
    % Port coefficients are returned in info.cPlus and info.cMinus.
    % Request only [a,info] to skip full-field reconstruction.
    % start the timer
    assemblytimer = tic;

    %Interprete the reduced operator
    [Ar,info] = online.evaluateReducedOperator(POD,w,parameters);
    
    b = info.frequencyBin;
    d = info.localRegion;
    V = POD.modes{b,d};

    % Use all available modes unless a count is supplied.
    if nargin < 4 || isempty(nModes)
        nModes = size(V,2);
    end

    
    validateattributes(nModes,{'numeric'}, ...
        {'scalar','integer','positive','<=',size(V,2)});

    % Interpolate the reduced source at this frequency.
    br = interp1(POD.wpoints{b}(:), ...
        POD.br{b,d}.',w,'linear').';

    

    % Restrict to the requested POD modes.
    Ar = Ar(1:nModes,1:nModes);
    br = br(1:nModes);
    info.assemblyTime=toc(assemblytimer);
    % Solve for the reduced field coefficients.
    solveTimer=tic;
    a = Ar \ br;
    info.solveTime=toc(solveTimer);

    info.numModes = nModes;

    portTimer = tic;
    [info.cPlus,info.cMinus] = online.evaluatePorts(POD,b,d,w,a);
    info.portTime = toc(portTimer);
    info.reconstructionTime = 0; % this is if the user does not do this

    % Reconstruct the interior field only if requested.
    if nargout >= 3
        reconstructionTimer = tic;
        u = V(:,1:nModes)*a;
        info.reconstructionTime = toc(reconstructionTimer);
    end
    info.totalTime = info.assemblyTime+info.solveTime + info.portTime + info.reconstructionTime;
end