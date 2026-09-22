% This is a finite difference code for FDFD
function [A]=assembleTM(gridData,w,e,sgmx,sgmy)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % COEFICIENTS
    C = physics.constants();
    e0=C.e0;
    mu0=C.mu0;
    c0=C.c0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute the number of grid points
    Nx=gridData.Nx; % number of points including boundaries
    Ny=gridData.Ny; % number of points including boundaries
    % get interior points
    NxI = gridData.NxInterior;
    NyI = gridData.NyInterior;
    N   =  gridData.numberOfUnknowns;
    % Compute grid points
    dx=gridData.dx;
    dy=gridData.dy;
 
   
    k0=w/c0;
    % CREATE THE PML
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sx=1+1i.*sgmx./(w.*e0);
    sxinv=1./sx;
    sy=1+1i.*sgmy./(w.*e0);
    syinv=1./sy;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % create needed matricies for solver
    % get permativities at average points for both Y and X

    sxinvX=(sxinv(:,2:end)+sxinv(:,1:end-1))./2;
    sxinvXS=sxinvX/dx^2; % scale by the differential
    syinvY=(syinv(2:end,:)+syinv(1:end-1,:))./2;
    syinvYS=syinvY/dy^2;
    
    sxCenter=sxinv(2:end-1,2:end-1);
    syCenter = syinv(2:end-1,2:end-1);


    aW = sxCenter .* sxinvXS(2:end-1,1:NxI);
    aE = sxCenter .* sxinvXS(2:end-1,2:NxI+1);

    aS = syCenter .* syinvYS(1:NyI,2:end-1);
    aN = syCenter .* syinvYS(2:NyI+1,2:end-1);
    eInterior = e(2:end-1,2:end-1);
    aC = -(aW + aE + aS + aN) + k0^2*eInterior;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Column-major global indices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G = reshape(1:N,NyI,NxI);

hasW = false(NyI,NxI);
hasE = false(NyI,NxI);
hasS = false(NyI,NxI);
hasN = false(NyI,NxI);

hasW(:,2:end)   = true;
hasE(:,1:end-1) = true;
hasS(2:end,:)   = true;
hasN(1:end-1,:) = true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assemble sparse triplet arrays
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assemble sparse triplet arrays
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rows = [ ...
    G(:);
    G(hasW);
    G(hasE);
    G(hasS);
    G(hasN)];

cols = [ ...
    G(:);
    G(hasW)-NyI;
    G(hasE)+NyI;
    G(hasS)-1;
    G(hasN)+1];

data = [ ...
    aC(:);
    aW(hasW);
    aE(hasE);
    aS(hasS);
    aN(hasN)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = sparse(rows,cols,data,N,N); 
end
