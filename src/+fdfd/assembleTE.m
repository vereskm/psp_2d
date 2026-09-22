function A = assembleTE(gridData,w,e,sgmx,sgmy)
%TE_FDFD Assemble the TE FDFD operator using a five-point stencil.
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
 
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PML stretching
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sxinv = 1 ./ (1 + 1i*sgmx/(w*e0));
syinv = 1 ./ (1 + 1i*sgmy/(w*e0));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Face-centered coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
einv = 1./e;

% Quantities that appear inside the spatial derivatives
einvSx = einv .* sxinv;
einvSy = einv .* syinv;

% Average to the x- and y-directed faces and include grid scaling
einvX = (einvSx(:,1:end-1) + einvSx(:,2:end))/(2*dx^2);
einvY = (einvSy(1:end-1,:) + einvSy(2:end,:))/(2*dy^2);

% Outer PML factors evaluated at interior grid points
sxCenter = sxinv(2:end-1,2:end-1);
syCenter = syinv(2:end-1,2:end-1);

% Five-point stencil coefficients
aW = sxCenter .* einvX(2:end-1,1:NxI);
aE = sxCenter .* einvX(2:end-1,2:NxI+1);

aS = syCenter .* einvY(1:NyI,2:end-1);
aN = syCenter .* einvY(2:NyI+1,2:end-1);

aC = -(aW + aE + aS + aN) + k0^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Column-major global indices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G = reshape(1:N,NyI,NxI);

% Points having the indicated neighbor
gW = G(:,2:end);
gE = G(:,1:end-1);
gS = G(2:end,:);
gN = G(1:end-1,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assemble sparse matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rows = [
    G(:);
    gW(:);
    gE(:);
    gS(:);
    gN(:)];

cols = [
    G(:);
    gW(:)-NyI;
    gE(:)+NyI;
    gS(:)-1;
    gN(:)+1];

data = [
    aC(:);
    reshape(aW(:,2:end),[],1);
    reshape(aE(:,1:end-1),[],1);
    reshape(aS(2:end,:),[],1);
    reshape(aN(1:end-1,:),[],1)];

A = sparse(rows,cols,data,N,N);

end
