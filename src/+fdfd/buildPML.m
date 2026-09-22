function [sgmx,sgmy] = buildPML(problem,dataGrid,e)
%BUILDPML: Creates the pml 
%
% Inputs:
% problem: settings concerning the problem you are solving
% dataGrid: information regarding the mesh 
% e relative permativity
% Outpus: 
% sgmx, sgmy: conductivy profiles in x and y directions

Nx=dataGrid.Nx;
Ny=dataGrid.Ny;
Lx=dataGrid.Lx;
Ly=dataGrid.Ly;
dx=dataGrid.dx;
dy=dataGrid.dy;
X=dataGrid.X;
Y=dataGrid.Y;
d=problem.pml.thickness;
R=problem.pml.reflection;
m=problem.pml.order;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the size of the PML ( this needs to be an integer number of grid
% points)
NxD=ceil(d/dx);
dax=NxD*dx;
NyD=ceil(d/dy);
day=NyD*dy;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get parameters needed for PML
% Have the PML line exactly on points
C=physics.constants();
e0=C.e0;
mu0=C.mu0;

eta0=sqrt(mu0/e0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sgmx=zeros(size(X)); % initialize the boundaries (zeros)
sgmy=zeros(size(Y));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE LEFT SIDE sgmx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% put the PMl on an actual grid point
Yq=linspace(0,Ly,Ny)';
XqLeft=dax*ones(Ny,1);
eLeft = interp2(X,Y,e,XqLeft,Yq);
etaL=eta0./sqrt(eLeft);
sigMaxL=-(m+1)*log(R)./(2*etaL*dax);
sgmx(:,1:NxD)=sigMaxL.*(abs(dax-X(:,1:NxD))./dax).^m;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE RIGHT SIDE sgmx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
XqRight=(Lx-dax)*ones(Ny,1);
eRight=interp2(X,Y,e,XqRight,Yq);
etaR=eta0./sqrt(eRight);
sigMaxR=-(m+1)*log(R)./(2*etaR*dax);
sgmx(:,end-(NxD-1):end)=sigMaxR.*(abs(X(:,end-(NxD-1):end)-(Lx-dax))./dax).^m;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% determine sgmy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE Bottom SIDE of sgmy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Xq=linspace(0,Lx,Nx);
YqBottom=(day)*ones(1,Nx);
eBottom=interp2(X,Y,e,Xq,YqBottom);
etaB=eta0./sqrt(eBottom);
sigMaxB=-(m+1)*log(R)./(2*etaB*day);
sgmy(1:NyD,:)=sigMaxB.*(abs(day-Y(1:NyD,:))./day).^m;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE TOP SIDE of sgmy
Xq=linspace(0,Lx,Nx);
YqTop=(Ly-day)*ones(1,Nx);
eTop=interp2(X,Y,e,Xq,YqTop);
etaT=eta0./sqrt(eTop);
sigMaxT=-(m+1)*log(R)./(2*etaT*day);
sgmy(end-(NyD-1):end,:)=sigMaxT.*(abs(Y(end-(NyD-1):end,:)-(Ly-day))./day).^m;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

