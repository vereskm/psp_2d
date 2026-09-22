function [eSmooth] = generatePermittivity(problem,gridData)
%GENERATEPERMITTIVITY generate a mesh for the permativity

nSub=problem.mesh.subpixelSamples;
Nx=gridData.Nx;
Ny=gridData.Ny;
dx=gridData.dx;
dy=gridData.dy;
X=gridData.X;
Y=gridData.Y;
eSmooth = zeros(Ny,Nx);

% Subpixel locations centered within the control volume
ox = ((1:nSub)-0.5)/nSub - 0.5;
oy = ((1:nSub)-0.5)/nSub - 0.5;

for jy = 1:nSub
    for ix = 1:nSub
        Xs = X + ox(ix)*dx;
        Ys = Y + oy(jy)*dy;
        eSample=problem.material.backgroundEr*ones(Ny,Nx);
        for n=1:numel(problem.geometry) % loop through devices
            isInside=geometry.inside(problem.geometry(n),Xs,Ys);
            eSample(isInside)=problem.geometry(n).er;
        end
        eSmooth=eSmooth+eSample;
    end
end
eSmooth = eSmooth/nSub^2;
end

