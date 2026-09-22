function gridData = create(problem)
%CREATE the grid the problem is solved on
 % define the grid
    gridData.x=linspace(0,problem.domain.Lx,problem.mesh.Nx);
    gridData.y=linspace(0,problem.domain.Ly,problem.mesh.Ny);
    [gridData.X,gridData.Y]=meshgrid(gridData.x,gridData.y);

    % get the grid spaceing
    gridData.dx=gridData.x(2)-gridData.x(1);
    gridData.dy=gridData.y(2)-gridData.y(1);


    gridData.Nx=problem.mesh.Nx;
    gridData.Ny=problem.mesh.Ny;
    gridData.Lx=problem.domain.Lx;
    gridData.Ly=problem.domain.Ly;
    gridData.NxInterior = gridData.Nx - 2;
    gridData.NyInterior = gridData.Ny - 2;


    gridData.numberOfUnknowns=gridData.NxInterior.*gridData.NyInterior;

    gridData.interiorX = gridData.X(2:end-1,2:end-1);
    gridData.interiorY = gridData.Y(2:end-1,2:end-1);
end
