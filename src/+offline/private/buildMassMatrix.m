function [Q] = buildMassMatrix(gridData,e,problem)
    %BUILDMASSMATRIX: This function will build the mass matrix
    N=gridData.numberOfUnknowns;
    if problem.polarization=="TE"
        Q=speye(N);
    elseif problem.polarization=="TM"
        eint=e(2:end-1,2:end-1);
        Q=spdiags(eint(:),0,N,N);
    else
        error("Polarization must be TE or TM");
    end
end