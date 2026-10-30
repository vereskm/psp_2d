function [A] = buildReferenceOperator(gridData,w,e,sgmx,sgmy,problem)
    %BUILDREFERENCEOPERATOR: this function will build the A matrix as part of
    %the MDEIM.
    if problem.polarization=="TE"
        A=fdfd.assembleTE(gridData,w,e,sgmx,sgmy);
    elseif problem.polarization=="TM"
        A=fdfd.assembleTM(gridData,w,e,sgmx,sgmy);
    end
end