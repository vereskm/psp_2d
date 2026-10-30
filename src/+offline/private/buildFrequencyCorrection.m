function [L] = buildFrequencyCorrection(gridData,w,wref,e,sgmx,sgmy,problem)
    %BUILDFREQUENCYCORRECTION: this function computes L
    C=physics.constants();
    c0=C.c0;
    Aw=buildReferenceOperator(gridData,w,e,sgmx,sgmy,problem);
    Aref=buildReferenceOperator(gridData,wref,e,sgmx,sgmy,problem);
    Q=buildMassMatrix(gridData,e,problem);
    L=Aw-Aref -(w.^2-wref.^2)/c0.^2*Q;
end