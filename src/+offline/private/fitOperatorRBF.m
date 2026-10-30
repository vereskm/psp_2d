function rbf = fitOperatorRBF(matrices,indices,projM,geometrySamples,lower,upper,rom)
    %FITOPERATORRBF Fit geometry directly to DEIM coefficients.
    columns=cellfun(@(M) M(:),matrices(:),'UniformOutput',false);
    snapshots=[columns{:}];
    coefficients=(projM*snapshots(indices,:)).';
    options=struct();
    if isfield(rom,'rbf'), options=rom.rbf; end
    rbf=offline.fitRBF(geometrySamples,coefficients,lower,upper,options);
    rbf.outputSpace='deimCoefficients';
end
