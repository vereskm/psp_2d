function [U,lambda] = generatePODModes(snapShots,maxModes)
    %BUILDPODBASES builds the POD basis

    [U, sigma, ~] = svds([snapShots{:}],maxModes);
    lambda=diag(sigma).^2; % get the POD eigenvalues



end