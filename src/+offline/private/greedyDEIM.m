function [U,projM,pi] = greedyDEIM(Mats,maxModes,energytol)
    %GREEDYDEIM Summary of this function goes here
    columns = cellfun(@(M) M(:),Mats(:), ...
    'UniformOutput',false);
    S = [columns{:}];
    totalEnergy=norm(S,'fro').^2;
    k = min([maxModes,size(S)]);

    [U, sigma, ~] = svds(S,k);
    lambda=diag(sigma).^2;

    [lambda,order] = sort(lambda,'descend');
    U = U(:,order);

    retainedEnergy = cumsum(lambda)/totalEnergy;
    nModes = find(retainedEnergy >= 1-energytol,1,'first');

    if isempty(nModes)
        nModes=k;
        warning("Mode limit reached before energy tolerance.");
    end
    U=U(:,1:nModes);
    pi = zeros(nModes,1);
    [~,pi(1)] = max(abs(U(:,1)));

    for mI = 2:nModes
        previous = 1:mI-1;
        Uold = U(:,previous);

        c = Uold(pi(previous),:) \ U(pi(previous),mI);
        residual = U(:,mI) - Uold*c;

        scores = abs(residual);
        scores(pi(previous)) = -Inf;
        [~,pi(mI)] = max(scores);
    end

    projM = U(pi,:) \ eye(nModes);

end