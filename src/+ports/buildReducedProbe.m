function [LibraryCP,LibraryCN] = buildReducedProbe(e,sgmx,sgmy,w,gridData,Modes,problem)
    % library of probs to solve POD quickly in
    nModes = size(Modes,2);
    nProbes = numel(problem.probes);

    LibraryCP = complex(zeros(nModes,numel(w),nProbes));
    LibraryCN = complex(zeros(nModes,numel(w),nProbes));


    fprintf('[Ports] %d frequencies x %d modes x %d probes.\n',numel(w),nModes,nProbes);
    parfor i=1:length(w)
        for j=1:nModes
            [cp,cn]=ports.SolveProb(e,sgmx,sgmy,w(i),gridData,Modes(:,j),problem);
            LibraryCP(j,i,:) = reshape(cp,1,1,nProbes);
            LibraryCN(j,i,:) = reshape(cn,1,1,nProbes);

        end
        fprintf('[Ports] Frequency %d/%d finished.\n',i,numel(w));
    end



end

