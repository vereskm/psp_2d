function [LibraryCP,LibraryCN] = ProbLibPOD(e,sgmx,sgmy,w,X,Y,Modes,Probes,settings,pol)
    % library of probs to solve POD quickly in
    nModes = size(Modes,2);
    nProbes = numel(Probes);

    LibraryCP = complex(zeros(nModes,numel(w),nProbes));
    LibraryCN = complex(zeros(nModes,numel(w),nProbes));

   
    for i=1:length(w)
        for j=1:nModes
            [cp,cn]=SolveProb(e,sgmx,sgmy,w(i),X,Y,Modes(:,j),Probes,settings,pol);
            LibraryCP(j,i,:) = reshape(cp,1,1,nProbes);
            LibraryCN(j,i,:) = reshape(cn,1,1,nProbes);

        end
    end
   
    

end

