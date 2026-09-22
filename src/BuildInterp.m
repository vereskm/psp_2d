function [FAr,Fbr] = BuildInterp(Ar,br,b,rangesM,wTp,wBinEdges)

    % Geometry parameter interpolation grids.
    % rangesM is a cell array: {p1Grid, p2Grid, ..., pNGrid}.
    % Ar/br are stored with a single linear geometry index gi produced by
    % ndgrid(rangesM{:}) in Train.m.
    nfeat = numel(rangesM);
    paramSizes = reshape(cellfun(@numel,rangesM),1,[]);
    interpGrid = [rangesM(:).' {[]}];

    % Frequencies belonging to this POD bin
    Iw = find(wTp >= wBinEdges(b) & ...
              wTp <= wBinEdges(b+1));

    wBin = wTp(Iw);
    interpGrid{end} = wBin;

    assert(numel(wBin) >= 2 && ...
           abs(wBin(1)-wBinEdges(b)) <= 10*eps(wBinEdges(b)) && ...
           abs(wBin(end)-wBinEdges(b+1)) <= 10*eps(wBinEdges(b+1)), ...
        "Frequency samples for POD bin %d must span both bin edges.",b);

    % Number of parameter/frequency samples
    np = prod(paramSizes);
    nw = numel(Iw);

    % Reduced dimension
    nr = size(Ar{b,1,Iw(1)},1);

    % One interpolant for each reduced-matrix/vector entry
    FAr = cell(nr,nr);
    Fbr = cell(nr,1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Reduced matrix interpolation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for r = 1:nr
        for c = 1:nr

            % Values of A_r(r,c) over parameter grid x frequency
            V = complex(zeros([paramSizes,nw]));

            for pi = 1:np
                paramSub = cell(1,nfeat);
                if nfeat == 1
                    paramSub{1} = pi;
                else
                    [paramSub{:}] = ind2sub(paramSizes,pi);
                end
                for jj = 1:nw

                    wi = Iw(jj);

                    Ared = Ar{b,pi,wi};

                    V(paramSub{:},jj) = Ared(r,c);
                end
            end

            FAr{r,c} = griddedInterpolant( ...
                interpGrid, ...
                V, ...
                'linear', ...
                'none');
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Reduced RHS interpolation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for r = 1:nr

        V = complex(zeros([paramSizes,nw]));

        for pi = 1:np
            paramSub = cell(1,nfeat);
            if nfeat == 1
                paramSub{1} = pi;
            else
                [paramSub{:}] = ind2sub(paramSizes,pi);
            end
            for jj = 1:nw

                wi = Iw(jj);

                bred = br{b,pi,wi};

                V(paramSub{:},jj) = bred(r);
            end
        end

        Fbr{r} = griddedInterpolant( ...
            interpGrid, ...
            V, ...
            'linear', ...
            'none');
    end

end
