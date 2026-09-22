function [b]=applyQAAQ(dataGrid,A,phi,beta,problem)
    %Does the QAAQ method Remember Q is positive in the scattering fiel
    Ns=numel(problem.sources); % get the number of sources

    Xin=dataGrid.X(2:end-1,2:end-1);
    Yin=dataGrid.Y(2:end-1,2:end-1); % get interior points
    sz=numel(Yin);
    b=complex(zeros(sz,1));
    for si=1:Ns
        sc=problem.sources(si);
        p1=sc.p1;
        cphi=phi{si};
        cbeta=beta(si);
        if sc.dir=="left" || sc.dir=="right"
            xl=p1(1);
           if sc.dir=="left"
                f=reshape(cphi(2:end-1).*exp(-1i.*cbeta.*(Xin-xl)),[],1);
                Q=reshape(Xin>xl,[],1);
           else
               f=reshape(cphi(2:end-1).*exp(1i.*cbeta.*(Xin-xl)),[],1);
               Q=reshape(Xin<xl,[],1);
           end
        elseif sc.dir=="down" || sc.dir=="up"
            yl=p1(2);
            if sc.dir=="down"
                f=reshape(cphi(2:end-1).*exp(-1i.*cbeta.*(Yin-yl)),[],1);
                Q=reshape(Yin>yl,[],1);
            else
                f=reshape(cphi(2:end-1).*exp(1i.*cbeta.*(Yin-yl)),[],1);
                Q=reshape(Yin<yl,[],1);
            end
            
        else
            
            error("Source direction must left, right, down or up");
        end
        b=b+(Q.*(A*f)-A*(Q.*f));
    end

end