function [phi,beta] = buildSource(e,dataGrid,w,problem,sgmx,sgmy)
    % solves for the sources (Does not do the QAAQ method just ye)
    assert(problem.polarization=="TE"||problem.polarization ...
        =="TM","Polarization needs to be either TE or TM");
    C=physics.constants();
    e0=C.e0;
    mu0=C.mu0;
    sources=problem.sources;
    Ns=numel(sources); % get the number of sources
    phi=cell(Ns,1);
    beta=zeros(Ns,1);
    if nargin < 6 || isempty(sgmx) || isempty(sgmy)
        [sgmx,sgmy] = fdfd.buildPML(problem,dataGrid,e);
    end
    x=dataGrid.x';
    y=dataGrid.y';
    X=dataGrid.X;
    Y=dataGrid.Y;
    for si=1:Ns
        sc=sources(si);
        p1=sc.p1;
        p2=sc.p2;
        np=sc.np;
        if sc.dir=="left" || sc.dir=="right"
            assert(abs(p1(1)-p2(1))<10E-16,"x components of source should match")
            xq=p1(1)*ones(np,1);
            yq=linspace(p1(2),p2(2),np)';
            dy=yq(2)-yq(1);
        
            eq=interp2(X,Y,e,xq,yq);
            sgmq=interp2(X,Y,sgmy,xq,yq); 
            if problem.polarization=="TE"
                    [phi{si},beta(si)] = ports.solveTEPort(eq,yq,w,sgmq);
                    P=1/2*beta(si)/(w)*sum(1./(e0.*eq).*phi{si}.*conj(phi{si}))*dy;
                    P=ports.validateModePower(P,"TE source "+si);
                    phi{si}=phi{si}/sqrt(P);
              
            else
                    [phi{si},beta(si)] = ports.solveTMPort(eq,yq,w,sgmq);
                    P=1/2*beta(si)/(w.*mu0)*sum(phi{si}.*conj(phi{si}))*dy;
                    P=ports.validateModePower(P,"TM source "+si);
                    phi{si}=phi{si}/sqrt(P);
               
            end
            phi{si}=interp1(yq,phi{si},y,"linear",0);
            
            phi{si}=reshape(phi{si},[],1);
        elseif sc.dir=="down" || sc.dir=="up"
              assert(abs(p1(2)-p2(2))<10E-16,"y components of source should match")
              yq=p1(2)*ones(np,1);
              xq=linspace(p1(1),p2(1),np)';
              dx=xq(2)-xq(1);
              eq=interp2(X,Y,e,xq,yq);
              sgmq=interp2(X,Y,sgmx,xq,yq);
              if problem.polarization=="TE"
                    [phi{si},beta(si)] = ports.solveTEPort(eq,xq,w,sgmq);
                    P=1/2*beta(si)/(w)*sum(1./(eq*e0).*phi{si}.*conj(phi{si}))*dx;
                    P=ports.validateModePower(P,"TE source "+si);
                    phi{si}=phi{si}/sqrt(P);
                  
              else
                    [phi{si},beta(si)] = ports.solveTMPort(eq,xq,w,sgmq);
                    P=1/2*beta(si)/(w.*mu0)*sum(phi{si}.*conj(phi{si}))*dx;
                    P=ports.validateModePower(P,"TM source "+si);
                    phi{si}=phi{si}/sqrt(P);
                    
              end
              phi{si}=interp1(xq,phi{si},x,"linear",0);
              phi{si}=reshape(phi{si},1,[]);
        else
            error("Source direction must left, right, down or up");
        end

    end
end

