function [cPlus,cMinus] = SolveProb(e,sgmx,sgmy,w,X,Y,fieldI,Probes,settings,pol)
    % solves prob coefficients
    assert(pol=="TE"||pol=="TM","Polarization needs to be either TE or TM");
    e0=8.8541878188E-12;
    mu0=1.256637061E-6;
    Lx=settings.Lx;
    Ly=settings.Ly;
    Nx=settings.Nx;
    Ny=settings.Ny;
    x=linspace(0,Lx,Nx)';
    y=linspace(0,Ly,Ny)';
    Np=numel(Probes); % get the total number of probes
    cPlus=zeros(Np,1); % get positive coefficients
    cMinus=zeros(Np,1); % get negative coefficients
    % field onl contain interior points so they need to be exteded
    field=complex(zeros(Ny,Nx));
    field(2:end-1,2:end-1)=reshape(fieldI,Ny-2,Nx-2);

    if pol=="TE"
         [Hzq_dx,Hzq_dy]=gradient(field,x,y);
    elseif pol=="TM"
        [Ezq_dx,Ezq_dy]=gradient(field,x,y);
    end
    for pi=1:Np
        pc=Probes(pi);
        p1=pc.p1;
        p2=pc.p2;
        np=pc.np;
        if pc.dir=="left" || pc.dir=="right"
            assert(abs(p1(1)-p2(1))<10E-16,"x components of source should match")
            xq=p1(1)*ones(np,1);
            yq=linspace(p1(2),p2(2),np)';
            dy=yq(2)-yq(1);
            eq=interp2(X,Y,e,xq,yq);
            sgmq=interp2(X,Y,sgmy,xq,yq);
            if pol=="TE"
                    [Hmz,beta] =TE_Source(eq,yq,w,sgmq);
                    P=1/2*beta/(w)*sum(1./(e0.*eq).*Hmz.*conj(Hmz))*dy;
                    P=validateModePower(P,"TE probe "+pi);
                    Hmz=Hmz/sqrt(P);
                    Hmz=reshape(Hmz,[],1);
                    Hm=[zeros(np,1),zeros(np,1),Hmz];
                    if pc.dir=="left"
                        Em=[1i./(w.*eq.*e0).*gradient(Hmz,dy),-beta./(w.*eq.*e0).*Hmz,zeros(np,1)];
                        pHat=[-1,0,0];
                    elseif pc.dir=="right"
                        Em=[1i./(w.*eq.*e0).*gradient(Hmz,dy),beta./(w.*eq.*e0).*Hmz,zeros(np,1)];
                        pHat=[1,0,0];
                    end
                    Hzq=reshape(interp2(X,Y,field,xq,yq),[],1);
                    %[Hzq_dx,Hzq_dy]=gradient(field,x,y); % TODO don't evalaute the gradient on the entire grid
                    Ex=1i./(w.*eq.*e0).*reshape(interp2(X,Y,Hzq_dy,xq,yq),[],1);
                    Ey=-1i./(w.*eq.*e0).*reshape(interp2(X,Y,Hzq_dx,xq,yq),[],1);
                    H=[zeros(np,1),zeros(np,1),Hzq];
                    E=[Ex,Ey,zeros(np,1)];
                    [cPlus(pi),cMinus(pi)] = ModeOverlap(yq,E,H,Em,Hm,pHat);
            else % TM
                    [Emz,beta] =TM_Source(eq,yq,w,sgmq);
                    P=1/2*beta/(w.*mu0)*sum(Emz.*conj(Emz))*dy;
                    P=validateModePower(P,"TM probe "+pi);
                    Emz=Emz/sqrt(P);
                    Emz=reshape(Emz,[],1);
                    Em=[zeros(np,1),zeros(np,1),Emz];
                    if pc.dir=="left"
                        Hm=[-1i./(w.*mu0).*gradient(Emz,dy),beta./(w.*mu0).*Emz,zeros(np,1)];
                        pHat=[-1,0,0];
                    elseif pc.dir=="right"
                        Hm=[-1i./(w.*mu0).*gradient(Emz,dy),-beta./(w.*mu0).*Emz,zeros(np,1)];
                        pHat=[1,0,0];
                    end
                    Ezq=reshape(interp2(X,Y,field,xq,yq),[],1);
                    %[Ezq_dx,Ezq_dy]=gradient(field,x,y); % TODO don't evalaute the gradient on the entire grid
                    Hx=-1i./(w.*mu0).*reshape(interp2(X,Y,Ezq_dy,xq,yq),[],1);
                    Hy=1i./(w.*mu0).*reshape(interp2(X,Y,Ezq_dx,xq,yq),[],1);
                    E=[zeros(np,1),zeros(np,1),Ezq];
                    H=[Hx,Hy,zeros(np,1)];
                    [cPlus(pi),cMinus(pi)] = ModeOverlap(yq,E,H,Em,Hm,pHat);
            end
         
        elseif pc.dir=="down" || pc.dir=="up"
              assert(abs(p1(2)-p2(2))<10E-16,"y components of source should match")
              yq=p1(2)*ones(np,1);
              xq=linspace(p1(1),p2(1),np)';
              dx=xq(2)-xq(1);
              eq=interp2(X,Y,e,xq,yq);
              sgmq=interp2(X,Y,sgmx,xq,yq);
              if pol=="TE"
                [Hmz,beta] =TE_Source(eq,xq,w,sgmq);
                P=1/2*beta/(w)*sum(1./(eq*e0).*Hmz.*conj(Hmz))*dx;
                P=validateModePower(P,"TE probe "+pi);
                Hmz=reshape(Hmz/sqrt(P),[],1);
                Hm=[zeros(np,1),zeros(np,1),Hmz];
                if pc.dir=="down"
                     Em=[beta/(w.*eq.*e0).*Hmz,-1i/(w.*eq.*e0).*gradient(Hmz,xq),zeros(np,1)];
                    pHat=[0,-1,0];
                elseif pc.dir=="up"
                    Em=[-beta/(w.*eq.*e0).*Hmz,-1i/(w.*eq.*e0).*gradient(Hmz,xq),zeros(np,1)];
                    pHat=[0,1,0];
                end
                 Hzq=reshape(interp2(X,Y,field,xq,yq),[],1);
                 %[Hzq_dx,Hzq_dy]=gradient(field,x,y); % TODO don't evalaute the gradient on the entire grid
                 Ex=1i./(w.*eq.*e0).*reshape(interp2(X,Y,Hzq_dy,xq,yq),[],1);
                 Ey=-1i./(w.*eq.*e0).*reshape(interp2(X,Y,Hzq_dx,xq,yq),[],1);
                 H=[zeros(np,1),zeros(np,1),Hzq];
                 E=[Ex,Ey,zeros(np,1)];
                 [cPlus(pi),cMinus(pi)] = ModeOverlap(xq,E,H,Em,Hm,pHat);
                  
              else
                [Emz,beta] =TM_Source(eq,xq,w,sgmq);
                P=1/2*beta/(w.*mu0)*sum(Emz.*conj(Emz))*dx;
                P=validateModePower(P,"TM probe "+pi);
                Emz=reshape(Emz/sqrt(P),[],1);
                Em=[zeros(np,1),zeros(np,1),Emz];
                if pc.dir=="down"
                     Hm=[-beta/(w.*mu0).*Emz,1i/(w.*mu0).*gradient(Emz,xq),zeros(np,1)];
                    pHat=[0,-1,0];
                elseif pc.dir=="up"
                    Hm=[beta/(w.*mu0).*Emz,1i/(w.*mu0).*gradient(Emz,xq),zeros(np,1)];
                    pHat=[0,1,0];
                end
                Ezq=reshape(interp2(X,Y,field,xq,yq),[],1);
                %[Ezq_dx,Ezq_dy]=gradient(field,x,y); % TODO don't evalaute the gradient on the entire grid
                Hx=-1i./(w.*mu0).*reshape(interp2(X,Y,Ezq_dy,xq,yq),[],1);
                Hy=1i./(w.*mu0).*reshape(interp2(X,Y,Ezq_dx,xq,yq),[],1);
                E=[zeros(np,1),zeros(np,1),Ezq];
                H=[Hx,Hy,zeros(np,1)];
                [cPlus(pi),cMinus(pi)] = ModeOverlap(xq,E,H,Em,Hm,pHat);
              end
           
        else
            error("prob direction must left, right, down or up");
        end

    end


end

