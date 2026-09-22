function C = constants()
%CONSTANTS Defines physical constants throught the program
    persistent phyConsts
    if isempty(phyConsts)
        phyConsts.e0=8.8541878188E-12; %
        phyConsts.mu0=1.256637061E-6;
        phyConsts.c0=1./sqrt(phyConsts.e0.*phyConsts.mu0);
    end
    C=phyConsts;
end