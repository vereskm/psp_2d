function [cPlusDNS,cMinusDNS,cPlusPOD,cMinusPOD,FieldLSE,ModeLSE] = Solve(settings,PData,G,w,Sources,Probes,um,options,nworkers)
    % umare the mode counts used for the reported solution.
    % maxM are the maximum mode counts used in the LSE sweep.
    %umte are the maximum number of MODES
    arguments
        settings
        PData
        G
        w (:,1) double
        Sources
        Probes
        um   (1,1) double {mustBeInteger,mustBePositive}
        options
        nworkers
    end



    if ~isfolder("./Results")
        mkdir("./Results");
    end
    pool = gcp("nocreate");
    if isempty(pool)
        pool = parpool("Processes",nworkers);
    elseif pool.NumWorkers ~= nworkers
        delete(pool);
        pool = parpool("Processes",nworkers);
    end


    e0=8.8541878188E-12;
    mu0=1.256637061E-6;
    c0=1./sqrt(e0.*mu0);

    % G is a cell of structures
    Lx=settings.Lx;
    Ly=settings.Ly;
    Nx=settings.Nx;
    Ny=settings.Ny;
    pol=settings.pol;
    x=linspace(0,Lx,Nx);
    y=linspace(0,Ly,Ny);
    [X,Y]=meshgrid(x,y);
    nP=numel(Probes); % get the total number of probes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create data arrays
    % Both DNS names must exist because the probe parfor references both
    % polarization branches during variable classification.
    OutputTEDNS=cell(1,numel(w));
    OutputTMDNS=cell(1,numel(w));
    aOutputTEPOD=[];
    aOutputTMPOD=[];
    OutputTEPOD=[]; % complete real-space TE fields
    OutputTMPOD=[]; % complete real-space TM fields
    
    MTimeTE=zeros(1,numel(w));
    DNSTESolveTime=zeros(1,numel(w));
    DNSProbeTimeTE=zeros(1,numel(w));
    PODSolveTimeTE=zeros(1,numel(w));
    PODRecoverTimeTE=zeros(1,numel(w));
    PODProbeTimeTE=zeros(1,numel(w));
    MTimeTM=zeros(1,numel(w));
    DNSTMSolveTime=zeros(1,numel(w));
    DNSProbeTimeTM=zeros(1,numel(w));
    PODSolveTimeTM=zeros(1,numel(w));
    PODRecoverTimeTM=zeros(1,numel(w));
    PODProbeTimeTM=zeros(1,numel(w));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve DNS first
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nw=numel(w);
    cPlusTE  = cell(1,numel(w));
    cMinusTE = cell(1,numel(w));
    cPlusTM  = cell(1,numel(w));
    cMinusTM = cell(1,numel(w));
    e=gen_e(settings,X,Y,G); % get the permativity
    [sgmx,sgmy]=PML(settings,X,Y,e);
    fprintf("Solve DNS\n");
    parfor wi=1:numel(w)
        fprintf("Solving frequency %d/%d\n",wi,nw);
        fprintf("Solving Freq %d out of %d\n",wi,numel(w));
        if pol=="TE" 
             tic
             [phi,beta]=MakeSource( ...
                 e,w(wi),X,Y,Sources,settings,"TE",sgmx,sgmy);
             ATE=TE_FDFD(settings,w(wi),e,sgmx,sgmy);
             bTE=QAAQ(X,Y,ATE,phi,beta,Sources);
             MTimeTE(wi)=toc;
             %SOLVE
             tic
             H=ATE\bTE;
             DNSTESolveTime(wi)=toc;
             OutputTEDNS{wi}=H;
        end
        if pol=="TM" 
             tic
             [phi,beta]=MakeSource( ...
                 e,w(wi),X,Y,Sources,settings,"TM",sgmx,sgmy);
             ATM=TM_FDFD(settings,w(wi),e,sgmx,sgmy);
             bTM=QAAQ(X,Y,ATM,phi,beta,Sources);
             MTimeTM(wi)=toc;
             % SOLVE
             tic
             E=ATM\bTM;
             DNSTMSolveTime(wi)=toc;
             OutputTMDNS{wi}=E;
             
        end
    end
    % Process DNS probes after all complete-field solves.
    parfor wi = 1:numel(w)
        if pol == "TE" 
            tic
            [cPlusTE{wi},cMinusTE{wi}] = SolveProb( ...
                e,sgmx,sgmy,w(wi),X,Y,OutputTEDNS{wi},Probes,settings,"TE");
            DNSProbeTimeTE(wi)=toc;
        end
        if pol == "TM" 
            tic
            [cPlusTM{wi},cMinusTM{wi}] = SolveProb( ...
                e,sgmx,sgmy,w(wi),X,Y,OutputTMDNS{wi},Probes,settings,"TM");
            DNSProbeTimeTM(wi)=toc;
        end
    end
   

cPlusDNS  = struct;
cMinusDNS = struct;

if pol == "TE" 
    cPlusDNS.TE  = [cPlusTE{:}];
    cMinusDNS.TE = [cMinusTE{:}];
end

if pol == "TM" 
    cPlusDNS.TM  = [cPlusTM{:}];
    cMinusDNS.TM = [cMinusTM{:}];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POD SOLVE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Load POD modes
    fprintf("Loading In Data\n");
    Modes = PData.Modes;
    if pol=="TE" 
        nrTE=PData.nrTE;
    end
    if pol=="TM" 
        nrTM=PData.nrTM;
    end
    % Load required interpolation libraries
    if pol == "TE"
        TEData = PData.TEData;
    end
    
    if pol == "TM" 
        TMData = PData.TMData;
    end
    
    % TrainingParams are the same information used to build
    % the parameter interpolation grid
    if pol == "TE" 
        TrainingParams = TEData.TrainingParams;
    else
        TrainingParams = TMData.TrainingParams;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get current geometry parameter values from G
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nfeat = numel(TrainingParams);
    p = cell(1,nfeat);
    for i = 1:nfeat
        d = TrainingParams(i).device;
        fieldname = TrainingParams(i).name;
        p{i} = G(d).(fieldname);
    end
    fprintf("Solving for POD\n");
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Solve each requested frequency
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    cPlusTE  = cell(1,numel(w));
    cMinusTE = cell(1,numel(w));
    cPlusTM  = cell(1,numel(w));
    cMinusTM = cell(1,numel(w));
    nInterior = (Nx-2)*(Ny-2);
    if pol == "TE" 
        aOutputTEPOD = complex(zeros(um,numel(w)));
        OutputTEPOD = complex(zeros(nInterior,numel(w)));
        binTE = zeros(1,numel(w));
    end
    if pol == "TM" 
        aOutputTMPOD = complex(zeros(um,numel(w)));
        OutputTMPOD = complex(zeros(nInterior,numel(w)));
        binTM = zeros(1,numel(w));
    end
    
    for wi = 1:numel(w)
        wq = w(wi);
        %=============================================================
        % TE
        %=============================================================
        if pol == "TE" 
            [PODSolveTimeTE(wi),aOutputTEPOD(:,wi),binTE(wi)] = ...
                solvePOD(TEData,wq,p,nrTE,um,"TE");
            tic
            OutputTEPOD(:,wi) = ...
                Modes.etaTE{binTE(wi)}(:,1:um)*aOutputTEPOD(:,wi);
            PODRecoverTimeTE(wi)=toc;
            tic
            [cPlusTE{wi},cMinusTE{wi}] = SolvePODProb( ...
                TEData,binTE(wi),wq,aOutputTEPOD(:,wi),"TE");
            PODProbeTimeTE(wi)=toc;
        end
        %=============================================================
        % Solves for TMTM
        %=============================================================
        if pol == "TM" 
            [PODSolveTimeTM(wi),aOutputTMPOD(:,wi),binTM(wi)] = ...
            solvePOD(TMData,wq,p,nrTM,um,"TM");
             tic
            OutputTMPOD(:,wi) = ...
                Modes.etaTM{binTM(wi)}(:,1:um)*aOutputTMPOD(:,wi);
            PODRecoverTimeTM(wi)=toc;
            tic
            [cPlusTM{wi},cMinusTM{wi}] = SolvePODProb( ...
                TMData,binTM(wi),wq,aOutputTMPOD(:,wi),"TM");
            PODProbeTimeTM(wi)=toc;
        end
    end

    % Expand all reduced coefficient vectors into complete real-space
    % fields. Each bin is one dense matrix-matrix multiplication.
 

    cPlusPOD  = struct;
    cMinusPOD = struct;

    if pol == "TE" 
        cPlusPOD.TE  = [cPlusTE{:}];
        cMinusPOD.TE = [cMinusTE{:}];
    end

    if pol == "TM" 
        cPlusPOD.TM  = [cPlusTM{:}];
        cMinusPOD.TM = [cMinusTM{:}];
    end
 fprintf("Done solving POD\n");
    % Selected-mode full-field LSE at every requested frequency. This uses
    % the same definition as the optional LSE-versus-mode-count sweep.
FieldLSE = struct();
if pol == "TE"
    FieldLSE.TE = zeros(1,numel(w));
    for wi = 1:numel(w)
        difference = OutputTEDNS{wi}-OutputTEPOD(:,wi);
        numerator = sum(difference.*conj(difference),"all");
        denominator = sum( ...
            OutputTEDNS{wi}.*conj(OutputTEDNS{wi}),"all");
        FieldLSE.TE(wi) = real(sqrt(numerator/denominator)*100);
    end
end
if pol == "TM" 
    FieldLSE.TM = zeros(1,numel(w));
    for wi = 1:numel(w)
        difference = OutputTMDNS{wi}-OutputTMPOD(:,wi);
        numerator = sum(difference.*conj(difference),"all");
        denominator = sum( ...
            OutputTMDNS{wi}.*conj(OutputTMDNS{wi}),"all");
        FieldLSE.TM(wi) = real(sqrt(numerator/denominator)*100);
    end
end
%%%%%%%%
   
if options.contour
    if options.contour
        assert(~isnan(options.contourfreq), ...
        "contourfreq must be specified when contour=true.");
    end
      % plot DNS
      f=options.contourfreq;
      wq=f*2*pi;
      [~,wi]=min(abs(w-wq));
      selectedFrequency = w(wi)/(2*pi);

      if pol=="TE"
        DNSC=figure;
        H=zeros(Ny,Nx);
        H(2:end-1,2:end-1)=reshape(OutputTEDNS{wi},Ny-2,Nx-2);
        contourf(X/1e-6,Y/1e-6,abs(H));
        title("DNS TE","FontSize",16);
        xlabel("X ($\mu$m)","FontSize",16,"Interpreter","latex");
        ylabel("Y ($\mu$m)","FontSize",16,"Interpreter","latex");
        box on;
        ax=gca;
        set(ax,'LineWidth',2);
        set(ax,"FontSize",15);
        savefig(DNSC,"./Results/DNS_TE_Contour_at_"+string(selectedFrequency)+".fig");

        PODC=figure;
        H=zeros(Ny,Nx);
        H(2:end-1,2:end-1)=reshape(OutputTEPOD(:,wi),Ny-2,Nx-2);
        contourf(X/1e-6,Y/1e-6,abs(H));
        title("POD TE","FontSize",16);
        xlabel("X ($\mu$m)","FontSize",16,"Interpreter","latex");
        ylabel("Y ($\mu$m)","FontSize",16,"Interpreter","latex");
        box on;
        ax=gca;
        set(ax,'LineWidth',2);
        set(ax,"FontSize",15);
        savefig(PODC,"./Results/POD_TE_Contour_at_"+string(selectedFrequency)+".fig");
      end
      if pol=="TM"
        DNSC=figure;
        E=zeros(Ny,Nx);
        E(2:end-1,2:end-1)=reshape(OutputTMDNS{wi},Ny-2,Nx-2);
        contourf(X/1e-6,Y/1e-6,abs(E));
        title("DNS TM","FontSize",16);
        xlabel("X ($\mu$m)","FontSize",16,"Interpreter","latex");
        ylabel("Y ($\mu$m)","FontSize",16,"Interpreter","latex");
        box on;
        ax=gca;
        set(ax,'LineWidth',2);
        set(ax,"FontSize",15);
        savefig(DNSC,"./Results/DNS_TM_Contour_at_"+string(selectedFrequency)+".fig");

        PODC=figure;
        E=zeros(Ny,Nx);
        E(2:end-1,2:end-1)=reshape(OutputTMPOD(:,wi),Ny-2,Nx-2);
        contourf(X/1e-6,Y/1e-6,abs(E));
        title("POD TM","FontSize",16);
        xlabel("X ($\mu$m)","FontSize",16,"Interpreter","latex");
        ylabel("Y ($\mu$m)","FontSize",16,"Interpreter","latex");
        box on;
        ax=gca;
        set(ax,'LineWidth',2);
        set(ax,"FontSize",15);
        savefig(PODC,"./Results/POD_TM_Contour_at_"+string(selectedFrequency)+".fig");
      end
end

 % COMPUTE LSE per mode
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE LSE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ModeLSE = struct();
if options.LSE
    if pol == "TE" 
        maxM=options.maxM;
        LSETEMax = zeros(1,maxM);
        fprintf('\nStarting TE mode-convergence sweep: %d modes, %d frequencies\n', ...
            maxM,numel(w));
        for m=1:maxM
            modeTimer = tic;
            fprintf('  TE mode %d/%d: solving...\n',m,maxM);
            LSEw= zeros(1,numel(w));
            parfor wi = 1:numel(w)
                wq=w(wi);
                [~,aLSE,bLSE] = solvePOD(TEData,wq,p,nrTE,m,"TE");
                PODsol = Modes.etaTE{bLSE}(:,1:m)*aLSE;
                num=sum((OutputTEDNS{wi} - PODsol).*conj(OutputTEDNS{wi} - PODsol),"all");
                denom=sum(OutputTEDNS{wi}.*conj(OutputTEDNS{wi}),'all');        
                LSEw(wi)=sqrt(num/denom)*100; % convert to a percent
            end
            LSETEMax(m)=max(LSEw);
            fprintf('  TE mode %d/%d: max LSE = %.6g%% (%.1f s)\n', ...
                m,maxM,LSETEMax(m),toc(modeTimer));
        end 
        ModeLSE.TE = LSETEMax;
    end
    
    if pol == "TM" 
        maxM=options.maxM;
        LSETMMax = zeros(1,maxM);
        fprintf('\nStarting TM mode-convergence sweep: %d modes, %d frequencies\n', ...
            maxM,numel(w));
        for m=1:maxM
            modeTimer = tic;
            fprintf('  TM mode %d/%d: solving...\n',m,maxM);
            LSEw= zeros(1,numel(w));
            parfor wi = 1:numel(w)
                wq=w(wi);
                [~,aLSE,bLSE] = solvePOD(TMData,wq,p,nrTM,m,"TM");
                PODsol = Modes.etaTM{bLSE}(:,1:m)*aLSE;
                num=sum((OutputTMDNS{wi} - PODsol).*conj(OutputTMDNS{wi} - PODsol),"all");
                denom=sum(OutputTMDNS{wi}.*conj(OutputTMDNS{wi}),'all');        
                LSEw(wi)=sqrt(num/denom)*100; % convert to a percent
            end
            LSETMMax(m)=max(LSEw);
            fprintf('  TM mode %d/%d: max LSE = %.6g%% (%.1f s)\n', ...
                m,maxM,LSETMMax(m),toc(modeTimer));
        end 
        ModeLSE.TM = LSETMMax;
    end
    
    LSEFig=figure;
    if pol=="TE"
        semilogy(1:maxM,LSETEMax,"LineWidth",2);
    elseif pol=="TM"
        semilogy(1:maxM, LSETMMax,"LineWidth",2);
    
    end
    xlabel("Number of Modes","FontSize",16);
    ylabel("Maximum LSE over frequency (%)","FontSize",16);
    legend("location","best","FontSize",16);
    box on;
    ax=gca;
    set(ax,'YScale','log');
    set(ax,"LineWidth",2);
    set(ax,"FontSize",14);
    savefig(LSEFig,"./Results/LSE.fig")

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create and save timing table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

omega = w(:);
frequency = omega/(2*pi);
wavelength = 2*pi*c0./omega;

if pol == "TE"

    PODTotalTimeTE = PODSolveTimeTE + PODRecoverTimeTE + PODProbeTimeTE;
    DNSTotalTimeTE = MTimeTE + DNSTESolveTime + DNSProbeTimeTE;
    speedupTE = DNSTotalTimeTE./PODTotalTimeTE;

    Polarization = repmat("TE",numel(w),1);
    NumModes = repmat(um,numel(w),1);
    TimingTable = table( ...
    Polarization, ...
    NumModes, ...
    omega, ...
    frequency, ...
    wavelength, ...
    MTimeTE(:), ...
    DNSTESolveTime(:), ...
    DNSProbeTimeTE(:), ...
    DNSTotalTimeTE(:), ...
    PODSolveTimeTE(:), ...
    PODRecoverTimeTE(:), ...
    PODProbeTimeTE(:), ...
    PODTotalTimeTE(:), ...
    speedupTE(:), ...
    'VariableNames',{ ...
    'Polarization', ...
    'NumModes', ...
    'AngularFrequency_rad_s', ...
    'Frequency_Hz', ...
    'Wavelength_m', ...
    'DNS_AssemblyTime_s', ...
    'DNS_SolveTime_s', ...
    'DNS_ProbeTime_s', ...
    'DNS_TotalTime_s', ...
    'POD_SolveTime_s', ...
    'POD_ReconstructionTime_s', ...
    'POD_ProbeTime_s', ...
    'POD_TotalTime_s', ...
    'Speedup'});

elseif pol == "TM"

    PODTotalTimeTM = PODSolveTimeTM + PODRecoverTimeTM + PODProbeTimeTM;
    DNSTotalTimeTM = MTimeTM + DNSTMSolveTime + DNSProbeTimeTM;
    speedupTM = DNSTotalTimeTM./PODTotalTimeTM;

    Polarization = repmat("TM",numel(w),1);
    NumModes = repmat(um,numel(w),1);
    TimingTable = table( ...
    Polarization, ...
    NumModes, ...
    omega, ...
    frequency, ...
    wavelength, ...
    MTimeTM(:), ...
    DNSTMSolveTime(:), ...
    DNSProbeTimeTM(:), ...
    DNSTotalTimeTM(:), ...
    PODSolveTimeTM(:), ...
    PODRecoverTimeTM(:), ...
    PODProbeTimeTM(:), ...
    PODTotalTimeTM(:), ...
    speedupTM(:), ...
    'VariableNames',{ ...
    'Polarization', ...
    'NumModes', ...
    'AngularFrequency_rad_s', ...
    'Frequency_Hz', ...
    'Wavelength_m', ...
    'DNS_AssemblyTime_s', ...
    'DNS_SolveTime_s', ...
    'DNS_ProbeTime_s', ...
    'DNS_TotalTime_s', ...
    'POD_SolveTime_s', ...
    'POD_ReconstructionTime_s', ...
    'POD_ProbeTime_s', ...
    'POD_TotalTime_s', ...
    'Speedup'});
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save the timing table
if ~isfolder("./Results")
    mkdir("./Results");
end

writetable(TimingTable,"./Results/TimingResults.xlsx");
writetable(TimingTable,"./Results/TimingResults.csv");
save("./Results/TimingResults.mat","TimingTable");


end



    
