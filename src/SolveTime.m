function [cPlusDNS,cMinusDNS,cPlusPOD,cMinusPOD] = SolveTime(settings,PData,G,w,Sources,Probes,um)
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
    end



    if ~isfolder("./Results")
        mkdir("./Results");
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
    for wi=1:numel(w)
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
    for wi = 1:numel(w)
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



    
