function MonteCarloR = runMonteCarlo(SolveSettings,I,fields,CenterValues,P,Nrand,Seed)
    %RUNMONTECARLO Compare DNS and ROM for randomly perturbed geometries.
    % I are the geometry indices;
    % fields are there names
    % centervalues are the center values to try
    %P are the percentages in each direction
    % Nrand is the total number of random tests
    problemSample=SolveSettings.problem;
    POD=SolveSettings.POD;
    frequencies=SolveSettings.w;
    nModes=SolveSettings.nModes;
    MonteCarloR.Nrand=Nrand; % set the Nrand
    MonteCarloR.fieldLSE = zeros(nModes,numel(frequencies),Nrand);

    MonteCarloR.uROM = cell(nModes,1);
  

    % define the random seed
    rng(Seed);
    MonteCarloR.cpDNS=cell(Nrand,1);
    MonteCarloR.cnDNS=cell(Nrand,1);
    MonteCarloR.cpROM=cell(Nrand,1);
    MonteCarloR.cnROM=cell(Nrand,1);
    % Full interior fields: first index is sample, second is frequency.
    MonteCarloR.uROM = cell(Nrand,1);
    MonteCarloR.uDNS=cell(Nrand,1);
    MonteCarloR.geometry=cell(Nrand,1);
    MonteCarloR.gridData=mesh.create(SolveSettings.problem);
    MonteCarloR.frequencies=frequencies(:);
    % initialize the total times
    romTotalTime = zeros(numel(frequencies),1);
    romPortTime = zeros(numel(frequencies),1);
    romFieldTime = zeros(numel(frequencies),1);
    MonteCarloR.Nrand=Nrand;
    MonteCarloR.devices=I;
    MonteCarloR.fields=fields;
    MonteCarloR.values=zeros(Nrand,numel(I));
    for m = 1:Nrand
        % define a sample problem
        problemSample = SolveSettings.problem;
        % for each of the devices
        for k = 1:numel(I)
            % get the
            value = CenterValues(k) *(1+P(k)*(2*rand-1));
            MonteCarloR.values(m,k)=value;
            % put in this new value
            problemSample.geometry(I(k)).(fields(k)) = value;
        end
    
         % get the total number of parameters
         parameters = zeros(1,numel(POD.trainingParameters));
         % loop through the parameters
        for k = 1:numel(POD.trainingParameters)
            feature = POD.trainingParameters(k); % struct defining the kth trained parameter

            parameters(k) = ...
                problemSample.geometry(feature.device).(feature.field);
        end
        %solve DNS
        [uDNS,dnsInfo,cpDNS,cnDNS] = ...
            online.runDNS(problemSample,frequencies);
        MonteCarloR.cpDNS{m}=cpDNS;
        MonteCarloR.cnDNS{m}=cnDNS;
        MonteCarloR.uDNS{m}=uDNS;
        MonteCarloR.geometry{m}=problemSample.geometry;
        % put in the DNS times
        MonteCarloR.time.dnsSetup(m) = dnsInfo.setupTime;
        MonteCarloR.time.dnsAssembly(:,m) = dnsInfo.assemblyTime;
        MonteCarloR.time.dnsSolve(:,m) = dnsInfo.solveTime;
        MonteCarloR.time.dnsPorts(:,m) = dnsInfo.portTime;
        MonteCarloR.time.dnsTotal(:,m)=dnsInfo.assemblyTime+dnsInfo.solveTime+dnsInfo.portTime;
        MonteCarloR.time.dnsSweepTime(m) = dnsInfo.totalTime;
        for j = 1:numel(frequencies)
            for nm=1:nModes
                [a,romInfo,uROM] = online.solveROM( ...
                    POD,frequencies(j),parameters,nm);
                MonteCarloR.fieldLSE(nm,j,m) = ...
                    100*norm(uDNS{j}-uROM)/norm(uDNS{j});

                  
            end
            
         
             MonteCarloR.uROM{m}{j}=uROM;
            % Port coefficients:
            cpROM(:,j) = romInfo.cPlus;
            cnROM(:,j) = romInfo.cMinus;
            
            
            MonteCarloR.time.romAssembly(j,m) = romInfo.assemblyTime;
            MonteCarloR.time.romSolve(j,m) = romInfo.solveTime;
            MonteCarloR.time.romPorts(j,m) = romInfo.portTime;
            MonteCarloR.time.romReconstruction(j,m) = ...
            romInfo.reconstructionTime;
            MonteCarloR.time.romTotal(j,m) = romInfo.totalTime;
        end
        MonteCarloR.time.romSweepTime(m)=sum(MonteCarloR.time.romTotal(:,m));
        MonteCarloR.cpROM{m}=cpROM;
        MonteCarloR.cnROM{m}=cnROM;
    end


   
end





   

