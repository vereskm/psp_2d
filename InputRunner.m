% Run after defining or loading problem, POD, and frequencies.
% Keep these input variables in the workspace.
clc; close all;
projectRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(projectRoot,'src'));
device;
plotContours=true;
resultsDir = fullfile(projectRoot,'Results');
if ~isfolder(resultsDir)
    mkdir(resultsDir);
end
% define frequencies
C = physics.constants();
wmin = 2*pi*C.c0/lambda1;
wmax = 2*pi*C.c0/lambda2;
nfreq=120;
frequencies=linspace(wmin,wmax,nfreq);
nworkers=30;


data = load(fullfile(projectRoot,'Problem.mat'),'problem');
problem = data.problem;

dataPOD = load(fullfile(projectRoot,'Model','POD.mat'),'POD');
POD = dataPOD.POD;

S.problem = problem;
S.POD = POD;
S.w = frequencies;
S.nModes = 20;
Nrand=10;
Seed=1;
R = online.runMonteCarlo( ...
    S,2,"er",problem.geometry(2).er,0.007,Nrand,Seed);
save(fullfile(resultsDir,'MonteCarlo.mat'),'R','-v7.3');
createTimingTables(R,resultsDir);

%% Plot the scattering output S_21.
C = physics.constants();
[wavelength_nm,order] = sort(2*pi*C.c0 ./ R.frequencies(:) * 1e9);
S21 = figure;
hold on;
for m = 1:numel(R.cpDNS)
    baseColor = rand(1,3);
    dnsColor = 0.6*baseColor + 0.4; % Mix with white.
    podColor = 0.65*baseColor;      % Darken toward black.

    DNS_S21 = R.cpDNS{m}(2,:) ./ R.cpDNS{m}(1,:);
    POD_S21 = R.cpROM{m}(2,:) ./ R.cpROM{m}(1,:);
    plot(wavelength_nm,20*log10(abs(DNS_S21(order))), ...
        'LineWidth',2,'Color',dnsColor, ...
        'DisplayName',sprintf('DNS sample %d',m));
    plot(wavelength_nm,20*log10(abs(POD_S21(order))),'--', ...
        'LineWidth',2,'Color',podColor, ...
        'DisplayName',sprintf('POD sample %d',m));
end
xlabel('Wavelength (nm)','FontSize',18);
ylabel('|S_{21}| (dB)','FontSize',18);
ax = gca;
box on;
set(ax,'Linewidth',2);
grid on;
legend('show','Location','best',"FontSize",18);
if plotContours==true;
    gridData=R.gridData;
    for m = 1:numel(R.cpDNS)
        DNS_S21 = R.cpDNS{m}(2,:) ./ R.cpDNS{m}(1,:);
        POD_S21 = R.cpROM{m}(2,:) ./ R.cpROM{m}(1,:);
        [~,minIndex]=min(abs(DNS_S21).^2);
        DNSC=figure;
        contourDataDNS = R.uDNS{m}{minIndex};
        contourDataDNS=padarray(reshape(contourDataDNS,gridData.Ny-2,gridData.Nx-2),1);
        contourf(gridData.X/1E-6,gridData.Y/1E-6,abs(contourDataDNS));
        box on;
        ax=gca;
        set(ax,"FontSize",14);
        xlabel("X $(\mu m)$","Interpreter","latex","FontSize",18);
        ylabel("Y $(\mu m)$","Interpreter","latex","FontSize",18);
        axis equal tight;
        PODC=figure;
        contourDataPOD=R.uROM{m}{minIndex};
        contourDataPOD=padarray(reshape(contourDataPOD,gridData.Ny-2,gridData.Nx-2),1);
        contourf(gridData.X/1E-6,gridData.Y/1E-6,abs(contourDataPOD));
        box on;
        ax=gca;
        set(ax,"FontSize",14);
        xlabel("X $(\mu m)$","Interpreter","latex","FontSize",18);
        ylabel("Y $(\mu m)$","Interpreter","latex","FontSize",18);
        axis equal tight;
        savefig(DNSC,fullfile(resultsDir, ...
            sprintf('DNS_Contour_Sample_%d.fig',m)));

        savefig(PODC,fullfile(resultsDir, ...
        sprintf('POD_Contour_Sample_%d.fig',m)));
    end
end

%% Get the Q factors and their errors for POD and DNS.
% Half-prominence widths estimate loaded Q for resolved Lorentzian dips.
% Work in linear power, not dB. Wavelength ratios approximate frequency Q.
nSamples = numel(R.cpDNS);
lambdaDNS = cell(nSamples,1);
lambdaPOD = cell(nSamples,1);
widthDNS = cell(nSamples,1);
widthPOD = cell(nSamples,1);
QDNS = cell(nSamples,1);
QPOD = cell(nSamples,1);
Qerror = cell(nSamples,1);
matchedPOD = cell(nSamples,1);
assert(numel(wavelength_nm)>=3 && all(diff(wavelength_nm)>0), ...
    'Q extraction needs at least three distinct wavelength points.');
% Adjust this maximum allowed DNS/POD resonance-center difference if needed.
matchTolerance_nm = 2*median(diff(wavelength_nm));

for m = 1:nSamples
    DNS_T = abs(R.cpDNS{m}(2,:) ./ R.cpDNS{m}(1,:)).^2;
    POD_T = abs(R.cpROM{m}(2,:) ./ R.cpROM{m}(1,:)).^2;
    DNS_T = DNS_T(:);
    POD_T = POD_T(:);

    [~,lambdaDNS{m},widthDNS{m}] = findpeaks( ...
        -DNS_T(order),wavelength_nm,'WidthReference','halfprom');
    [~,lambdaPOD{m},widthPOD{m}] = findpeaks( ...
        -POD_T(order),wavelength_nm,'WidthReference','halfprom');
    QDNS{m} = lambdaDNS{m} ./ widthDNS{m};
    QPOD{m} = lambdaPOD{m} ./ widthPOD{m};

    % Match mutual nearest resonances. Leave missing/ambiguous matches as NaN.
    Qerror{m} = nan(size(QDNS{m}));
    matchedPOD{m} = nan(size(QDNS{m}));
    tolerance = matchTolerance_nm;
    if numel(lambdaDNS{m})>1
        tolerance = min(tolerance,min(diff(lambdaDNS{m}))/2);
    end
    for q = 1:numel(QDNS{m})
        if isempty(lambdaPOD{m}), continue; end
        [distance,p] = min(abs(lambdaPOD{m}-lambdaDNS{m}(q)));
        [~,nearestDNS] = min(abs(lambdaDNS{m}-lambdaPOD{m}(p)));
        if nearestDNS==q && distance<=tolerance
            matchedPOD{m}(q) = p;
            Qerror{m}(q) = 100*abs(QPOD{m}(p)-QDNS{m}(q))/QDNS{m}(q);
        end
    end
end


%% Create one Q-factor table for each sample.
QTables = cell(nSamples,1);
qTableFile = fullfile(resultsDir,'QFactors.xlsx');
for m = 1:nSamples
    Resonance = (1:numel(QDNS{m})).';
    DNS_Wavelength_nm = lambdaDNS{m}(:);
    DNS_Width_nm = widthDNS{m}(:);
    DNS_Q = QDNS{m}(:);
    POD_Wavelength_nm = nan(size(DNS_Q));
    POD_Width_nm = nan(size(DNS_Q));
    POD_Q = nan(size(DNS_Q));
    Q_Error_percent = Qerror{m}(:);

    for q = 1:numel(DNS_Q)
        p = matchedPOD{m}(q);
        if ~isnan(p)
            POD_Wavelength_nm(q) = lambdaPOD{m}(p);
            POD_Width_nm(q) = widthPOD{m}(p);
            POD_Q(q) = QPOD{m}(p);
        end
    end

    QTables{m} = table(Resonance,DNS_Wavelength_nm,POD_Wavelength_nm, ...
        DNS_Width_nm,POD_Width_nm,DNS_Q,POD_Q,Q_Error_percent);
    writetable(QTables{m},qTableFile, ...
        'Sheet',sprintf('Sample_%d',m),'WriteMode','overwritesheet');
end

%% Compute the relative field L2 error.
nSamples = numel(R.uDNS);
nFrequencies = numel(R.frequencies);
nSamples=R.Nrand;
fieldLSE = R.fieldLSE;

% Plot in the same wavelength order as the S_21 curves.

figure;
hold on;
for m=1:nSamples
    semilogy(1:S.nModes,max(fieldLSE(:,:,m),[],2),'LineWidth',2);
end

xlabel('Number of POD Modes',"FontSize",18);
ylabel('LSE (%)',"FontSize",18);
box on;
ax=gca;
set(ax,'Yscale','log');
set(ax,'FontSize',14);
grid on;
