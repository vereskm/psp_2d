function [MonteCarloR] = RunMonteCarlo(SolveSettings,I,fields,CenterValue,P,Nrand,Seed)
%RUNMONTECARLO Gets port amplitudes after changing some device by a certain
%range
% Inputs are G (the device)
% I: the element being changed can be a vector if more than one
% fields: name of the field being changed
% CenterValues: values you are changing by a percentage
% Change ( Percentage change in both directions symetric)
% P: (percentage change in both directions)
assert(P>0,"Percentage must be possitive");
npc=numel(I); % get number of parameters which are changed
rng(seed); % set up the random number generator
randomPar=zeros(length(I),Nrand); % get the random parameters
for i=1:length(I)
    S=-P(i);
    E=P(i);
    RandomP=S+(E-S)*rand(1,Nrand); % get the random probabilities
    randomPar(i,:)=(1+RandomP).*CenterValues(i);
end
% Initialize the output cell arrays
cPlusDNS = cell(1, Nrand);
cMinusDNS = cell(1, Nrand);
cPlusPOD = cell(1, Nrand);
cMinusPOD = cell(1, Nrand);
FieldLSE = cell(1, Nrand);
ModeLSE = cell(1, Nrand);
for m=1:Nrand
    caseTimer=tic; % run the timer
    for p=1:npc
        Value=randomPar(p,m); % get the parameter value
        SolveSettings.G(p).(fields(p))=Value;
    end
    fprintf('\n=== Running Monte Carlo %d out of %d ===\n',m,Nrand);
    [cPlusDNS{m},cMinusDNS{m},cPlusPOD{m},cMinusPOD{m},FieldLSE{m},ModeLSE{m}]= ...
        Solve(SolveSettings.settings,SolveSettings.PData,-+, ...
        SolveSettings.w,SolveSettings.Sources,SolveSettings.Probes, ...
        SolveSettings.um,SolveSettings.maxM,SolveSettings.nworkers);

     fprintf('=== Finished Monte Carlo run %d/%d in %.1f minutes ===\n', ...
        m,Nrand,toc(caseTimer)/60);
end
% get a summary
MonteCarloR.cPlusDNS=cPlusDNS;
MonteCarloR.cMinusDNS=cMinusDNS;
MonteCarloR.cPlusPOD=cPlusPOD;
MonteCarloR.cMinusPOD=cMinusPOD;
MonteCarloR.FieldLSE=FieldLSE;
MonteCarloR.ModeLSE=ModeLSE;
MonteCarloR.fields=fields;
MonteCarloR.randomPar=randomPar;

end