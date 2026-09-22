function [PData] = LoadLib(pol)
     
    fprintf("Loading In Data\n");
    PData.Modes = load('./Lib/PODModes.mat');
    if pol=="TE" 
        PData.nrTE=size(PData.Modes.etaTE{1},2);
    end
    if pol=="TM" 
        PData.nrTM=size(PData.Modes.etaTM{1},2);
    end
    % Load required interpolation libraries
    if pol == "TE"
        PData.TEData = load('./Lib/InterpTE.mat');
    end
    
    if pol == "TM" 
        PData.TMData = load('./Lib/InterpTM.mat');
    end
    
    % TrainingParams are the same information used to build
    % the parameter interpolation grid
    if pol == "TE" 
        PData.TrainingParams = PData.TEData.TrainingParams;
    else
        PData.TrainingParams = PData.TMData.TrainingParams;
    end
    
   
end
