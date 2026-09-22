function []=saveSnapshots(problem,training,snapshots)
    %SAVESNAPSHOTS save training snapshots to a file
    pol=problem.polarization;
    outputDirectory =training.outputDir;
    if training.deletePrevious && isfolder(outputDirectory)
        mkdir(outputDirectory);
    end

    if ~isfolder(outputDirectory)
        mkdir(outputDirectory);
    end
   
    dateString = string(datetime( ...
        "now","Format","yyyy-MM-dd_HH-mm-ss"));
    baseName = string(training.randomSeed) + "_Training_" + dateString;
    if pol == "TE"
        save( ...
            fullfile(outputDirectory,"TE_" + baseName + ".mat"), ...
            "snapshots","training","problem","-v7.3");
    end
    if pol == "TM"
        save( ...
            fullfile(outputDirectory,"TM_" + baseName + ".mat"), ...
            "snapshots","training","problem","-v7.3");
    end

end
