function [] = createTimingTables(R,outputDir)
    %CREATETIMINGTABLES Tabulate recorded Monte Carlo DNS/ROM component times.
    for m=1:numel(R.cpDNS)
        frequencies=R.frequencies(:);
        DNS_Setup=R.time.dnsSetup(m)*ones(size(frequencies));
        DNS_Assembly=R.time.dnsAssembly(:,m);
        DNS_Solve=R.time.dnsSolve(:,m);
        DNS_Ports=R.time.dnsPorts(:,m);
        DNS_Total=R.time.dnsTotal(:,m);


        ROM_Assembly=R.time.romAssembly(:,m);
        ROM_Solve=R.time.romSolve(:,m);
        ROM_Ports=R.time.romPorts(:,m);
        ROM_Reconstruction=R.time.romReconstruction(:,m);
        ROM_Total=R.time.romTotal(:,m);


        T = table(frequencies,DNS_Setup, ...
        DNS_Assembly,DNS_Solve,DNS_Ports,DNS_Total, ...
        ROM_Assembly,ROM_Solve,ROM_Ports, ...
        ROM_Reconstruction,ROM_Total);

        sheetName = sprintf('Sample_%d',m);
        if ~isfolder(outputDir)
            mkdir(outputDir);
        end

        filename = fullfile(outputDir,'TimingResults.xlsx');
        writetable(T,filename,'Sheet',sheetName);
    end


end