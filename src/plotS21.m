function fig = plotS21(R,inputProbe,outputProbe)
    %PLOTS21 Compare DNS and ROM forward transmission for every random sample.
    % Port modes are power normalized. Defaults: input probe 1, output probe 2.
    if nargin<2, inputProbe=1; end
    if nargin<3, outputProbe=2; end
    validateattributes(inputProbe,{'numeric'},{'scalar','integer','positive'});
    validateattributes(outputProbe,{'numeric'},{'scalar','integer','positive'});
    w=R.frequencies(:);
    validateattributes(w,{'numeric'},{'real','finite','positive','nonempty'});
    assert(numel(R.cpDNS)==numel(R.cpROM),'DNS and ROM sample counts must match.');
    C=physics.constants();
    [wavelength,order]=sort(2*pi*C.c0./w*1e9);

    fig=figure;
    hold on;
    for m=1:numel(R.cpDNS)
        dns=R.cpDNS{m};
        rom=R.cpROM{m};
        assert(size(dns,2)==numel(w) && isequal(size(dns),size(rom)), ...
            'Port matrices must have probes as rows and frequencies as columns.');
        assert(max(inputProbe,outputProbe)<=size(dns,1), ...
            'Requested probe index is not present in the results.');
        assert(all(abs(dns(inputProbe,:))>0) && all(abs(rom(inputProbe,:))>0), ...
            'S21 is undefined when the input forward amplitude is zero.');
        s21DNS=dns(outputProbe,:)./dns(inputProbe,:);
        s21ROM=rom(outputProbe,:)./rom(inputProbe,:);
        h=plot(wavelength,20*log10(abs(s21DNS(order))),'-', ...
            'LineWidth',2,'DisplayName',sprintf('DNS sample %d',m));
        plot(wavelength,20*log10(abs(s21ROM(order))),'--', ...
            'Color',h.Color,'LineWidth',2,'DisplayName',sprintf('ROM sample %d',m));
    end
    xlabel('Wavelength (nm)');
    ylabel('|S_{21}| (dB)');
    grid on;
    box on;
    legend('show','Location','best');
end
