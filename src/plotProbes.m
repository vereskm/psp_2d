function plotProbes(Probes,Lx,Ly,unit)

if nargin < 4
    unit = "m";
end

% Unit scaling
switch unit
    case "mm"
        S = 1e-3;
    case "um"
        S = 1e-6;
    case "nm"
        S = 1e-9;
    case "m"
        S = 1;
    otherwise
        error("Unknown unit: %s",unit);
end

hold on;

for pi = 1:numel(Probes)

    p1 = Probes(pi).p1 / S;
    p2 = Probes(pi).p2 / S;

    % Plot probe plane
    plot( ...
        [p1(1),p2(1)], ...
        [p1(2),p2(2)], ...
        "-.", ...
        "Color",[0.10,0.35,0.90], ...
        "LineWidth",2, ...
        "DisplayName","Probe "+pi);

    % Probe center
    pc = 0.5*(p1+p2);

    % Label probe
    text( ...
        pc(1), ...
        pc(2), ...
        "  P_"+pi, ...
        "Color",[0.05,0.20,0.75], ...
        "FontSize",12, ...
        "FontWeight","bold", ...
        "VerticalAlignment","top");

end

% Match domain axes
switch unit
    case "mm"
        axis([0,Lx/1e-3,0,Ly/1e-3]);

    case "um"
        axis([0,Lx/1e-6,0,Ly/1e-6]);

    case "nm"
        axis([0,Lx/1e-9,0,Ly/1e-9]);

    case "m"
        axis([0,Lx,0,Ly]);
end

end