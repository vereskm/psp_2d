function plotSources(Sources,Lx,Ly,unit)

lw = 3;

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

for n = 1:numel(Sources)

    p1 = Sources(n).p1 / S;
    p2 = Sources(n).p2 / S;

    % Plot source plane
    plot([p1(1),p2(1)], ...
         [p1(2),p2(2)], ...
         'k-', ...
         'LineWidth',lw);

    % Center of source
    pc = (p1+p2)/2;

    % Arrow length
    arrowLength = 0.05*max(Lx,Ly)/S;

    % Propagation direction
    switch Sources(n).dir

        case "right"
            dx = arrowLength;
            dy = 0;

        case "left"
            dx = -arrowLength;
            dy = 0;

        case "up"
            dx = 0;
            dy = arrowLength;

        case "down"
            dx = 0;
            dy = -arrowLength;

        otherwise
            error("Unknown source direction: %s",Sources(n).dir);
    end

    % Plot propagation arrow
    quiver(pc(1),pc(2),dx,dy,0, ...
           'k', ...
           'LineWidth',2, ...
           'MaxHeadSize',0.8);

end

% Match plotGeom axes
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