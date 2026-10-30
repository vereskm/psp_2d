function plotGeom(g,Lx,Ly,unit)
lw=2;
if nargin < 4
    unit = "m";
end

switch unit

    case "mm"
        S=1E-3;
    case "um"
        S=1E-6;
    case "nm"
        S=1E-9;
    case "m"
        S=1;
end

for n = 1:numel(g)
      % pick a color for each component
      if isfield(g(n),"color") && ~isempty(g(n).color)
         color = g(n).color;
      else
         color = [0.7 0.7 0.7];   % default gray
      end

    switch g(n).type
      

        case "rod"
            rectangle( ...
                'Position',[ ...
                g(n).x/S-g(n).r/S, ...
                g(n).y/S-g(n).r/S, ...
                2*g(n).r/S, ...
                2*g(n).r/S], ...
                'Curvature',[1 1],'LineWidth',lw,"FaceColor",color);
        case "waveg"
            rectangle( ...
                'Position',[ ...
                g(n).x/S-g(n).w/2/S, ...
                g(n).y/S-g(n).h/2/S, ...
                g(n).w/S, ...
                g(n).h/S],'LineWidth',lw,"FaceColor",color);

        case "ring"
            rin  = g(n).r - g(n).w/2;
            rout = g(n).r + g(n).w/2;

            rectangle( ...
                'Position',[ ...
                g(n).x/S-rout/S, ...
                g(n).y/S-rout/S, ...
                2*rout/S, ...
                2*rout/S], ...
                'Curvature',[1 1],'LineWidth',lw,"FaceColor",color);

            rectangle( ...
                'Position',[ ...
                g(n).x/S-rin/S, ...
                g(n).y/S-rin/S, ...
                2*rin/S, ...
                2*rin/S], ...
                'Curvature',[1 1],'LineWidth',lw,"FaceColor","w");

        otherwise
            error("Unknown geometry type: %s",g(n).type)

    end

end
fontsize=16;
axfontsize=14;
ax=gca;
box on;
set(ax,'Fontsize',axfontsize);
set(ax,'LineWidth',2);
axis equal;
switch  unit
    case "mm"
        axis([0,Lx/1E-3,0,Ly/1E-3]);
        xlabel("X (mm)","FontSize",fontsize);
        ylabel("Y (mm)","FontSize",fontsize);
    case "um"
       axis([0,Lx/1E-6,0,Ly/1E-6]);
       xlabel("X ($\mu$m)","FontSize",fontsize,'Interpreter','latex');
       ylabel("Y ($\mu$m)","FontSize",fontsize,'Interpreter','latex');
    case "nm"
        axis([0,Lx/1E-9,0,Ly/1E-9]);
        xlabel("X (nm)","FontSize",fontsize);
        ylabel("Y (nm)","FontSize",fontsize);
    otherwise


        axis([0,Lx,0,Ly]);
        xlabel("X","FontSize",fontsize);
        ylabel("Y","FontSize",fontsize);
end

end