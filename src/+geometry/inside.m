%TODO: Expand to more geometries in the future if need be
function inside=inside(g,X,Y)
    % determine if you are inside a geometry
    switch g.type
        case "rod"
            inside=(X-g.x).^2+(Y-g.y).^2<=g.r.^2;
        case "waveg"
            inside=abs(X-g.x)<g.w/2&abs(Y-g.y)<=g.h/2;
        case "ring"
             rin=g.r-g.w/2;
             rout=g.r+g.w/2;
             inside=((X-g.x).^2+(Y-g.y).^2>=rin.^2)&(X-g.x).^2+(Y-g.y).^2<=rout.^2;
        otherwise
            error("Unknow geometry type")
    end

end
