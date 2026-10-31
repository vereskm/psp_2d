function mask = insideRod(g,X,Y)
    mask=(X-g.x).^2+(Y-g.y).^2<=g.r.^2;
end