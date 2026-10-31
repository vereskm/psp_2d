function mask = insideWaveguide(g,X,Y)
    mask=abs(X-g.x)<g.w/2&abs(Y-g.y)<=g.h/2;
end