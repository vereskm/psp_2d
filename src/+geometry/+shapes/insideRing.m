function mask = insideRing(g,X,Y)
    rin=g.r-g.w/2; 
    rout=g.r+g.w/2;
    mask=((X-g.x).^2+(Y-g.y).^2>=rin.^2)&(X-g.x).^2+(Y-g.y).^2<=rout.^2;
end