function [cPlus,cMinus] = ModeOverlap(q,E,H,Em,Hm,pHat)
% This part of the program was writen with ChatGPT (Verified by Martin
% though)
% q     : N-by-1 probe coordinate
% E,H   : N-by-3 computed fields
% Em,Hm : N-by-3 positive-direction mode fields
% pHat  : 1-by-3 positive port direction

pHat = pHat(:);  % 3-by-1

% Cross products at every probe point
C1 = cross(E,conj(Hm),2); % compute the first term ExH_m^*
C2 = cross(conj(Em),H,2); % compute the second term E_m^*times H
Cm = cross(Em,conj(Hm),2); % compute the denominator

% Dot each cross product with the port direction (this grabs the
% corresponding component
f1 = C1*pHat; 
f2 = C2*pHat;
fm = Cm*pHat;

% Port normalization
Nmode = trapz(q,fm);

% Forward and backward modal amplitudes
cPlus  = trapz(q,f1 + f2)/(2*Nmode);
cMinus = trapz(q,f1 - f2)/(2*Nmode);
end