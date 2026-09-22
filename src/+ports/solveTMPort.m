function [phi,beta] = solveTMPort(e,x,w,sgm)
% This function compute the TM source (It compute the largest eigenvalue)
%INPUTS: 
% e is the permitivity
% x is the grid
% w is the angular frequency
% sgm is the optional transverse PML conductivity on the same 1D grid
C=physics.constants();
e0=C.e0;
c0=C.c0;
if nargin < 4 || isempty(sgm)
    sgm = zeros(size(e));
end
e = reshape(e,[],1);
x = reshape(x,[],1);
sgm = reshape(sgm,[],1);
n=numel(x); % get the number of points
dx=x(2)-x(1);
row=zeros(3*(n-2)-2,1);
col=zeros(3*(n-2)-2,1);
data=complex(zeros(3*(n-2)-2,1));
k_0=w/c0;
sinv = 1./(1+1i.*sgm./(w.*e0));
I=1;
for i=2:n-1
    p=i-1;
    ec=e(i);
    aE=sinv(i).*(sinv(i+1)+sinv(i))/2;
    aW=sinv(i).*(sinv(i-1)+sinv(i))/2;
    row(I)=p;
    col(I)=p;
    data(I)=-(aE+aW)/dx^2+k_0.^2*ec;
    I=I+1;
    if i==2
        row(I)=p;
        col(I)=p+1;
        data(I)=aE/dx^2;
        I=I+1;
    elseif i==n-1
        row(I)=p;
        col(I)=p-1;
        data(I)=aW/dx^2;
        I=I+1;
    else
        row(I)=p;
        col(I)=p+1;
        data(I)=aE/dx^2;
        I=I+1;

        row(I)=p;
        col(I)=p-1;
        data(I)=aW/dx^2;
        I=I+1;
    end 
end
    A=sparse(row,col,data,(n-2),(n-2));
    [phinB,beta2]=eigs(A,1,"largestreal");
    phi=[0;phinB;0]; % put back the boundary
    beta=sqrt(beta2);
    if real(beta) < 0
        beta = -beta;
    end
    % normalize the eigenvector
    [~,iMax] = max(abs(phi));
    if max(abs(sgm)) < eps
        phi = real(phi);
    else
        phi = phi.*exp(-1i*angle(phi(iMax)));
    end
    if real(phi(iMax)) < 0
        phi = -phi;
    end
    phi = phi/max(abs(phi));

end

