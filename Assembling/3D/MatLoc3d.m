function [A] = MatLoc3d(T1,T2,k)

% Returns the matrix Aij, 1 <= i,j <= 2, where
% Aij = int_{[A,B]x[C,D]} -1/(2*pi)log(||x-y||) (phi_i)'(x) (psi_j)'(y)dxdy
% Where phi_i (resp. psi_i) is the shape function on [A,B] (resp. [C,D])
% attached to the vertex nb. i, and A = [xA,yA], ..., D = [xD,yD].


% Cheating using Gypsilab
m1 = msh(T1,[1,2,3]);
m2 = msh(T2,[1,2,3]);

Vh1 = fem(m1,'P1');
Vh2 = fem(m2,'P1');

Gamma1 = dom(m1,7);
Gamma2 = dom(m2,7);

if (k == 0)
    Gxy = @(X,Y)(femGreenKernel(X,Y,'[1/r]',0));
    c = 1/(4*pi);
else
    Gxy = @(X,Y)(femGreenKernel(X,Y,'[exp(ikr/r)]',k));
    c = 1/(4*pi);
end

A = c*integral(Gamma1,Gamma2,nxgrad(Vh1),Gxy,nxgrad(Vh2)) ...
    - k^2*c*integral(Gamma1,Gamma2,ntimes(Vh1),Gxy,ntimes(Vh2));
A = A + c*regularize(Gamma1,Gamma2,nxgrad(Vh1),'[1/r]',nxgrad(Vh2))...
    - k^2*(-1)/(2*pi)*regularize(Gamma1,Gamma2,ntimes(Vh1),'[1/r]',ntimes(Vh2));


end

