function [A] = MatLoc2d(xA,yA,xB,yB,xC,yC,xD,yD,k)

% Returns the matrix Aij, 1 <= i,j <= 2, where
% Aij = int_{[A,B]x[C,D]} -1/(2*pi)log(||x-y||) (phi_i)'(x) (psi_j)'(y)dxdy
% Where phi_i (resp. psi_i) is the shape function on [A,B] (resp. [C,D])
% attached to the vertex nb. i, and A = [xA,yA], ..., D = [xD,yD].


% Cheating using Gypsilab
V = [xA, yA 0; xB, yB 0; xC yC 0; xD, yD 0];
E = [1 2; 3 4];
m1 = msh(V(1:2,:),[1,2]);
m2 = msh(V(3:4,:),[1,2]);

Vh1 = fem(m1,'P1');
Vh2 = fem(m2,'P1');

Gamma1 = dom(m1,7);
Gamma2 = dom(m2,7);

if (k == 0)
    Gxy = @(X,Y)(femGreenKernel(X,Y,'[log(r)]',0));
    c = -1/(2*pi);
else
    Gxy = @(X,Y)(femGreenKernel(X,Y,'[H0(kr)]',k));
    c = 1i/4;
end

A = c*integral(Gamma1,Gamma2,nxgrad(Vh1),Gxy,nxgrad(Vh2)) ...
    - k^2*c*integral(Gamma1,Gamma2,ntimes(Vh1),Gxy,ntimes(Vh2));
A = A + -1/(2*pi)*regularize(Gamma1,Gamma2,nxgrad(Vh1),'[log(r)]',nxgrad(Vh2))...
    - k^2*(-1)/(2*pi)*regularize(Gamma1,Gamma2,ntimes(Vh1),'[log(r)]',ntimes(Vh2));


end

