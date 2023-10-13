function [A] = MatLoc3d(T1,T2,k)


m1 = msh(T1,[1,2,3]);
m2 = msh(T2,[1,2,3]);

Vh1 = fem(m1,'P1');
Vh2 = fem(m2,'P1');

Gamma1 = dom(m1,12);
Gamma2 = dom(m2,12);

if (k == 0)
    Gxy = @(X,Y)(femGreenKernel(X,Y,'[1/r]',0));
    c = 1/(4*pi);
else
    Gxy = @(X,Y)(femGreenKernel(X,Y,'[exp(ikr)/r]',k));
    c = 1/(4*pi);
end

A = c*integral(Gamma1,Gamma2,nxgrad(Vh1),Gxy,nxgrad(Vh2)) ...
    - k^2*c*integral(Gamma1,Gamma2,ntimes(Vh1),Gxy,ntimes(Vh2));
A = A + c*regularize(Gamma1,Gamma2,nxgrad(Vh1),'[1/r]',nxgrad(Vh2))...
    - k^2*c*regularize(Gamma1,Gamma2,ntimes(Vh1),'[1/r]',ntimes(Vh2));


end

