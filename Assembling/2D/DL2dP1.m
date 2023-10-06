function [r] = DL2dP1(M,X,k)


[J,Nf] = localToGlobal(M,0);
Nelt = M.nelt;
r = zeros(size(X,1),Nf);

if (k == 0)
    Gxy = @(X,Y)(femGreenKernel(X,Y,'[log(r)]',0));
    gradyGxy{1} = @(X,Y)(femGreenKernel(X,Y,'grady[log(r)]1',0));
    gradyGxy{2} = @(X,Y)(femGreenKernel(X,Y,'grady[log(r)]2',0));
    gradyGxy{3} = @(X,Y)(femGreenKernel(X,Y,'grady[log(r)]3',0));
    c = -1/(2*pi);
else
    Gxy = @(X,Y)(femGreenKernel(X,Y,'[H0(kr)]',k));
    gradyGxy{1} = @(X,Y)(femGreenKernel(X,Y,'grady[H0(kr)]1',k));
    gradyGxy{2} = @(X,Y)(femGreenKernel(X,Y,'grady[H0(kr)]2',k));
    gradyGxy{3} = @(X,Y)(femGreenKernel(X,Y,'grady[H0(kr)]3',k));
    c = 1i/4;
end


for el = 1:Nelt
    genF1 = J(el,1);
    genF2 = J(el,2);
    V1 = M.vtx(M.elt(el,1),:);
    V2 = M.vtx(M.elt(el,2),:);
    E = [1,2];
    m = msh([V1;V2],E);
    Vh = fem(m,'P1');
    Gamma = dom(m,3);
    DLel = c*integral(X,Gamma,gradyGxy,ntimes(Vh));
    DLel = DLel + (-1)/(2*pi)*regularize(X,Gamma,'grady[log(r)]',ntimes(Vh));
    
    r(:,genF1) = r(:,genF1) + DLel(:,1);
    r(:,genF2) = r(:,genF2) + DLel(:,2);
end



end

