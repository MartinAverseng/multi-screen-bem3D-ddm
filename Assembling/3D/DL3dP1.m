function [r] = DL3dP1(M,X,k)


[J,Nf] = localToGlobal(M,0);
Nelt = M.nelt;
r = zeros(size(X,1),Nf);
c = 1/(4*pi);
if (k == 0)
%     Gxy = @(X,Y)(femGreenKernel(X,Y,'[1/r]',0));
    gradyGxy{1} = @(X,Y)(femGreenKernel(X,Y,'grady[1/r]1',0));
    gradyGxy{2} = @(X,Y)(femGreenKernel(X,Y,'grady[1/r]2',0));
    gradyGxy{3} = @(X,Y)(femGreenKernel(X,Y,'grady[1/r]3',0));
else
%     Gxy = @(X,Y)(femGreenKernel(X,Y,'[exp(ikr)/r]',k));
    gradyGxy{1} = @(X,Y)(femGreenKernel(X,Y,'grady[exp(ikr)/r]1',k));
    gradyGxy{2} = @(X,Y)(femGreenKernel(X,Y,'grady[exp(ikr)/r]2',k));
    gradyGxy{3} = @(X,Y)(femGreenKernel(X,Y,'grady[exp(ikr)/r]3',k));
end







for el = 1:Nelt
    genF1 = J(el,1);
    genF2 = J(el,2);
    genF3 = J(el,3);
    V1 = M.vtx(M.elt(el,1),:);
    V2 = M.vtx(M.elt(el,2),:);
    V3 = M.vtx(M.elt(el,3),:);
    E = [1,2,3];
    m = msh([V1;V2;V3],E);
    Vh = fem(m,'P1');
    Gamma = dom(m,7);
    DLel = c*integral(X,Gamma,gradyGxy,ntimes(Vh));
    DLel = DLel + c*regularize(X,Gamma,'grady[1/r]',ntimes(Vh));
    
    r(:,genF1) = r(:,genF1) + DLel(:,1);
    r(:,genF2) = r(:,genF2) + DLel(:,2);
    r(:,genF3) = r(:,genF3) + DLel(:,3);
end



end

