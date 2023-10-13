function [r] = DL3dP1_vectorized(M,X,k)


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

Gamma = dom(M,3);
Vh = GenFem(M,'P1');
r = c*integral(X,Gamma,gradyGxy,ntimes(Vh));
r = r+c*regularize(X,Gamma,'grady[1/r]',ntimes(Vh));






end
