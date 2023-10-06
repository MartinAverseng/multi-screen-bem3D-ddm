function [r] = rhs2dW(M,func)



[J,Nf] = localToGlobal(M,0);
Nelt = M.nelt;
r = zeros(Nf,1);


for el = 1:Nelt
    genF1 = J(el,1);
    genF2 = J(el,2);
    V1 = M.vtx(M.elt(el,1),:);
    V2 = M.vtx(M.elt(el,2),:);
    E = [1,2];
    m = msh([V1;V2],E);
    Vh = fem(m,'P1');
    Gamma = dom(m,3);
    Iel = integral(Gamma,func,ntimes(Vh));
    
    r(genF1) = r(genF1) + Iel(1);
    r(genF2) = r(genF2) + Iel(2);
end




end

