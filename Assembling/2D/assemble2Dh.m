function [A] = assemble2Dh(M,waveNum)


[J,Nf] = localToGlobal(M,0);

Nelt = M.nelt;
n = M.n;
Kd = n+1;
V = M.vtx;
A = zeros(Nf,Nf);
for el1 = 1:Nelt
        T1 = M.elt(el1,:);
    
    for el2 = 1:Nelt
        T2 = M.elt(el2,:);
        xA = V(T1(1),1);
        yA = V(T1(1),2);
        xB = V(T1(2),1);
        yB = V(T1(2),2);
        xC = V(T2(1),1);
        yC = V(T2(1),2);
        xD = V(T2(2),1);
        yD = V(T2(2),2);
        Al = MatLoc2d(xA,yA,xB,yB,xC,yC,xD,yD,waveNum);
        
        for k = 1:Kd % subsimplices of element el1
            for l = 1:Kd % subsimplices of element el2
                genFk = J(el1,k);
                genFl = J(el2,l);
                A(genFk,genFl) = A(genFk,genFl) + Al(k,l);
            end
        end
    end
    
end



end

