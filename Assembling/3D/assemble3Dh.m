function [A] = assemble3Dh(M,waveNum)


[J,Nf] = localToGlobal(M,0);

Nelt = M.nelt;
n = M.n;
Kd = n+1;
V = M.vtx;
A = zeros(Nf,Nf);
for el1 = 1:Nelt
        T1 = M.vtx(M.elt(el1,:),:);
    
    for el2 = 1:Nelt
        T2 = M.vtx(M.elt(el2,:),:);
        
        Al = MatLoc3d(T1,T2,waveNum);
        
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

