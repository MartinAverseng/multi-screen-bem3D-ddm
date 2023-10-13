function [P,JumpOp,AvOp,F,gamma,I] = jumpSpaceP1(M)
% P defined by
% yk = sum_{(i,j)} P_{(i,j),k} phi_{(i,j)}
% J defined by
% [phi_{(i,j)}]_Gamma = sum_{k} J_{k,(i,j)} [y_k]_Gamma

% P is the matrix of the inclusion operator p : Psi_h -> \V_h
% J is the matrix of the projection j: \V_h -> Psi_h
% with respect to the direct sum \Vh = Psi_h + V_h

% One has the property j(p(y))  = y
% i.e. JP = Id_{QxQ}
% where Q is the dimension of the Psi_h




[F,gamma,I] = M.generalizedSubfacets(0);



Nvtx = max(F);


Nd = size(F,1);
memv = cell(Nvtx,1);
for gf = 1:Nd
    memv{F(gf)} = [memv{F(gf)};gf];
end

Ntilde = Nd - Nvtx;

nnzJ = 0;
for v = 1:Nvtx
    L = length(memv{v});
    nnzJ = nnzJ + L*(L-1);
end

nnzAv = 0;
for v = 1:Nvtx
    L = length(memv{v});
    nnzAv = nnzAv + L^2;
end

k = 0; % number of lines of P filled
l = 0;
tic;
iP = zeros(2*Ntilde,1);
jP = zeros(2*Ntilde,1);
vP = zeros(2*Ntilde,1);


for v = 1:Nvtx
    gfs = memv{v};
%     gfs = find(idx);
    for i = 1:(size(gfs)-1)
        iP([l+1,l+2],1) = [k+1;k+1];
        jP([l+1,l+2],1) = [gfs(i);gfs(end)];
        vP([l+1,l+2],1) = [1;-1];
        k = k+1;
        l = l+2;
    end
end






iJ = zeros(nnzJ,1);
jJ = zeros(nnzJ,1);
vJ = zeros(nnzJ,1);
% size of iJ is sum(L(L-1)) where L is the multiplicity of a genfacet



k = 0; % number of lines of JumpOp filled
l = 0;
for v = 1:Nvtx
    gfs = memv{v};
    L = length(gfs);
    for i = 1:L
        iJ((l+1):(l+(L-1)),1) = gfs(i)*ones(L-1,1);
        jJ((l+1):(l+(L-1)),1) = k+(1:(L-1))';
        vJtmp = ones(L-1,1)*(-1/L);
        if (i<=L-1)
            vJtmp(i) = 1 - 1/L;
        end
        vJ((l+1):(l+(L-1)),1) = vJtmp;
        l = l+L-1;
    end
    k = k+(L-1);
end


iAv = zeros(nnzAv,1);
jAv = zeros(nnzAv,1);
vAv = zeros(nnzAv,1);
k = 0;
l = 0;

for v = 1:Nvtx
    gfs = memv{v};
    L = length(gfs);
    for i = 1:L
        for j = 1:L
            iAv(l+1) = gfs(i);
            jAv(l+1) = k+j;
            vAv(l+1) = 1/L;
            l = l+1;
        end
    end
    k = k+L;
end




JumpOp = sparse(jJ,iJ,vJ,Ntilde,Nd);
P = sparse(jP,iP,vP,Nd,Ntilde);
AvOp = sparse(iAv,jAv,vAv,Nd,Nd);

end
