function[S] = coarseSpaceP1(MH,Mh,parentElt)

% Returns the matrix S defined by
% phi_{H,(i,j)} = \sum_{(i',j')} S((i',j'),(i,j)) phi_{h,(i',j')}
% i.e. S is the matrix of the embedding
% V_H(Omega \ Gamma) -> V_h(Omega \ Gamma)


[l2gH,NH] = localToGlobal(MH,0);
[l2gh,Nh] = localToGlobal(Mh,0);



% Loop over coarse elements
S = zeros(Nh,NH);
for kH = 1:MH.nelt
    VH = MH.vtx(MH.elt(kH,:),:); % Coordinates of the simplex
    khs = (parentElt==kH);
    l = sum(khs);
    ts = zeros((Mh.n+1)*l,1);
    Vhs = zeros((Mh.n+1)*l,3);
    for sigma = 1:Mh.n+1
        ind = (sigma-1)*l + (1:l);
        ts(ind) = l2gh(khs,sigma);
        Vhs(ind,:) = Mh.vtx(Mh.elt(khs,sigma),:);
    end
    bary = barycentricCoords(VH,Vhs);
    for beta= 1:(MH.n+1)
        s = l2gH(kH,beta);
        S(ts,s) = bary(:,beta);
    end
end




end