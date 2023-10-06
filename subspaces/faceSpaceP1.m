function[S] = faceSpaceP1(MH,Mh,parentElt)


[l2gH,NH] = localToGlobal(MH,0);
[l2gh,Nh] = localToGlobal(Mh,0);

elt2 = MH.nelt/2;
S = cell(elt2,1);

% Loop over coarse elements
for kH = 1:elt2
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
    ind = abs(prod(bary,2)) > 1e-10; % Only keep vertices in the relative interior
    idx = unique(ts(ind));
    NH = length(idx);
    jdx = 1:NH;
    S{kH} = sparse(idx,jdx,1,Nh,NH);
end


end