function[S,Mult] = wireBasketSpaceP1(MH,Mh,parentElt)


[l2gh,Nh] = localToGlobal(Mh,0);

idx = [];
mult = [];


[F,~,~] = Mh.generalizedSubfacets(0);
A = accumarray(F,1);
mult = A(F);

lastMult = 0*mult;
i = 1;

while i <= length(lastMult)
    j = mult(i);
    lastMult(i+j - 1) = 1;
    i = i + j;
end
    

toKeep = lastMult==0;


% Loop over coarse elements
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
    ind = abs(prod(bary,2)) < 1e-10; % Only keep vertices on the relative boundary
    idx = [idx;ts(ind)];
end
idx = unique(idx);
idx = idx(toKeep(idx));
NW = length(idx); 
S = sparse(idx,[1:NW]',1,Nh,NW);
Mult = sparse(idx,[1:NW],mult(idx),Nh,NW);

end