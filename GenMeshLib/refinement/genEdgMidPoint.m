function[Mfine,parentElt,parentVtx] = genEdgMidPoint(Mcoarse)

% old2new(i,:) = [k,l]: element number i gave rise to elements k and l,
% with k the "left half" and l the "right half"

Elt = Mcoarse.elt;
Nelt = size(Elt,1);

Nvtx = size(Mcoarse.vtx,1);
[edg2vtx,~,~] = subsimplices(Mcoarse,1);
Nedg = size(edg2vtx,1);
idx = [edg2vtx(:,1); edg2vtx(:,2)];
jdx = [edg2vtx(:,2); edg2vtx(:,1)];
v = [(1:size(edg2vtx,1))';(1:size(edg2vtx,1))'];

V = Mcoarse.vtx;
Nv = size(V,1);
edgMid = (V(edg2vtx(:,1),:) + V(edg2vtx(:,2),:))/2;
NewV = [V; edgMid]; % Set of vertices of the new mesh.

VTE = sparse(idx,jdx,v+Nv,Nvtx,Nvtx);
% VTE is a matrix of size Nvtx x Nvtx which maps a pair of indices of 
% coarse vertices and returns the index of their middle, if this middle
% exists.


ind = 1:Nelt;
NewElem = zeros(Nelt*2,2);

% We access VTE(i,j) via VTE((i-1)*Nvtx + j)
% "Left halves"
NewElem(1:Nelt,:) = [Elt(:,1) VTE((Elt(:,1)-1)*Nvtx + Elt(:,2))];
% "Right halves"
NewElem(Nelt + (1:Nelt),:) = [VTE((Elt(:,1)-1)*Nvtx + Elt(:,2)) Elt(:,2)];


neiElt = Mcoarse.nei_elt;
neiFct = Mcoarse.nei_fct;

NewNeiElt = zeros(2*Nelt,2);
NewNeiFct = zeros(2*Nelt,2);

% Newly created let halves are neighbors to corresponding right halves. 
NewNeiElt(1:Nelt,1) = Nelt + (1:Nelt);
NewNeiElt(Nelt + (1:Nelt),2) = 1:Nelt;
NewNeiFct(1:Nelt,1) = 2;
NewNeiFct(Nelt + (1:Nelt),2) = 1;

for el = 1:Nelt
    left_el = neiElt(el,2);
    right_el = neiElt(el,1);
    NewNeiFct(el,2) = neiFct(el,2);
    NewNeiFct(el+Nelt,1) = neiFct(el,1);
    
    if neiFct(el,2) == 1
        NewNeiElt(el,2) = left_el + Nelt;
    else
        NewNeiElt(el,2) = left_el;
    end
    if neiFct(el,1) == 1
        NewNeiElt(el+Nelt,1) = right_el + Nelt;
    else
        NewNeiElt(el+Nelt,1) = right_el;
    end
end

Mfine = GeneralizedMesh(NewV,NewElem,NewNeiElt,NewNeiFct);
parentElt = [1:Nelt,1:Nelt]';
parentVtx = [1:Nvtx,zeros(1,Nvtx)]';

end