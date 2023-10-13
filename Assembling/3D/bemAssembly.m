% element-wise assembly of hypersingular operator on a generalized mesh.
function [A] = bemAssembly(M,k)
% R restriction matrix.

if ~exist('k','var')||isempty(k)
    k = 0;
end

MatlabInd = 1; % Recall C++ is O-indexed and Matlab 1-indexed.
Ne = M.nelt;
[J,Nf] = localToGlobal(M,0);
J = J - MatlabInd;
Nvtx = size(M.vtx,1);
Melt = M.elt - MatlabInd;
if k == 0
    A = CppAssembly(Nf,Ne,Nvtx,J,M.vtx,Melt);
else
    [Ar,Ai] = CppAssemblyHelmholtz(Nf,Ne,Nvtx,J,M.vtx,Melt,k);
    A = Ar + 1i*Ai;
end

end

