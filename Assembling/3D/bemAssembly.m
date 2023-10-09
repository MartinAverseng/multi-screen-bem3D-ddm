% element-wise assembly of hypersingular operator on a generalized mesh.
function [A] = bemAssembly(M)
% R restriction matrix.

MatlabInd = 1; % Recall C++ is O-indexed and Matlab 1-indexed.
Ne = M.nelt;
[J,Nf] = localToGlobal(M,0);
J = J - MatlabInd;
Nvtx = size(M.vtx,1);
Melt = M.elt - MatlabInd;
A = CppAssembly(Nf,Ne,Nvtx,J,M.vtx,Melt);


end

