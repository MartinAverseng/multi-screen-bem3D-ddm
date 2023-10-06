% element-wise assembly of hypersingular operator on a generalized mesh.
function [A] = bemAssembly(M)
% R restriction matrix.

MatlabInd = 1; % Recall C++ is O-indexed and Matlab 1-indexed.
Ne = M.nelt;
[J,Nf] = localToGlobal(M,0);
J = J - MatlabInd;
Nvtx = size(M.vtx,1);
% A = zeros(Nf);
Melt = M.elt - MatlabInd;
% A = CppSubAssembly(Nf,Ne,J,M.vtx,Melt,A,Klist,rows,cols,vals,n); % -1 for 0-indexing in C++.
A = CppAssembly(Nf,Ne,Nvtx,J,M.vtx,Melt);

% n = M.n;
% Kd = n+1;
%
% A2 = zeros(Nf,Nf);
% for el1 = 1:Nelt
%     T1 = M.elt(el1,:);
%
%     for el2 = 1:Nelt
%         T2 = M.elt(el2,:);
%         [b,I] = ismember(T2,T1); % Identifying common vertices
%
%         %communicating with C++
%         t2 = [4 5 6];
%         t2(b) = I(b);
%         tri = [1 2 3; t2];
%         V = M.vtx([T1;T2]',:);
%         Al = singularIntegral(V,tri-1); % Retrieving local matrices
%
%         for k = 1:Kd % subsimplices of element el1
%             for l = 1:Kd % subsimplices of element el2
%                 genFk = J(el1,k);
%                 genFl = J(el2,l);
%                 A2(genFk,genFl) = A2(genFk,genFl) + Al(k,l);
%             end
%         end
%     end
%
% end
%
% norm(A1 - A2)
% A = (A + A')/2;


end

