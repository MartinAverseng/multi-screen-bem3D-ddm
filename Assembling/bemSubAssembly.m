% element-wise assembly of hypersingular operator on a generalized mesh.
function [A] = bemSubAssembly(M,R,gamma)
% R restriction matrix.
% 
MatlabInd = 1; % Recall C++ is O-indexed and Matlab 1-indexed.
[J,Nf] = localToGlobal(M,0);
J = J - MatlabInd;
[rows,~,~] = find(R);
dofList = rows; % All dofs i,j for which we need the value of a([phi_i],[phi_j])
Nsub = length(dofList); % Submatrix that we compute using the sub-assembly
dictionnary = zeros(Nf,1);
dictionnary(dofList) = 1:length(dofList);% dictionnary(dofList(i)) = i. 
Isigma = sparse(1:Nsub,dofList,1,Nsub,Nf);
% dictionnary(dof) = 0 if dof not in the list -> we don't store the value. 
Klist = incidentElements(gamma,dofList) - MatlabInd; % All elements which contribute to those values 


Nsubelt = size(Klist,1);
Ne = M.nelt;
Nvtx = size(M.vtx,1);
Melt = M.elt - MatlabInd;
CppDic = dictionnary - MatlabInd;
Asub = CppSubAssembly(Nsub,Nsubelt,Nf,Ne,Nvtx,J,M.vtx,Melt,Klist,CppDic);
A = R'*(Isigma'*(Asub*(Isigma*R)));


end

