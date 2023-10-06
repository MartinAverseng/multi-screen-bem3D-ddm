mOmega = msh('mesh10.msh');
[F,~,M] = subsetsRows(mOmega.elt,3);
A = mOmega.col.*M;
A = sort(A,1);
L1 = A(end-1,:);
L2 = A(end,:);
ind = and(L1~=0,L2-L1~=0);
ind = find(and(ind,or(L1~=1,L2~=2)));
interface = clean(msh(mOmega.vtx,F(ind,:)),1e-10);
ball = sum(interface.ctr.^2,2)<0.5;
mGamma = interface.sub(ball);



