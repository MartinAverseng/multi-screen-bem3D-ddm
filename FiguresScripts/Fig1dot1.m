close all;
clear all; %#ok
cd ..
cd FiguresScripts
cd ..
addpath(genpath(pwd));
cd FiguresScripts

k = 20;

m = msh('junction2.msh');
m.vtx = 0.5*m.vtx;
m.vtx(:,1) = m.vtx(:,1)-0.25;
plot(m);


M = intrinsicInflation(m);
M = M.refine(2);
[Ph,Jh,Av,F,gamma,I] = jumpSpaceP1(M);
Wh = bemAssembly(M,k);

Whtilde = Ph'*Wh*Ph;


u(1,1) = -1/sqrt(3);
u(2,1) = 1/sqrt(3);
u(3,1) = -1/sqrt(3);

f = @(X)(exp(1i*k*X*u));
gradf = cell(3,1);
for j = 1:3
    gradf{j} = @(X)(1i*k*u(j)*f(X));
end
rhs = rhs3dW(M,gradf);

rhstilde = Ph'*rhs;

phi1 = Ph*(Whtilde\rhstilde);


%%
figure;
mSquare = mshSquare(20000,[10 10]);
mSquare.vtx(:,3) = -5;


DL = DL3dP1_vectorized(M,mSquare.vtx,k);
u1 = f(mSquare.vtx)+DL*phi1;
plot(mSquare,abs(u1).^2);
hold on


%%
mSquare = mshSquare(20000,[10 10]);
mSquare.vtx(:,3) = mSquare.vtx(:,2);
mSquare.vtx(:,2) = 5;

DL = DL3dP1_vectorized(M,mSquare.vtx,k);
u1 = f(mSquare.vtx)+DL*phi1;
plot(mSquare,abs(u1).^2);

%%
mSquare = mshSquare(20000,[10 10]);
mSquare.vtx(:,3) = mSquare.vtx(:,1);
mSquare.vtx(:,1) = -5;

DL = DL3dP1_vectorized(M,mSquare.vtx,k);
u1 = f(mSquare.vtx)+DL*phi1;
plot(mSquare,abs(u1).^2);
cc = caxis;
plot(m);
caxis(cc);
axis equal;
view(3)
arrow3([0 0 0]-4*u',[0 0 0] - 3*u');
plot(m.bnd,'r');




