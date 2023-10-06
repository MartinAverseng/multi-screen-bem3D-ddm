%% Comparison of "New" and "Exact" solutions with a 4-valued solution (Fig 2.5)

% run from the folder where the script is located
clear all; %#ok
close all;

% checking I'm in the right place
cd ..
cd FiguresScripts

% Adding library filed
cd ..
addpath(genpath(pwd));

% Coming back to the start
cd FiguresScripts

beta = 1;
Npoints = 10; % Number of points per branch. Ntot = 4*Npoints.


gradf{1} = @(X)(0*X(:,1) + 1);
gradf{2} = @(X)(0*X(:,2) + 2);
gradf{3} = @(X)(0*X(:,1));

%% Computing "New" approximations

h = 1/Npoints;
m1 = mshSegment(Npoints,1);
m1 = translate(m1,[0.5, 0 0]);
m1.vtx = 1 - (1 - m1.vtx).^beta;
m = m1;
m2 = m1;
for i = 1:3
    m2 = rotate(m2,[0 0 1],pi/2);
    m2 = swap(m2);
    m = union(m,m2);
end
m = mshClean(m,1e-5);


M = intrinsicInflation(m);

W = assemble2Dh(M,0);
[P,J] = jumpSpaceP1(M);
rhs = rhs2dW(M,gradf);

rhstilde = P'*rhs;
Wtilde = P'*W*P;
phiNew = P*(Wtilde\rhstilde);

N = 500;
x = linspace(-2,2,N);
y = linspace(-2,2,N);
[X,Y] = meshgrid(x,y);
VTZ = [X(:),Y(:),0*X(:)];

DL = DL2dP1(M,VTZ,0);
unew = -DL*phiNew;

unew = reshape(unew,N,N);



%% Computing naive approximations

h = 1/Npoints;
m1 = mshSegment(Npoints,1);
m1 = translate(m1,[0.5, 0 0]);
m1.vtx = 1 - (1 - m1.vtx).^beta;
m = m1;
m2 = m1;
for i = 1:3
    m2 = rotate(m2,[0 0 1],pi/2);
    m2 = swap(m2);
    m = union(m,m2);
end
m = mshClean(m,1e-5);


% Solve with Naive method
Vh = dirichlet(fem(m,'P1'),m.bnd);
Gamma = dom(m,3);
k = 0;
Gxy = @(X,Y)femGreenKernel(X,Y,'[log(r)]',k);
gradyGxy{1} = @(X,Y)(femGreenKernel(X,Y,'grady[log(r)]1',k));
gradyGxy{2} = @(X,Y)(femGreenKernel(X,Y,'grady[log(r)]2',k));
gradyGxy{3} = @(X,Y)(femGreenKernel(X,Y,'grady[log(r)]3',k));
W = -1/(2*pi)*integral(Gamma,Gamma,nxgrad(Vh),Gxy,nxgrad(Vh)) ...
    + (-1)/(2*pi)*regularize(Gamma,Gamma,nxgrad(Vh),'[log(r)]',nxgrad(Vh));

rhs = integral(Gamma,ntimes(Vh),gradf);
phiNaive = W\rhs;

N = 500;
x = linspace(-2,2,N);
y = linspace(-2,2,N);
[X,Y] = meshgrid(x,y);
VTZ = [X(:),Y(:),0*X(:)];

DL = -1/(2*pi)*integral(VTZ,Gamma,gradyGxy,ntimes(Vh),1e-2);
DL = DL + (-1)/(2*pi)*regularize(VTZ,Gamma,'grady[log(r)]',ntimes(Vh));
unaive = -DL*phiNaive;

unaive = reshape(unaive,N,N);


%% Plotting "New" approximation

figure
[M,c] = contourf(x,y,unew,30,'--');
c.LineWidth = 1;
axis equal
hold on
plot([-1,1],[0 0],'k','LineWidth',7)
plot([0,0],[-1 1],'k','LineWidth',7)
[a,b] =  caxis;
axis off
axis tight
set(gca,'FontSize',30);
colormap(jet)
colorbar
c = colorbar;
c.TickLabelInterpreter = 'latex';


%% Plotting "Naive" approximation

figure
[M,c] = contourf(x,y,unaive,30,'--');
c.LineWidth = 1;
axis equal
hold on
plot([-1,1],[0 0],'k','LineWidth',7)
plot([0,0],[-1 1],'k','LineWidth',7)
caxis([a,b]);
axis off
axis tight
set(gca,'FontSize',30);
colormap(jet)
colorbar
c = colorbar;
c.TickLabelInterpreter = 'latex';





