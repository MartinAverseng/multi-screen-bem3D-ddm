%% Comparison of "New" and "Exact" solutions (Fig 2.4)

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
Npoints = 20; % Number of points per branch. Ntot = 4*Npoints.



%% Conformal mapping U, and gradient of U


Z = @(X)(X(:,1) + 1i*X(:,2));
w = @(X)(w_conformal(Z(X)));
dw = @(X)(1 + Z(X)./(w_conformal(Z(X)) - Z(X)));
U = @(X)(real(-1./(2*1i*w(X))));
gradU{1} = @(X)(real(1./(2*1i*w(X).^2).*dw(X)));
gradU{2} = @(X)(real(1./(2*1i*w(X).^2).*1i.*dw(X)));
gradU{3} = @(X)(0*X(:,1));

%% Computing new approximations

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
rhs = rhs2dW(M,gradU);

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

% Compute exact solution

uexact = U(VTZ);
uexact = reshape(uexact,N,N);


%% Plotting exact solution

figure


[M,c] = contourf(x,y,uexact,30,'--');
c.LineWidth = 1;
axis equal
hold on
plot([-1,1],[0 0],'k','LineWidth',7)
plot([0,0],[-1 1],'k','LineWidth',7)
axis off
axis tight
set(gca,'FontSize',30);
colormap(jet)
[a,b] = caxis;
c = colorbar;
c.TickLabelInterpreter = 'latex';
% title("Exact solution",'Interpreter','latex')

%% Plotting "New" approximation

figure
[M,c] = contourf(x,y,unew,30,'--');
c.LineWidth = 1;
axis equal
hold on
plot([-1,1],[0 0],'k','LineWidth',7)
plot([0,0],[-1 1],'k','LineWidth',7)
% caxis([a,b]);
axis off
axis tight
set(gca,'FontSize',30);
caxis([a,b]);
colormap(jet)
colorbar
c = colorbar;
c.TickLabelInterpreter = 'latex';
% title("New approximation, $h = 0.1$",'Interpreter','latex')

%% Plotting error

figure

[M,c] = contourf(x,y,abs(unew - uexact),30,'--');
c.LineWidth = 1;
axis equal
hold on
plot([-1,1],[0 0],'k','LineWidth',7)
plot([0,0],[-1 1],'k','LineWidth',7)
% caxis([a,b]);
axis off
axis tight
set(gca,'FontSize',30);
caxis([1e-4,1]);
colorbar
c = colorbar;
c.TickLabelInterpreter = 'latex';
% title("Error for naive approximation, $h = 0.1$",'Interpreter','latex')

set(gca,'ColorScale','log')





