%% Comparison of "Naive" and "Exact" solution U (Fig 2.3)

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
Npoints = 40; % Number of points per branch. Ntot = 4*Npoints.



%% Conformal mapping U, and gradient of U


Z = @(X)(X(:,1) + 1i*X(:,2));
w = @(X)(w_conformal(Z(X)));
dw = @(X)(1 + Z(X)./(w_conformal(Z(X)) - Z(X)));
U = @(X)(real(-1./(2*1i*w(X))));
gradU{1} = @(X)(real(1./(2*1i*w(X).^2).*dw(X)));
gradU{2} = @(X)(real(1./(2*1i*w(X).^2).*1i.*dw(X)));
gradU{3} = @(X)(0*X(:,1));

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

rhs = integral(Gamma,ntimes(Vh),gradU);
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

%% Plotting "Naive" approximation

figure
[M,c] = contourf(x,y,unaive,30,'--');
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
% title("Naive approximation, $h = 0.025$",'Interpreter','latex')

%% Plotting error

figure

[M,c] = contourf(x,y,abs(unaive - uexact),30,'--');
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
% title("Error for naive approximation, $h = 0.025$",'Interpreter','latex')

set(gca,'ColorScale','log')





