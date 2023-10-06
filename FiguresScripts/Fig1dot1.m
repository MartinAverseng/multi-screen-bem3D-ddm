%% Sound-hard wave scattering by threefold junction (Fig 1.1, right panel)

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

%% Parameters:

Nrefine = 2; % Use value smaller than 3 to get somewhat quick results. 
thetaInc = pi/4; % Angle of incidence (pi/5 in the paper)

%% Geometry of the screen specified as a simple mesh. 

V = [0 0 0; 1 0 0; -cos(pi/3),sin(pi/3),0; -cos(pi/3) -sin(pi/3) 0];
E = [1 2; 1 3; 1 4];
m = msh(V,E);
plot(m);
axis equal;
title("The threefold junction")

%% Automatic generation of the inflated mesh

M = intrinsicInflation(m);
M = M.refine(Nrefine);
h = max(M.ndv);
hold on
plot(M.vtx(:,1),M.vtx(:,2),'ro');


%% Assembling of Helmholtz hypersingular operator

k = 0.7/h; % Pick a frequency resolved by the mesh. 
tic
W = assemble2Dh(M,k);
toc
% non-optimized routine for 2D assembling
% This can take some time (please refer to the 3D code for faster assembling). 
[P,J] = jumpSpaceP1(M); 
% J: jump operator multi-traces -> jumps
% P: extension operator : jumps -> multi-traces (choice of a
% representative)
disp('Done assembling')

%% Solving

f = @(X)(exp(1i*k*(cos(thetaInc)*X(:,1) + sin(thetaInc)*X(:,2))));
gradf{1} = @(X)(1i*k*cos(thetaInc)*f(X));
gradf{2} = @(X)(1i*k*sin(thetaInc)*f(X));
gradf{3} = @(X)(0*f(X));

rhs = rhs2dW(M,gradf);

rhstilde = P'*rhs;
Wtilde = P'*W*P;
phi = P*(Wtilde\rhstilde);
disp('Done solving');

%% Plotting 

% We plot separately in the 3 sectors 
% {r e^{i theta} : 2jpi/3 <= theta <= 2(j+1)pi/3}, j = 0,1,2
figure
side = 4;
mSquare = mshSquare(5000,[side,side]);

% First sector

mSquare1 = mSquare;
mSquare.vtx(:,1) = mSquare.vtx(:,1) + side/2;
mSquare.vtx(:,2) = mSquare.vtx(:,2) + side/2;
Z = mSquare.vtx(:,1) + 1i*mSquare.vtx(:,2);
Z = Z.^(4/3); % transformation into the first sector
mSquare1.vtx(:,1) = real(Z);
mSquare1.vtx(:,2) = imag(Z);


DL = DL2dP1(M,mSquare1.vtx,k); % Helmholtz Double-layer potential
u = f(mSquare1.vtx) + DL*phi;
plot(mSquare1,real(u));


% Second sector
mSquare2 = mSquare;
Z = exp(1i*2*pi/3)*Z;
mSquare2.vtx(:,1) = real(Z);
mSquare2.vtx(:,2) = imag(Z);
DL = DL2dP1(M,mSquare2.vtx,k);
u = f(mSquare2.vtx) + DL*phi;
hold on
plot(mSquare2,real(u));



% Third sector
mSquare3 = mSquare;
Z = exp(1i*2*pi/3)*Z;
mSquare3.vtx(:,1) = real(Z);
mSquare3.vtx(:,2) = imag(Z);
DL = DL2dP1(M,mSquare3.vtx,k);
u = f(mSquare3.vtx) + DL*phi;
plot(mSquare3,real(u));

% Overlay junction mesh
plot(m);
axis equal;
axis off
plot([0,1],[0,0],'k-');
plot([0,cos(2*pi/3)],[0,sin(2*pi/3)],'k-');
plot([0,cos(4*pi/3)],[0,sin(4*pi/3)],'k-');
title("Wave-Scattering by the threefold junction")


