%% Comparison of "Naive" and "New" Galerkin methods (Fig 2.2, right panel)

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

NpointsTab = [5 10 20 40 80];
betaTab = [1,2]; % Grading parameter

% Domain where we compute the error in each of the N sectors
Nsquare = 2000;
mSquare = mshSquare(Nsquare,[2 2]);
mSquare = translate(mSquare,[1 1 0]);

%% Conformal mapping U, and gradient of U


Z = @(X)(X(:,1) + 1i*X(:,2));
w = @(X)(w_conformal(Z(X)));
dw = @(X)(1 + Z(X)./(w_conformal(Z(X)) - Z(X)));
U = @(X)(real(-1./(2*1i*w(X))));
gradU{1} = @(X)(real(1./(2*1i*w(X).^2).*dw(X)));
gradU{2} = @(X)(real(1./(2*1i*w(X).^2).*1i.*dw(X)));
gradU{3} = @(X)(0*X(:,1));

%% Calculating new and naive approximations
% This takes a very long runtime (need to optimize the BEM assembling code 
% for the "New" algorithm in 2D). 

for np = 1:length(NpointsTab)
    for beta_num = 1:length(betaTab)
        beta = betaTab(beta_num);
        % Create mesh
        Npoints = NpointsTab(np);
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
        phi1 = W\rhs;
        
        
        err2 = 0;
        for sector = 0:3
            mz = mSquare;
            X = mz.vtx(:,1);
            Y = mz.vtx(:,2);
            Z = X+1i*Y;
            Z = exp(1i*2*pi*(sector/4)).*Z; % Rotate to go to next sector
            Z = Z +	0.01*exp(1i*pi/4)*exp(1i*2*pi*(sector/4)); % avoid the boundary
            VTZ = [real(Z) imag(Z) 0*real(Z)];
            mz.vtx = VTZ;
            DL = -1/(2*pi)*integral(VTZ,Gamma,gradyGxy,ntimes(Vh));
            DL = DL + (-1)/(2*pi)*regularize(VTZ,Gamma,'grady[log(r)]',ntimes(Vh));
            uNaive = -DL*phi1;
            err2 = err2 + 1/Nsquare*sum(abs(uNaive - U(VTZ)).^2);
        end
        errNaive(np,beta_num) = sqrt(err2);
        
        % solve with New method
        
        
        M = intrinsicInflation(m);
        
        W = assemble2Dh(M,k);
        [P,J] = jumpSpaceP1(M);
        
        
        rhs = rhs2dW(M,gradU);
        
        rhstilde = P'*rhs;
        Wtilde = P'*W*P;
        phi = P*(Wtilde\rhstilde);
        
        err2 = 0;
        for sector = 0:3
            mz = mSquare;
            X = mz.vtx(:,1);
            Y = mz.vtx(:,2);
            Z = X+1i*Y;
            Z = exp(1i*2*pi*sector/4).*Z;
            Z = Z +	0.01*exp(1i*pi/4)*exp(1i*2*pi*(sector/4));
            VTZ = [real(Z) imag(Z) 0*real(Z)];
            DL = DL2dP1(M,VTZ,k);
            uNew = -DL*phi;
            err2 = err2 + 1/Nsquare*sum(abs(uNew - U(VTZ)).^2);
        end
        errNew(np,beta_num) = sqrt(err2);
        
        
    end
    
end
% save errNew and errNaive to avoid recomputing them. 
save('errNew','errNew');
save('errNaive','errNaive');

%% Plotting

% load errNew and errNaive to avoid recomputing them. 
load('errNew','errNew');
load('errNaive','errNaive');

htab = 1./NpointsTab;
figure
loglog(htab(1:5),errNaive(1:5,1),'bo-','DisplayName','``Naive", uniform','LineWidth',2);
hold on
loglog(htab(1:5),errNaive(1:5,2),'bx--','DisplayName','``Naive", non-uniform','LineWidth',2);

loglog(htab(1:5),errNew(1:5,1),'ro-','DisplayName','``New", uniform','LineWidth',2);
loglog(htab(1:5),errNew(1:5,2),'rx--','DisplayName','``New", non-uniform','LineWidth',2);
loglog(htab(1:5),0.15*htab.^(1),'k-','HandleVisibility','off','LineWidth',2);
loglog(htab(1:5),0.15*htab.^(2),'k--','HandleVisibility','off','LineWidth',2);
xlabel("$h$",'Interpreter','latex');
ylabel('Discrete $\ell^2$ error on $[-2,2]^2$','Interpreter','latex');
axis equal
axis tight
ylim([1e-5,0.2])
legend show
legend boxoff
set(gca,'FontSize',30)
set(legend,'Interpreter','latex');

xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex'; % latex for x-axis
yaxisproperties= get(gca, 'YAXIS');
yaxisproperties.TickLabelInterpreter = 'latex'; % latex for y-axis
set(gca,'XTick',[0.0125 0.025 0.05 0.1 0.2]);


