%% Condition numbers for the 2D algorithm (Fig. 2.6, bottom)
% Getting all results with NrefineMax = 9 takes about 8 hours.
% The assembing routine for the hypersingular is embarassingly non-optimized.


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

mexExec = true; % switch to true for mex execution

%% Geometry

V = [0 0 0; 1 0 0; cos(2*pi/3) sin(2*pi/3) 0; cos(4*pi/3) sin(4*pi/3) 0; 1+cos(pi/3) 0 sin(pi/3); 1+cos(pi/3) 0 -sin(pi/3)];
E = [1 2 3; 1 3 4; 1 4 2;2 1 5;2 1 6;5 2 6];
m = msh(V,E);
plot(m);
axis equal;
%% Main loop to compute condition numbers


NrefineMin = 0;
NrefineMax = 2;

condNoPrec = NaN + zeros(NrefineMax,1);
condPrec = NaN + zeros(NrefineMax);

for Nrefine = NrefineMin:NrefineMax % Assembling operator
    disp(Nrefine)
    M0 = intrinsicInflation(m);
    if Nrefine> 0
        M = M0.refine(Nrefine);
    else
        M = M0;
    end
    
    
    
    [Ph,Jh,Av,F,gamma,I] = jumpSpaceP1(M);
    if mexExec
        Wh = bemAssembly(M); 
    else
        Wh = slowBemAssembly(M);
    end
    
    
    Whtilde = Ph'*Wh*Ph;
    
    
    condNoPrec(Nrefine+1) = cond(Whtilde);
    
    for Nprec = 0:Nrefine % Assembling Preconditioner
        
        disp(Nprec)
        % Splitting
        
        if Nprec == 0
            Mcoarse = M0;
            [~,parentElt] = Mcoarse.refine(Nrefine);
        elseif Nprec < Nrefine
            Mcoarse = M0.refine(Nprec);
            [~,parentElt] = Mcoarse.refine(Nrefine-Nprec);
        else
            Mcoarse = M;
            parentElt = 1:M.nelt;
        end
        if Nrefine==0
            parentElt = 1:M.nelt;
        end
        [PH,JH] = jumpSpaceP1(Mcoarse);
        if mexExec
            WH = bemAssembly(Mcoarse); 
        else
            WH = slowBemAssembly(Mcoarse);
        end
        
        
        S = coarseSpaceP1(Mcoarse,M,parentElt); % S : {phi_{i,j}^H} -> {phi_{i,j}^h}
        F = faceSpaceP1(Mcoarse,M,parentElt);
        W = wireBasketSpaceP1(Mcoarse,M,parentElt);
        
        
        RS = Jh*S*PH; % {yij^H} -> {yij^h}
        RW = Jh*W;
        
        
        WS = PH'*WH*PH;
        Prec = RS*(WS\(RS'));
        for i = 1:length(F)
            RFi = Jh*F{i};
            R = Ph*RFi;
            if mexExec
                WFi = bemSubAssembly(M,R,gamma); 
            else
                WFi = R'*Wh*R;
            end
            Prec = Prec + RFi*(WFi\(RFi'));
        end
        
        %         WW = bemSubAssembly(M,Ph*RW,gamma);
        if mexExec
            WW = bemSubAssembly(M,Ph*RW,gamma); 
        else
            WW = RW'*Whtilde*RW;    
        end
        
        Prec = Prec + RW*(WW\(RW'));
        
        [~,D] = eig(Prec*Whtilde);
        d = sort(real((diag(D))));
        kappa = max(d)/min(d);
        condPrec(Nrefine+1,Nprec+1) = kappa;
    end
    
end
% 
% save('condNoPrec','condNoPrec')
% save('condPrec','condPrec');

%%
close all;
% 
% load('condNoPrec','condNoPrec');
% load('condPrec','condPrec');

h = 1./2.^(0:NrefineMax);
level = 1+log(1./h)/log(2);


figure;
loglog(level,condNoPrec,'-o','LineWidth',2);
hold on
for j = 1:min(5,NrefineMax)
    loglog(level,condPrec(:,j),'-o','LineWidth',2);
end
loglog(level,0.7 + 0.7*log(1./h),'k--','LineWidth',2);
loglog(level,0.7 + 0.7*log(1./h).^2,'k--','LineWidth',2);
loglog(level,0.7*1./h,'k--','LineWidth',2);

annotation('textarrow',[0.77,0.7],[0.8 0.8],'String','No preconditioner','Interpreter','latex','Color',[0 0.4470 0.7410],'FontSize',20)
annotation('textarrow',[0.77,0.7],[0.55 0.55],'String','$H = 1$ ','Interpreter','latex','Color',[0.8500 0.3250 0.0980],'FontSize',20)
annotation('textarrow',[0.77,0.7],[0.52 0.52],'String','$H = \frac12$ ','Interpreter','latex','Color',[0.9290 0.6940 0.1250],'FontSize',20)
annotation('textarrow',[0.77,0.7],[0.47 0.47],'String','$H = \frac14$ ','Interpreter','latex','Color',[0.4940 0.1840 0.5560],'FontSize',20)
annotation('textarrow',[0.77,0.7],[0.4 0.4],'String','$H = \frac18$ ','Interpreter','latex','Color',[0.4660 0.6740 0.1880],'FontSize',20)
annotation('textarrow',[0.77,0.7],[0.3 0.3],'String','$H = \frac1{16}$ ','Interpreter','latex','Color',[0.3010 0.7450 0.9330],'FontSize',20)
% ylabel()
axis tight
xlim([1 15]);
ylim([1, 100]);
xlabel("Refinement level: $1+\log_2(1/h)$","Interpreter","latex");
ylabel("Condition number $\kappa$","Interpreter","latex");
set(gca,'FontSize',30);
% axis tight

xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex'; % latex for x-axis
yaxisproperties= get(gca, 'YAXIS');
yaxisproperties.TickLabelInterpreter = 'latex'; % latex for y-axis
set(gca,'XTick',[1 2 3 4 5 6 7]);
