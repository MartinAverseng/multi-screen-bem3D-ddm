%% Condition numbers for the 2D algorithm (Fig. 2.6, top)
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


%% Geometry

V = [0 0 0; 1 0 0; 2 0 0];%; -cos(pi/3),sin(pi/3),0; -cos(pi/3) -sin(pi/3) 0];
E = [1 2; 2 3]; %1 3; 1 4];
m = msh(V,E);
plot(m);
axis equal;

%% Main loop to compute condition numbers


NrefineMin = 0;
NrefineMax = 8;

condNoPrec = NaN + zeros(NrefineMax,1);
condPrec = NaN + zeros(NrefineMax);

for Nrefine = NrefineMin:NrefineMax-1 % Assembling operator
    disp(Nrefine)
    M0 = intrinsicInflation(m);
    if Nrefine> 0
        M = M0.refine(Nrefine);
    else
        M = M0;
    end
    
    
    
    Wh = assemble2Dh(M,0);
    
    [Ph,Jh] = jumpSpaceP1(M);
    Whtilde = Ph'*Wh*Ph;
    
    
    condNoPrec(Nrefine+1) = cond(Whtilde);
    
    for Nprec = 0:Nrefine % Assembling Preconditioner
        
        disp(Nprec)
        % Splitting
        
        if Nprec == 0
            if Nrefine == 0
                Mcoarse = M;
                parentElt = 1:M.nelt;
            else
                
                Mcoarse = M0;
                [~,parentElt] = Mcoarse.refine(Nrefine);
            end
        elseif Nprec < Nrefine
            Mcoarse = M0.refine(Nprec);
            [~,parentElt] = Mcoarse.refine(Nrefine-Nprec);
        else
            Mcoarse = M;
            parentElt = 1:M.nelt;
        end
        [PH,~] = jumpSpaceP1(Mcoarse);
        
        S = coarseSpaceP1(Mcoarse,M,parentElt);
        F = faceSpaceP1(Mcoarse,M,parentElt);
        W = wireBasketSpaceP1(Mcoarse,M,parentElt);
        
        
        
        RS = Jh*S*PH; % removing kernel of jump operator
        WS = RS'*Whtilde*RS; % Local solver on coarse space
        
        RW = Jh*W;
        WW = RW'*Whtilde*RW;
        
        Prec = RS*inv(WS)*RS'; %#ok
        for i = 1:length(F)
            RFi = Jh*F{i};
            WFi = RFi'*Whtilde*RFi;
            Prec = Prec + RFi*inv(WFi)*RFi'; %#ok
            if sum(sum(isnan(Prec))) > 0
                1;
            else
            end
        end
        Prec = Prec + RW*inv(WW)*RW'; %#ok
        
        [~,D] = eig(Prec*Whtilde);
        d = sort((diag(D)));
        kappa = max(d)/min(d);
        condPrec(Nrefine+1,Nprec+1) = kappa;
    end
    
    
end


%%
close all;
%

h = 1./2.^(0:NrefineMax-1);
level = 1+log(1./h)/log(2);

figure;
loglog(level,condNoPrec,'-o','LineWidth',2);
hold on
for j = 1:min(5,NrefineMax-1)
    loglog(level,condPrec(:,j),'-o','LineWidth',2);
end
loglog(level,0.7 + 0.2*log(1./h),'k--','LineWidth',2);
loglog(level,0.7 + 0.2*log(1./h).^2,'k--','LineWidth',2);
loglog(level,1./h,'k--','LineWidth',2);

annotation('textarrow',[0.77,0.7],[0.8 0.8],'String','No preconditioner','Interpreter','latex','Color',[0 0.4470 0.7410],'FontSize',20)
annotation('textarrow',[0.77,0.7],[0.55 0.55],'String','$H = 1$ ','Interpreter','latex','Color',[0.8500 0.3250 0.0980],'FontSize',20)
annotation('textarrow',[0.77,0.7],[0.52 0.52],'String','$H = \frac12$ ','Interpreter','latex','Color',[0.9290 0.6940 0.1250],'FontSize',20)
annotation('textarrow',[0.77,0.7],[0.47 0.47],'String','$H = \frac14$ ','Interpreter','latex','Color',[0.4940 0.1840 0.5560],'FontSize',20)
annotation('textarrow',[0.77,0.7],[0.4 0.4],'String','$H = \frac18$ ','Interpreter','latex','Color',[0.4660 0.6740 0.1880],'FontSize',20)
annotation('textarrow',[0.77,0.7],[0.3 0.3],'String','$H = \frac1{16}$ ','Interpreter','latex','Color',[0.3010 0.7450 0.9330],'FontSize',20)
% ylabel()
axis tight
xlim([1 18]);
ylim([1, 400]);
xlabel("Refinement level: $1+\log_2(1/h)$","Interpreter","latex");
ylabel("Condition number $\kappa$","Interpreter","latex");
set(gca,'FontSize',30);
% axis tight

xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex'; % latex for x-axis
yaxisproperties= get(gca, 'YAXIS');
yaxisproperties.TickLabelInterpreter = 'latex'; % latex for y-axis
set(gca,'XTick',1:NrefineMax+1);
