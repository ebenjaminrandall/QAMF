% Runs metabolite code

tspan = 0:.001:10; 

printon = 0; 

%% Constants

Ve = 2; 
Vx = 1; 

R = 8.3145;
T = 310.15;
F = 96485.33; %C/mole

pH_outer = 7.2; 
pH_inner = 7.4; 

Mg  = 1e-3; %M 
K   = 100e-3;  %M 

psi = 180 / 1000; %V
Cr_tot = 40e-3; 

consts = [Ve; Vx; 
    R; T; F; 
    pH_outer; pH_inner; 
    Mg; K; 
    psi; Cr_tot;
    ]; 

%% Parameters

% Dissociation constants
K_HATP = 2.757e-7;
K_KATP = 9.809e-2;
K_MATP = 8.430e-5;
K_HADP = 4.106e-7;
K_KADP = 1.319e-1;
K_MADP = 7.149e-4;
K_HPI  = 2.308e-7;
K_KPI  = 3.803e-1;
K_MPI  = 2.815e-2;

% J_F1F0
X_F1F0 = 812; 
n_H = 8/3;
DGr0_F1F0 = -4510;

% J_ANT 
E_ANT = 0.141;           %mol (L mito)^(-1)
k2_ANT_o = 9.54 / 60;    % s^(-1) 
k3_ANT_o = 30.05/60;     % s^(-1)
Ko_Do = 38.89e-6;    % microM
Ko_To = 56.05e-6;    %microM
alpha_1 = 0.2829;
alpha_2 = -0.2086;
alpha_3 = 0.2372;
delta_T = 0.0699; %Metelkin et al 2006
delta_D = 0.0167;

% J_ATPase 
X_ATPase = 0.01 / 10 / 4;
K_iADP = 2.41e-4;
DGr0_ATPase = 4510;

% J_PiC 
X_PiC = 3.34e7; %mol/s/M/lmito
k_PiC = 1.61; %mM

% J_CK
DGr0_CK = -7910; %j/mol
X_CK = 1e6; %M/s
kref = 7.14e8; % Wu et al 2008

pars = [K_HATP; K_KATP; K_MATP;
    K_HADP; K_KADP; K_MADP;
    K_HPI; K_KPI; K_MPI;
    X_F1F0; n_H; DGr0_F1F0;
    E_ANT; k2_ANT_o; k3_ANT_o; Ko_Do; Ko_To; 
    alpha_1; alpha_2; alpha_3; delta_T; delta_D;
    X_ATPase; K_iADP; DGr0_ATPase; 
    X_PiC; k_PiC; 
    DGr0_CK; X_CK; kref;
    ];

%% Initial conditions

% Concentrations (M)
ATP_x = 1e-3; 
ADP_x = 9e-3; 
Pi_x  = 1e-3; 
ATP_e = 7.5e-3; 
ADP_e = 0.05e-3; 
Pi_e  = 0.3e-3; 
Cr    = 0.6 * Cr_tot; 

x0 = [ATP_x; ADP_x; Pi_x; ATP_e; ADP_e; Pi_e; Cr]; 

%% Solve model 

options = odeset('RelTol',1e-8','AbsTol',1e-8); 
sol = ode15s(@model,tspan,x0,options,pars,consts); 
sols = deval(sol,tspan); 

ATP_x = sols(1,:) * 1e3; 
ADP_x = sols(2,:) * 1e3; 
Pi_x  = sols(3,:) * 1e3; 
ATP_e = sols(4,:) * 1e3; 
ADP_e = sols(5,:) * 1e3; 
Pi_e  = sols(6,:) * 1e3; 
Cr    = sols(7,:) * 1e3; 
CrP   = Cr_tot*1e3 - Cr; 

%% Plot 

hfig1 = figure(1);
clf
hold on 
plot(tspan,ATP_x,'Color',[0 1 .5],'linewidth',2)
plot(tspan,ADP_x,'Color',[.3 .6 .5],'linewidth',2)
plot(tspan,ATP_e,'Color',[.6 .3 .5],'linewidth',2)
plot(tspan,ADP_e,'Color',[1 0 .5],'linewidth',2)
legend('ATP_x','ADP_x','ATP_e','ADP_e')
set(gca,'FontSize',20)
xlim([tspan(1) tspan(end)])
ylim([-.01 10.01])
xlabel('Time (s)')
ylabel('Concentration (mM)')

hfig2 = figure(2);
clf
hold on 
plot(tspan,Pi_x, 'b','linewidth',2)
plot(tspan,Pi_e, 'r','linewidth',2)
legend('Pi_x','Pi_e')
set(gca,'FontSize',20)
xlim([tspan(1) tspan(end)])
ylim([-.01 10.01])
xlabel('Time (s)')
ylabel('Concentration (mM)')

hfig3 = figure(3);
clf
hold on 
plot(tspan,Cr,'b','linewidth',2)
plot(tspan,CrP,'r','linewidth',2)
legend('Cr','CrP')
set(gca,'FontSize',20)
xlabel('Time (s)')
ylabel('Concentration (mM)')
xlim([tspan(1) tspan(end)])
ylim([-.01 Cr_tot*1e3+.01])

if printon == 1
    print(hfig1,'-dpng','adenyls.png')
    print(hfig2,'-dpng','phosphates.png')
    print(hfig3,'-dpng','creatines.png')
end 





