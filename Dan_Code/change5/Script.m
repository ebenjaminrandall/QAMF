clear 

%% Flags 

flags   = 1;    % 1 - solves model for dxdt;    2 - solves for auxiliary equations
Hleakon = 0;    % 0 - does not use H+ leak;     1 - uses H+ leak
printon = 0;    % 0 - does not print figues;    1 - prints figures 

%% Parameters defining metabolite pools

% Volume fractions and water space fractions
V_c = 0.6601;       % cytosol volume fraction       % L cyto (L cell)^(-1)
V_m = 0.2882;       % mitochondrial volume fraction % L mito (L cell)^(-1)
R_m2c = V_m / V_c;  % mito to cyto volume ratio     % L mito (L cyto)^(-1)
W_c = 0.8425;       % cytosol water space           % L cyto water (L cyto)^(-1) 
W_m = 0.7238;       % mitochondrial water space     % L mito water (L mito)^(-1)
W_x = 0.9*W_m;      % matrix water space            % L matrix water (L mito)^(-1)
W_i = 0.1*W_m;      % intermembrane water space     % L IM water (L mito)^(-1)

fixedpars.fractions.V_c = V_c; 
fixedpars.fractions.V_m = V_m; 
fixedpars.fractions.R_m2c = R_m2c; 
fixedpars.fractions.W_c = W_c; 
fixedpars.fractions.W_m = W_m; 
fixedpars.fractions.W_x = W_x; 
fixedpars.fractions.W_i = W_i; 

% Total pool concentrations 
NAD_tot = 2.97e-3;  % NAD+ and NADH conc            % mol (L matrix water)^(-1)
Q_tot   = 1.35e-3;  % Q and QH2 conc                % mol (L matrix water)^(-1) 
c_tot   = 2.7e-3;   % cytochrome c ox and red conc  % mol (L IM water)^(-1)

Cr_tot  = 40e-3;    % creatine and creatine phosphate conc % mol (L cell)^(-1)
Cr_tot_c = Cr_tot / (V_c * W_c); % convert to mol (L cyto water)^(-1)

fixedpars.pools.NAD_tot = NAD_tot; 
fixedpars.pools.Q_tot = Q_tot; 
fixedpars.pools.c_tot = c_tot; 
fixedpars.pools.Cr_tot_c = Cr_tot_c; 

% Membrane capacitance 
Cm = 3e-6; 

fixedpars.capacitance.Cm = Cm; 

%% Set fixed pH, cation concentrations, and O2 partial pressure 

% pH 
pH_x = 7.35; 
pH_c = 7.15; 

% K+ concentrations 
K_x  = 100e-3;      % mol (L matrix water)^(-1)
K_c  = 140e-3;      % mol (L cyto water)^(-1)

% Mg2+ concentrations
Mg_x = 1e-3;        % mol (L matrix water)^(-1)
Mg_c = 1e-3;        % mol (L cyto water)^(-1)

% Oxygen partial pressure 
PO2 = 25; % mmHg

conc = [pH_x; pH_c; K_x; K_c; Mg_x; Mg_c; PO2]; 

%% Adjustable parameters 

X_DH   = 0.0866;
X_C1   = 32365;
X_C3   = 0.79081* 40 * 1e5;
X_C4   = 0.00010761 * 7e2;
X_F1F0 = 100; 
E_ANT  = 1.5*0.006762/7.2679e-003*(0.70e-1);
X_PiC  = 3.3356e+07 * 2e-2; 
X_CK   = 1e7;
X_AtC  = 0.50e-3;
X_AK   = 1e8; 

pars = [X_DH; X_C1; X_C3; X_C4; 
    X_F1F0; E_ANT; X_PiC; 
    X_CK; X_AtC; X_AK;
    ]; 

%% Initial conditions 

% Membrane potential
DPsi_0    = 175*1e-3;

% Matrix species 
ATP_x_0   = 0.5e-3;
ADP_x_0   = 9.5e-3;
Pi_x_0    = 0.3e-3;
NADH_x_0  = 2/3 * NAD_tot;
QH2_x_0   = Q_tot/2;

% IM species
cred_i_0  = c_tot/2;

% Cytosol species
ATP_c_0   = 9.95e-3;
ADP_c_0   = 0.05e-3;
Pi_c_0    = 1.0e-3;
AMP_c_0   = 1e-6; 
CrP_c_0   = .3 * Cr_tot_c;

x0 = [DPsi_0; 
    ATP_x_0; ADP_x_0; Pi_x_0; NADH_x_0; QH2_x_0; 
    cred_i_0; 
    ATP_c_0; ADP_c_0; Pi_c_0; AMP_c_0; CrP_c_0; 
    ]; 

%% Solve model 

[t,x] = ode15s(@Cell_dXdT,[0 60],x0,[],pars,conc,fixedpars,flags,Hleakon);

% Find steady-state fluxes and polynomials 
f = zeros(15,length(t)); 
for i = 1:length(t)
    [~,f1] = Cell_dXdT([],x(:,end),pars,conc,fixedpars,2,Hleakon); 
    f(:,i) = f1; 
end 

PATP_c = f(4,:)'; 

% convert to mV 
DPsi = x(:,1) * 1e3; 

% convert to mmol (L matrix water)^(-1)
sumATP_x = x(:,2)  * 1e3;
sumADP_x = x(:,3)  * 1e3; 
sumPi_x  = x(:,4)  * 1e3; 
NADH_x   = x(:,5) * 1e3;

% convert to mmol (L IM water)^(-1)
cred_i = x(:,7) * 1e3; 

% convert to mmol (L cyto water)^(-1)
sumATP_c = x(:,8)  * 1e3; 
sumADP_c = x(:,9)  * 1e3; 
sumPi_c  = x(:,10) * 1e3; 
sumAMP_c = x(:,11) * 1e3; 
CrP_c    = x(:,12) * 1e3; 
Cr_c     = Cr_tot_c * 1e3 - CrP_c;

% convert to mmol (L cell)^(-1)
sumATP_x_cell = sumATP_x * (V_m * W_x); 
sumADP_x_cell = sumADP_x * (V_m * W_x); 
sumPi_x_cell  = sumPi_x  * (V_m * W_x); 

sumATP_c_cell = sumATP_c * (V_c * W_c + V_m * W_i); 
sumADP_c_cell = sumADP_c * (V_c * W_c + V_m * W_i);
sumPi_c_cell  = sumPi_x  * (V_c * W_c + V_m * W_i); 
CrP_c_cell    = CrP_c * (V_c * W_c);
Cr_c_cell     = Cr_c  * (V_c * W_c);  


ATP_c_cell = sumATP_x_cell + sumATP_c_cell; 
CrP_ATP = CrP_c_cell ./ ATP_c_cell; 

%% Plot figures 

figure(1)
plot(t,DPsi)
title('\Delta\Psi (mV)');
set(gca,'FontSize',20)

% figure(2)
% plot(t,CrP_ATP)
% title('[CrP]/[ATP]');
% set(gca,'FontSize',20)

figure(4)
plot(t,NADH_x/(NAD_tot * 1e3))
title('NADH');
set(gca,'FontSize',20)

figure(5)
clf 
plot(t,sumATP_x,'b',t,sumADP_x,'r',t,sumPi_x,'c')
title('Matrix')
legend('ATP','ADP','Pi')
set(gca,'FontSize',20)

figure(6)
clf
plot(t,sumATP_c,'b',t,sumADP_c,'r',t,sumPi_c,'c',t,sumAMP_c,'k')
title('Cytosol')
legend('ATP','ADP','Pi','AMP')
set(gca,'FontSize',20)

figure(7)
clf
plot(t,CrP_c_cell,'b',t,Cr_c_cell,'r')
legend('[CrP]','[Cr]')
title('Creatine')
set(gca,'FontSize',20)

figure(8)
clf
plot(t,cred_i,'b')
title('Cytochrome c')
set(gca,'FontSize',20) 

%% CrP/ATP vs ATP consumption rate 

clear CrP_ATP

% Range of ATP hydr02olysis from 0.36 to 1.2e-3 mmol  / sec / (l cell)
X_AtC = (0:0.1:1.52)*1e-3;
p = pars; 

CrP_ATP = zeros(size(X_AtC)); 
Pi_c = zeros(size(X_AtC)); 
for i = 1:length(X_AtC)
    p(9) = X_AtC(i); 
    [t,x] = ode15s(@Cell_dXdT,[0 60],x0,[],p,conc,fixedpars,1,Hleakon);
    
    sumATP_x_cell = x(end,2) * (V_m * W_x); 
    sumATP_c_cell = x(end,8) * (V_c * W_c + V_m * W_i); 
    CrP_c_cell = x(end,12) * (V_c * W_c); 
    CrP_ATP(i) = CrP_c_cell/(sumATP_c_cell + sumATP_x_cell); 
    
    Pi_c(i) = x(end,10); 
end 

hfig9 = figure(9);
clf
plot(X_AtC,CrP_ATP,'b')
xlabel('ATP consumption rate')
ylabel('[CrP]/[ATP]')
set(gca,'FontSize',20)
ylim([0 3])
xlim([0 15e-4])

hfig10 = figure(10);
clf
plot(X_AtC,Pi_c*1e3,'b')
xlabel('ATP consumption rate')
ylabel('[Pi]_c')
set(gca,'FontSize',20)
%ylim([0 8])
xlim([0 15e-4])

if printon == 1
    savefig(hfig9,'CrP_ATP.fig')
    savefig(hfig10,'Pi_c.fig')
end 