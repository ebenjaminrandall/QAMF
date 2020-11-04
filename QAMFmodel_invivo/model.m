function [dxdt,f] = model(~,x,pars,conc,fixedpars,flags,Hleakon)
%{
This function calculates time derivatives of state variables of the 
cell-level cardaic energetics. 
Inputs: 
    x       vector of state variables
    pars    vector of adjustable parameters 
    conc    vector of pH, K+ and Mg2+ concentration, and O2 partial
            pressure
    flags   scalar that if set to 1, solves dxdt and otherwise solves
            auxiliary equations 
    Hleakon scalar that if set to 1, uses the H+ leak and otherwise doesn't
Outputs: 
    dxdt    vector of the time derivatives of model states 
    f       vector of binding polynomials and fluxes calculated 
%} 

%% Unpack adjustable parameter vector

X_DH  = pars(1);  % mol (s * L mito)^(-1)
X_C1  = pars(2);  % mol (s * L mito)^(-1)
X_C3  = pars(3);  % mol (s * L mito)^(-1)
X_C4  = pars(4);  % mol (s * L mito)^(-1)
X_F   = pars(5);  % mol (s * L mito)^(-1)
E_ANT = pars(6);  % mol (L mito)^(-1)
E_PiC = pars(7);  % (L cell) (s * L mito)^(-1)
X_CK  = pars(8);  % mol (s * L cyto)^(-1)
X_AtC = pars(9);  % mol (s * L cyto)^(-1)
X_AK  = pars(10); % mol (s * L cyto)^(-1)

%% Upack pH, cation concentrations and O2 partial pressure vector

pH_x = conc(1); 
pH_c = conc(2); 
K_x  = conc(3); % mol (L matrix water)^(-1)
K_c  = conc(4); % mol (L cyto water)^(-1)
Mg_x = conc(5); % mol (L matrix water)^(-1)
Mg_c = conc(6); % mol (L cyto water)^(-1)
PO2  = conc(7); % mmHg

% [H+] 
H_x = 10^(-pH_x); % mol (L matrix water)^(-1)
H_c = 10^(-pH_c); % mol (L cyto water)^(-1)

% [O2]_x 
a_3  = 1.74e-6;   % oxygen solubility in cell   % mol (L matrix water * mmHg)^(-1)
O2_x = a_3 * PO2;   % mol (L matrix water)^(-1)

%% Import fixed parameters and specifications 

% Volume fractions and water space fractions 
V_c   = fixedpars.fractions.V_c;        % (L cyto) (L cell)^(-1)
V_m   = fixedpars.fractions.V_m;        % (L mito) (L cell)^(-1)
V_m2c = fixedpars.fractions.V_m2c;      % (L mito) (L cyto)^(-1)
W_c   = fixedpars.fractions.W_c;        % (L cyto water) (L cyto)^(-1)
W_x   = fixedpars.fractions.W_x;        % (L matrix water) (L mito)^(-1)
W_i   = fixedpars.fractions.W_i;        % (L IM water) (L mito)^(-1)

% Total pool concentrations
NAD_tot  = fixedpars.pools.NAD_tot;     % mol (L matrix water)^(-1)
Q_tot    = fixedpars.pools.Q_tot;       % mol (L matrix water)^(-1)
c_tot    = fixedpars.pools.c_tot;       % mol (L matrix water)^(-1)
Cr_tot_c = fixedpars.pools.Cr_tot_c;    % mol (L cyto water)^(-1) 

% Membrane capacitance 
Cm = fixedpars.capacitance.Cm;          % mol (V * L mito)^(-1)

%%  Listing fixed model parameters

% Thermochemical constants
R = 8.314;          % J (mol K)^(-1)
T = 37 + 273.15;    % K 
F = 96485;          % C mol^(-1)

% Proton motive force parameters (dimensionless) 
n_F  = 8/3; 
n_C1 = 4; 
n_C3 = 2; 
n_C4 = 4; 

% Dissociation constants (M)
K_MgATP = 10^(-3.88);  
K_HATP  = 10^(-6.33);  
K_KATP  = 10^(-1.02); 
K_MgADP = 10^(-3.00);  
K_HADP  = 10^(-6.26);  
K_KADP  = 10^(-0.89);
K_MgAMP = 10^(-1.86); 
K_HAMP  = 10^(-6.22); 
K_KAMP  = 10^(-1.05); 
K_MgPi  = 10^(-1.66);  
K_HPi   = 10^(-6.62);  
K_KPi   = 10^(-0.42);

%% Loading values of state variables

% Membrane potential 
DPsi = x(1);        %V 

% Matrix species 
sumATP_x = x(2);    % mol (L matrix water)^(-1)
sumADP_x = x(3);    % mol (L matrix water)^(-1)
sumPi_x  = x(4);    % mol (L matrix water)^(-1)
NADH_x   = x(5);    % mol (L matrix water)^(-1)
QH2_x    = x(6);    % mol (L matrix water)^(-1)

% IM space species
cred_i = x(7);      % mol (L IM water)^(-1)

% Cytoplasmic species
sumATP_c = x(8);    % mol (L cyto water)^(-1)
sumADP_c = x(9);    % mol (L cyto water)^(-1)
sumPi_c  = x(10);   % mol (L cyto water)^(-1)
sumAMP_c = x(11);   % mol (L cyto water)^(-1)
CrP_c    = x(12);   % mol (L cyto water)^(-1)

% Other concentrations computed from the state variables:
NAD_x = NAD_tot - NADH_x;  % mol (L matrix water)^(-1)
Q_x   = Q_tot - QH2_x;     % mol (L matrix water)^(-1)
cox_i = c_tot - cred_i;    % mol (L matrix water)^(-1)
Cr_c  = Cr_tot_c - CrP_c;  % mol (L cyto water)^(-1)

%% Binding polynomials 

% Matrix species % mol (L mito water)^(-1)
PATP_x = 1 + H_x/K_HATP + Mg_x/K_MgATP + K_x/K_KATP;
PADP_x = 1 + H_x/K_HADP + Mg_x/K_MgADP + K_x/K_KADP;
PPi_x  = 1 + H_x/K_HPi  + Mg_x/K_MgPi  + K_x/K_KPi; 

% Cytosol species % mol (L cyto water)^(-1)
PATP_c = 1 + H_c/K_HATP + Mg_c/K_MgATP + K_c/K_KATP;
PADP_c = 1 + H_c/K_HADP + Mg_c/K_MgADP + K_c/K_KADP;
PPi_c  = 1 + H_c/K_HPi  + Mg_c/K_MgPi  + K_c/K_KPi;
PAMP_c = 1 + H_c/K_HAMP + Mg_c/K_MgAMP + K_c/K_KAMP; 

%% Unbound species 

ATP_x = sumATP_x / PATP_x; % [ATP4-]_x
ADP_x = sumADP_x / PADP_x; % [ADP3-]_x
Pi_x  = sumPi_x  / PPi_x;  % [HPO42-]_x

ATP_c = sumATP_c / PATP_c; % [ATP4-]_c
ADP_c = sumADP_c / PADP_c; % [ADP3-]_c
AMP_c = sumAMP_c / PAMP_c; % [AMP2-]_c
Pi_c  = sumPi_c  / PPi_c;  % [HPO42-]_c

%% NADH Dehydrogenase

% Constants
r      = 4.559;
k_Pi1  = 0.1553e-3; % mol (L matrix water)^(-1)
k_Pi2  = 0.8222e-3; % mol (L matrix water)^(-1)

% Flux (mol (s * L mito)^(-1))
J_DH = X_DH * (r * NAD_x - NADH_x) * ...
    ((1 + sumPi_x / k_Pi1) / (1+sumPi_x / k_Pi2));

%% Complex I
% NADH_x + Q_x + 5 H+_x <-> NAD+_x + QH2_x + 4 H+_c + 4 DPsi

% Gibb's energy (J mol^(-1))
DrGo_C1 = -109680; 

% Equilibrium constants (dimensionless)
Keq_C1  = exp(-(DrGo_C1 + n_C1 * F * DPsi) / (R * T));
Kapp_C1 = Keq_C1 * H_x^(n_C1 + 1) / H_c^n_C1;

% Flux (mol (s * L mito)^(-1))
J_C1 = X_C1 * (Kapp_C1 * NADH_x * Q_x - NAD_x * QH2_x);

%% Complex III
% QH2_x + 2 c3+_i + 2 H+_x <-> Q_x + 2 c2+_i + 4 H+_i + 2 DPsi

% Gibb's energy (J mol^(-1))
DrGo_C3 = 46690;     

% Equilibrium constants dimensionless
Keq_C3  = exp(-(DrGo_C3 + n_C3 * F * DPsi)/(R * T));
Kapp_C3 = Keq_C3 * H_x^n_C3 / H_c^(n_C3 + 2);

% Flux (mol (s * L mito)^(-1))
J_C3 = X_C3 * (Kapp_C3 * cox_i^2 * QH2_x - cred_i^2 * Q_x);

%% Complex IV
% 2 c2+_i + 0.5 O2_x + 4 H+_x <-> 2 c3+_x + H2O_x + 2 H+_i + 2 DPsi

% Constants 
k_O2 = 1.2e-4;      % mol (L matrix water)^(-1)

% Gibb's energy (J mol^(-1))
DrGo_C4 = -202160;   

% Equilibrium constant (dimensionless) 
Keq_C4  = exp(-(DrGo_C4 + n_C4 * F * DPsi) / (R * T));
Kapp_C4 = Keq_C4 * H_x^n_C4 / H_c^(n_C4 - 2);
    
% Flux (mol (s * L mito)^(-1))
J_C4 = X_C4 *(Kapp_C4 * cred_i^2 * O2_x^0.5 - cox_i^2) * ...
    (1 / (1 + k_O2 / O2_x));

%% F0F1-ATPase
% ADP3-_x + HPO42-_x + H+_x + 8/3 H+_i <-> ATP4- + H2O + 8/3 H+_x

% Gibb's energy (J mol^(-1))
DrGo_F = -4510; 

% Equilibrium constants (dimensionless)
Keq_F  = exp(-(DrGo_F - n_F * F * DPsi) / (R * T));
Kapp_F = Keq_F * H_c^n_F / H_x^(n_F-1) * PATP_x / (PADP_x * PPi_x) * K_MgATP / K_MgADP;

% Flux (mol (s * L mito)^(-1))
J_F = X_F * (Kapp_F * sumADP_x * sumPi_x - sumATP_x);

%% ANT
% ATP4-_x + ADP3-_i <-> ATP4-_i + ADP3-_x

% Constants
del_D   = 0.0167;
del_T   = 0.0699;
k2o_ANT = 9.54/60;      % s^(-1)
k3o_ANT = 30.05/60;     % s^(-1)
Ko_D0   = 38.89e-6;     % mol (L cyto water)^(-1) 
Ko_T0   = 56.05e-6;     % mol (L cyto water)^(-1)
A       = +0.2829;
B       = -0.2086;
C       = +0.2372;

phi = F * DPsi / (R * T); 

% Reaction rates (s^(-1)) 
k2_ANT = k2o_ANT * exp((A * (-3) + B * (-4) + C) * phi);
k3_ANT = k3o_ANT * exp((A * (-4) + B * (-3) + C) * phi);

% Dissociation constants (mol (L cyto water)^(-1))
Ko_D = Ko_D0 * exp(3 * del_D * phi);
Ko_T = Ko_T0 * exp(4 * del_T * phi);

q     = k3_ANT * Ko_D * exp(phi) / (k2_ANT * Ko_T);
term1 = k2_ANT * ATP_x * ADP_c * q / Ko_D;
term2 = k3_ANT * ADP_x * ATP_c / Ko_T;
num   = term1 - term2;
den   = (1 + ATP_c/Ko_T + ADP_c/Ko_D) * (ADP_x + ATP_x * q);

% Flux (mol (s * L mito)^(-1))
J_ANT = E_ANT * num / den; 

%% H+-PI2 cotransporter
% H2PO42-_x + H+_x = H2PO42-_c + H+_c

% Constant
k_PiC = 0.0016125;  % mol (L cell)^(-1)

% [H2P04-] 
HPi_c = Pi_c * (H_c / K_HPi);
HPi_x = Pi_x * (H_x / K_HPi);

% convert to mol (L cell)^(-1)
H_c_cell = H_c * (V_c * W_c + V_m * W_i); 
H_x_cell = H_x * V_m * W_x; 
HPi_c_cell = HPi_c * (V_c * W_c + V_m * W_i); 
HPi_x_cell = HPi_x * V_m * W_x; 

% Flux (mol (s * L mito)^(-1))
J_PiC = E_PiC * (H_c_cell * HPi_c_cell - H_x_cell * HPi_x_cell ) / ... 
    (k_PiC * (1 + HPi_c_cell/k_PiC) * (1 + HPi_x_cell/k_PiC));

%% Creatine kinase reaction
% ADP3- + CrP2- + H+ = ATP4- + Cr

% Equilibrium constants (dimensionless)
Keq_CK  = 3.5e8; 
Kapp_CK = Keq_CK * H_c * PATP_c / PADP_c;

% Flux (mol (s * L cyto)^(-1))
J_CK = X_CK * (Kapp_CK * ADP_c * CrP_c - ATP_c * Cr_c);

%% ATPase
% ATP4- + H2O = ADP3- + HPO42- + H+

%Flux (mol (s * L cyto)^(-1))
J_AtC = X_AtC; 

%% Adenylate kinase 
% 2 ADP3- = ATP4- + AMP2- 

% Equilibrium constants (dimensionless)
Keq_AK  = 3.97e-1; 
Kapp_AK = Keq_AK * PATP_c * PAMP_c / PADP_c^2; 

% Flux (mol (s * L cyto)^(-1))
J_AK = X_AK * (Kapp_AK * ADP_c^2 - AMP_c * ATP_c); 

%% H+ Leak

if Hleakon == 1
    % Constants
    X_H = 1e2;  % mol (s * L mito)^(-1))
    
    % Flux (mol (s * L mito)^(-1))
    J_H = X_H * DPsi * (H_c * exp(phi) - H_x) / (exp(phi) - 1); 
else
    J_H = 0;
end 

%% Computing time derivatives of state variables

% Membrane potential 
dDPsi = (n_C1 * J_C1 + n_C3 * J_C3 + n_C4 * J_C4 - n_F * J_F - J_ANT - J_H) / Cm;

% Matrix species
dATP_x  = (J_F  - J_ANT) / W_x; 
dADP_x  = (-J_F + J_ANT) / W_x;
dPi_x   = (-J_F + J_PiC) / W_x;
dNADH_x = (J_DH  - J_C1) / W_x; 
dQH2_x  = (J_C1  - J_C3) / W_x; 

% IM space species
dcred_i = 2 * (J_C3 - J_C4) / W_i; 

% Buffer species
dATP_c = ( V_m2c * J_ANT - J_AtC + J_CK + J_AK) / W_c;  
dADP_c = (-V_m2c * J_ANT + J_AtC - J_CK - 2*J_AK) / W_c; 
dPi_c  = (-V_m2c * J_PiC + J_AtC) / W_c; 
dAMP_c = J_AK / W_c; 
dCrP_c = -J_CK / W_c;

if flags == 1
    dxdt = [dDPsi; 
        dATP_x; dADP_x; dPi_x; dNADH_x; dQH2_x; 
        dcred_i; 
        dATP_c; dADP_c; dPi_c; dAMP_c; dCrP_c]; 
else
    dxdt = []; 
    f = [PATP_x; PADP_x; PPi_x;
        PATP_c; PADP_c; PPi_c; 
        J_DH; J_C1; J_C3; J_C4; 
        J_F; J_ANT; J_PiC; J_CK; J_AK; 
        ]; 
end 


