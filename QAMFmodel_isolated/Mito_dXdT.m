function [dxdt,f] = Mito_dXdT(~,x,pars,conc,fixedpars,flags)
% This function is used to calculate time derivatives of state
% variables of the cuvette-level cardaic energetics.

%% Setting adjustable parameter values

X_DH   = pars(1); % mol (s * L mito)^(-1)
X_C1   = pars(2); % mol (s * L mito)^(-1)
X_C3   = pars(3); % mol (s * L mito)^(-1)
X_C4   = pars(4); % mol (s * L mito)^(-1)
X_F1F0 = pars(5); % mol (s * L mito)^(-1)
E_ANT  = pars(6); % (L matrix water) (L mito)^(-1)
X_PiC  = pars(7); 
X_AtC  = pars(8);

%% Import pH, cation concentrations and O2 partial pressure 

pH_x = conc(1); 
pH_c = conc(2); 
K_x  = conc(3); % mol (L matrix water)^(-1)
K_c  = conc(4); % mol (L cuvette water)^(-1)
Mg_x = conc(5); % mol (L matrix water)^(-1)
Mg_c = conc(6); % mol (L cuvette water)^(-1)
PO2  = conc(7); % mmHg

% Hydrogen ion concentration 
H_x = 10^(-pH_x); % mol (L matrix water)^(-1)
H_c = 10^(-pH_c); % mol (L cuvette water)^(-1)

% Oxygen concentration 
a_3  = 1.74e-6;   % oxygen solubility in cuvette   % mol (L matrix water * mmHg)^(-1)
O2_x = a_3*PO2;   % mol (L matrix water)^(-1)

%% Import fixed parameter values 

% Volume fractions and water space fractions 
V_c   = fixedpars.fractions.V_c;        % (L buffer) (L cuvette)^(-1)
V_m   = fixedpars.fractions.V_m;        % (L mito) (L cuvette)^(-1)
R_m2c = fixedpars.fractions.R_m2c;      % (L mito) (L buffer)^(-1)
W_c   = fixedpars.fractions.W_c;        % (L cuvette water) (L buffer)^(-1)
W_x   = fixedpars.fractions.W_x;        % (L matrix water) (L mito)^(-1)
W_i   = fixedpars.fractions.W_i;        % (L IM water) (L mito)^(-1)

% Total pool concentrations
NAD_tot  = fixedpars.pools.NAD_tot;     % mol (L matrix water)^(-1)
Q_tot    = fixedpars.pools.Q_tot;       % mol (L matrix water)^(-1)
c_tot    = fixedpars.pools.c_tot;       % mol (L matrix water)^(-1)

% Membrane capacitance 
Cm = fixedpars.capacitance.Cm; 

%%  Listing fixed model parameters
%  (i) Thermochemical constants
R = 8.314;          % J (mol K)^(-1)
T = 37 + 273.15;    % K 
F = 96484;          % C mol^(-1)

% Proton motive force parameters (dimensionless) 
n_F1F0 = 8/3; 
n_C1   = 4; 
n_C3   = 2; 
n_C4   = 4; 

% Dissociation constants
K_HATP = 10^(-6.59);  
K_MATP = 10^(-3.82);  
K_KATP = 10^(-1.013); 
K_HADP = 10^(-6.42);  
K_MADP = 10^(-2.79);  
K_KADP = 10^(-0.882);
K_HPi  = 10^(-6.71);  
K_MPi  = 10^(-1.69);  
K_KPi  = 10^(+0.0074);

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
cred_i   = x(7);      % mol (L IM water)^(-1)

% Cytoplasmic species
sumATP_c = x(8);    % mol (L cuvette water)^(-1)
sumADP_c = x(9);    % mol (L cuvette water)^(-1)
sumPi_c  = x(10);   % mol (L cuvette water)^(-1)

% Other concentrations computed from the state variables:
NAD_x = NAD_tot - NADH_x;  % mol (L matrix water)^(-1)
Q_x   = Q_tot - QH2_x;     % mol (L matrix water)^(-1)
cox_i = c_tot - cred_i;    % mol (L matrix water)^(-1)

%% Binding polynomials 

% Matrix species % mol (L mito water)^(-1)
PATP_x = 1 + H_x/K_HATP + Mg_x/K_MATP + K_x/K_KATP;
PADP_x = 1 + H_x/K_HADP + Mg_x/K_MADP + K_x/K_KADP;
PPi_x  = 1 + H_x/K_HPi  + Mg_x/K_MPi  + K_x/K_KPi; 

% Cytosol species % mol (L cuvette water)^(-1)
PATP_c = 1 + H_c/K_HATP + Mg_c/K_MATP + K_c/K_KATP;
PADP_c = 1 + H_c/K_HADP + Mg_c/K_MADP + K_c/K_KADP;
PPi_c  = 1 + H_c/K_HPi  + Mg_c/K_MPi  + K_c/K_KPi; 

%% Unbound species 

ADP_x = sumADP_x / PADP_x; 
ATP_x = sumATP_x / PATP_x;
Pi_x  = sumPi_x  / PPi_x; 

ADP_c = sumADP_c / PADP_c; 
ATP_c = sumATP_c / PATP_c; 
Pi_c  = sumPi_c  / PPi_c; 

%% NADH Dehydrogenase

% Constants
r      = 4.559;
k_Pi1  = 0.1553e-3; % mol (L matrix water)^(-1)
k_Pi2  = 0.8222e-3; % mol (L matrix water)^(-1)

% Flux 
J_DH = X_DH * (r * NAD_x - NADH_x) * ...
    ((1 + sumPi_x / k_Pi1) / (1+sumPi_x / k_Pi2));

%% Complex I
% NADH_x + Q_x + 5H+_x <-> NAD+_x + QH2_x + 4H+_i + 4dPsi

% Gibb's energy 
dGro_C1 = -109680;  % J mol^(-1) 

% Equilibrium constants
Keq_C1  = exp(-(dGro_C1 + n_C1 * F * DPsi) / (R * T));
Kapp_C1 = Keq_C1 * H_x^(n_C1 + 1) / H_c^n_C1;

% Flux 
J_C1 = X_C1 * (Kapp_C1 * NADH_x * Q_x - NAD_x * QH2_x);

%% Complex III
% QH2_x + 2cuvetteC(ox)3+_i + 2H+_x <-> Q_x + 2cuvetteC(red)2+_i + 4H+_i + 2dPsi

% Gibb's energy 
dGro_C3 = 46690;    % J mol^(-1) 

% Equilibrium constants
Keq_C3  = exp(-(dGro_C3 + n_C3 * F * DPsi)/(R * T));
Kapp_C3 = Keq_C3 * H_x^n_C3 / H_c^(n_C3 + 2);

% Flux 
J_C3 = X_C3 * (Kapp_C3 * cox_i^2 * QH2_x - cred_i^2 * Q_x);

%% Complex IV
% 2cuvetteC(red)2+_i + 0.5O2_x + 4H+_x <-> 2cuvetteC(ox)3+_x + H2O_x + 2H+_i +
% 2dPsi

% Constants 
k_O2 = 1.2e-4;      % mol (L matrix water)^(-1)

% Gibb's energy 
dGro_C4 = -202160;  % J mol^(-1) 

% Equilibrium constant 
Keq_C4  = exp(-(dGro_C4 + n_C4 * F * DPsi) / (R * T));
Kapp_C4 = Keq_C4 * H_x^n_C4 / H_c^(n_C4 - 2);
    
% Flux
J_C4 = X_C4 *(Kapp_C4 * cred_i^2 * O2_x^0.5 - cox_i^2) * ...
    (1 / (1 + k_O2 / O2_x));

%% F1F0-ATPase
% ADP3-_x + HPO42-_x + H+_x + n_A*H+_i <-> ATP4- + H2O + n_A*H+_x

% Gibb's energy 
dGro_F1 = -4510;    % J mol^(-1) 

% Equilibrium constants
Keq_F1F0  = exp(-(dGro_F1 - n_F1F0 * F * DPsi) / (R * T));
Kapp_F1F0 = Keq_F1F0 * H_c^n_F1F0 / H_x^(n_F1F0-1) * PATP_x / (PADP_x * PPi_x);

% Flux
J_F1F0 = X_F1F0 * (Kapp_F1F0 * sumADP_x * sumPi_x - sumATP_x);

%% ANT
% ATP4-_x + ADP3-_i <-> ATP4-_i + ADP3-_x

del_D   = 0.0167;
del_T   = 0.0699;
k2o_ANT = 9.54/60;      % s^(-1)
k3o_ANT = 30.05/60;     % s^(-1)
K0o_D   = 38.89e-6;     % mol (L cuvette water)^(-1) 
K0o_T   = 56.05e-6;     % mol (L cuvette water)^(-1)
A       = +0.2829;
B       = -0.2086;
C       = +0.2372;

phi = F * DPsi / (R * T); 

% Reaction rates 
k2_ANT = k2o_ANT * exp((A * (-3) + B * (-4) + C) * phi);
k3_ANT = k3o_ANT * exp((A * (-4) + B * (-3) + C) * phi);

% Dissociation constants 
K0_D = K0o_D * exp(3 * del_D * phi);
K0_T = K0o_T * exp(4 * del_T * phi);

q     = k3_ANT * K0_D * exp(phi) / (k2_ANT * K0_T);
term1 = k2_ANT * ATP_x * ADP_c * q / K0_D;
term2 = k3_ANT * ADP_x * ATP_c / K0_T;
num   = term1 - term2;
den   = (1 + ATP_c/K0_T + ADP_c/K0_D) * (ADP_x + ATP_x * q);

% Flux
J_ANT = E_ANT * num / den; 

%% H+-PI2 cotransporter

% Constant
k_PiC = 0.0016125;  % mol (L cuvette)^(-1)

% H2P04- species 
HPi_c = Pi_c * (H_c / K_HPi);
HPi_x = Pi_x * (H_x / K_HPi);

% Flux
J_PiC = X_PiC * (H_c * HPi_c  - H_x * HPi_x  ) / ... 
    (k_PiC * (1 + HPi_c/k_PiC) * (1 + HPi_x/k_PiC) );

%% H LEAK (from Wu et al.)

X_H = 2e3; %1e3;
J_H = X_H * (H_c * exp(phi/2) - H_x * exp(-phi/2));


%% ATPase
% ATP4- + H2O = ADP3- + PI2- + H+

%Flux
J_AtC = X_AtC; 

%% Computing time derivatives of state variables

% Membrane potential 
dDPsi = (n_C1 * J_C1 + n_C3 * J_C3 + n_C4 * J_C4 - n_F1F0 * J_F1F0 - J_ANT - J_H) / Cm;

% Matrix species
dATP_x  = (J_F1F0  - J_ANT) / W_x; 
dADP_x  = (-J_F1F0 + J_ANT) / W_x;
dPi_x   = (-J_F1F0 + J_PiC) / W_x;
dNADH_x = (J_DH  - J_C1)  / W_x; 
dQH2_x  = (J_C1  - J_C3)  / W_x; 

% IM space species
dcred_i = 2 * (J_C3 - J_C4) / W_i; 

% Buffer species
dATP_c = ( R_m2c * J_ANT - J_AtC ) / W_c;  
dADP_c = (-R_m2c * J_ANT + J_AtC ) / W_c; 
dPi_c  = (-R_m2c * J_PiC + J_AtC) / W_c; 


if flags == 1
    dxdt = [dDPsi; 
        dATP_x; dADP_x; dPi_x; dNADH_x; dQH2_x; 
        dcred_i; 
        dATP_c; dADP_c; dPi_c]; 
else
    dxdt = []; 
    f = [PATP_x; PADP_x; PPi_x;
        PATP_c; PADP_c; PPi_c; 
        J_DH; J_C1; J_C3; J_C4; 
        J_F1F0; J_ANT; J_PiC; 
        ]; 
end 


