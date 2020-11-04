function dxdt = model(~,x,X_ATPase)


%% Constants 

% Water volumes 
W_x = 0.650;  % L maxtrix water / L mito 
W_i = W_x / .9 * .1; % L intermembrane space water / L mito
W_e = 0.8425; % L cyto water / L cyto 

% Compartment volumes 
Vmito = 0.2882; % L / L tissue
Vcyto = 0.6801; % L / L tissue

% Membrane capacitance 
C_i = 3.11e-3; %3.11e-6; % M V^(-1)

% Constants 
R      = 8.3145; % J K^(-1) mol^(-1)
T      = 310.15; % K 
F      = 96485.33; % C mol^(-1)
pH_e   = 7.2; 
pH_x   = 7.4; 
Mg     = 1e-3; % M 
K      = 100e-3;  % M 
O2     = 2.6e-5; % M

% Stoichiometric coefficients for proton motive forces
n_H  = 8/3; % F1F0_ATPase transfer 
n_C1 = 4;   % Complex I 
n_C3 = 2;   % Complex III
n_C4 = 4;   % Complex IV

% Total concentration pools (M)
Cr_tot   = 40e-3;   % total creatine 
Q_tot    = 1.35e-3; % total ubiquinol  
NAD_tot  = 2.97e-3; % total NAD 
cytc_tot = 2.7e-3;   %total cytochrome c

% Dissociation constants
K_HATP = 2.757e-7;
K_KATP = 9.809e-2;
K_MATP = 8.430e-5;
K_HADP = 4.106e-7;
K_KADP = 1.319e-1;
K_MADP = 7.149e-4;
K_HPi  = 2.308e-7;
K_KPi  = 3.803e-1;
K_MPi  = 2.815e-2;

%% Parameters

X_F1F0   = 812; %8120.4 in Dan's code 
X_ANT    = 0.141;           %mol (L mito)^(-1)
X_ATPase = ; 
X_PiH    = 3.34e7; %mol/s/M/lmito
X_CK     = 2e9; %M/s % 1e7 in Dan's code 
X_C1     = 4405; % 32365 in Dan's code 
X_C3     = 4.887; 
X_C4     = 6.766e-5; 
X_DH     = 0.1099; 

%% States 

ATP_tot_x = x(1); 
ADP_tot_x = x(2); 
Pi_tot_x  = x(3); 
ATP_tot_e = x(4); 
ADP_tot_e = x(5); 
Pi_tot_e  = x(6); 
Cr        = x(7); 
QH2       = x(8); 
NADH_x    = x(9); 
cytc_red  = x(10); 
DPsi      = x(11); 

CrP     = Cr_tot - Cr; 
Q       = Q_tot - QH2; 
NAD_x   = NAD_tot - NADH_x; 
cytc_ox = cytc_tot - cytc_red; 

%% Hydrogen ion concentration 

H_x = 10^(-pH_x);
H_e = 10^(-pH_e);

%% Binding polynominals

P_ATP_x = 1 + H_x/K_HATP + Mg/K_MATP + K/K_KATP;
P_ADP_x = 1 + H_x/K_HADP + Mg/K_MADP + K/K_KADP;
P_Pi_x  = 1 + H_x/K_HPi  + Mg/K_MPi  + K/K_KPi; 

P_ATP_e = 1 + H_e/K_HATP + Mg/K_MATP + K/K_KATP;
P_ADP_e = 1 + H_e/K_HADP + Mg/K_MADP + K/K_KADP;
P_Pi_e  = 1 + H_e/K_HPi  + Mg/K_MPi  + K/K_KPi;

%% J_F1F0

DGro_F1F0 = -4.51 * 1e3; % J mol^(-1)

Keq_F1F0  = exp(-(DGro_F1F0 - F * DPsi * n_H)/(R * T)); 
Kapp_F1F0 =  Keq_F1F0 * (P_ATP_x / (P_ADP_x * P_Pi_x)) ...
    * (H_e^n_H / H_x^(n_H - 1)); 

J_F1F0 = X_F1F0 * (Kapp_F1F0 * Pi_tot_x * ADP_tot_x - ATP_tot_x); 

%% J_ANT 

del_D   = 0.0167;
del_T   = 0.0699;
k2o_ANT = 9.54 / 60;
k3o_ANT = 30.06 / 60;
K0o_D   = 38.89 * 1e-6;
K0o_T   = 56.05 * 1e-6;
A       = +0.2829;
B       = -0.2086;
C       = +0.2372;

phi = F * DPsi / (R * T);

ADP_e = ADP_tot_e / P_ADP_e; % [ADP^3-]_e;
ATP_e = ATP_tot_e / P_ATP_e; % [ATP^4-]_e;
ADP_x = ADP_tot_x / P_ADP_x; % [ADP^3-]_x;
ATP_x = ATP_tot_x / P_ATP_x; % [ATP^4-]_x;

k2_ANT = k2o_ANT * exp((A * (-3) + B * (-4) + C) * phi);
k3_ANT = k3o_ANT * exp((A * (-4) + B * (-3) + C) * phi);

K0_D = K0o_D * exp(3 * del_D * phi);
K0_T = K0o_T * exp(4 * del_T * phi);

q     = k3_ANT * K0_D * exp(phi) / (k2_ANT * K0_T);
term1 = k2_ANT * ATP_x * ADP_e * q / K0_D;
term2 = k3_ANT * ADP_x * ATP_e / K0_T;
num   = term1 - term2;
den   = (1 + ATP_e/K0_T + ADP_e/K0_D) * (ADP_x + ATP_x * q);

J_ANT = X_ANT*num/den;

%% J_ATPase 

J_ATPase = X_ATPase;

%% J_PiH

k_PiH = 1.61e-3;

H2PO4_e = (H_e / K_HPi) * (Pi_tot_e / P_Pi_e);
H2PO4_x = (H_x / K_HPi) * (Pi_tot_x / P_Pi_x);

J_PiH = X_PiH * (H_e * H2PO4_e - H_x * H2PO4_x) / ... 
    (k_PiH * (1 + H2PO4_e / k_PiH) * (1 + H2PO4_x / k_PiH));

%% J_CK

Keq_CK  = 3.57e8; %Wu et al 2007 %%5.18e8; %7.14e8; % Wu et al 2008
Kapp_CK = Keq_CK * H_e * P_ATP_e / P_ADP_e;
%{
Changed from equation in Bazil et al. to include kref from Wu et al 2008.
The original equation has a rate that is much too fast for the current
model
%} 

J_CK = X_CK * (Kapp_CK * ADP_tot_e * CrP - ATP_tot_e * Cr); 

%% J_C1 

DGro_C1 = -69.37 * 1000; %J mol^(-1) %% -109.68 in Dan's code 

Keq_C1  = exp(-(DGro_C1 + n_C1 * F * dPsi) / (R * T));
Kapp_C1 = Keq_C1 * H_x^(n_C1 + 1) / H_e^n_C1;

J_C1 = X_C1 * (Kapp_C1 * NADH_x * Q_x - NAD_x * QH2_x );

%% J_C3 

DGro_C3 = -32.53 * 1000; %J mol^(-1) %% +46.69 in Dan's Code 

Keq_C3  = exp(-(DGro_C3 + n_C3 * F * DPsi) / (R * T)); 
Kapp_C3 = Keq_C3 * H_x^n_C3 / H_e^(n_C3 + 2); 

J_C3 = X_C3 * (Kapp_C3 * cytc_ox^2 * QH2 - cytc_red^2 * Q); 

%% J_C4 

DGro_C4 = -122.94 * 1000; %J mol^(-1) %% -202.16 in Dan's code 
k_O2    = 120e-6;

Keq_C4  = exp(-(DGro_C4 + n_C4 * F * DPsi) / (R * T)); 
Kapp_C4 = Keq_C4 * H_x^n_C4 / H_e^(n_C4 - 2); 

f_C4   = (1 / (1 + k_O2 / O2)); 

J_C4 = X_C4 * f_C4 * ... 
    (Kapp_C4 * cytc_red^2 * O2^(1/2) - cytc_ox^2);

%% J_DH 

r = 4.253; % What is this constant? 

J_DH = X_DH * (r * NAD_x - NADH_x); 

%% Differential equations 

dATP_x   = (+J_F1F0 - J_ANT) / W_x; 
dADP_x   = (-J_F1F0 + J_ANT) / W_x; 
dPi_x    = (-J_F1F0 + J_PiH) / W_x; 
dATP_e   = (+J_ANT*Vmito/Vcyto - J_ATPase + J_CK) / W_e; 
dADP_e   = (-J_ANT*Vmito/Vcyto + J_ATPase - J_CK) / W_e; 
dPi_e    = (-J_PiH*Vmito/Vcyto + J_ATPase) / W_e; 
dCr      = J_CK / W_e;
dQH2     = (+J_C1 - J_C3) / W_x;
dNADH_x  = (J_DH - J_C1) / W_x; 
dcytCred = 2 * (J_C3 - J_C4) / W_i; 
dDPsi    = (n_C1 * J_C1 + n_C3 * J_C3 + n_C4 * J_C4 - n_H * J_F1F0 - J_ANT) / C_i; 

dxdt = [dATP_x; dADP_x; dPi_x; 
    dATP_e; dADP_e; dPi_e;
    dCr; dQH2; dNADH_x; dcytCred; 
    dDPsi
    ]; 

end 




