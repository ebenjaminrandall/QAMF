function dxdt = model(t,x,pars,consts)


%% Constants 

V_e = consts(1); 
V_x = consts(2); 

R = consts(3); 
T = consts(4);
F = consts(5);

pH_outer = consts(6); 
pH_inner = consts(7);

Mg = consts(8); 
K  = consts(9);

psi = consts(10); 

Cr_tot = consts(11);

%% Parameters

K_HATP = pars(1); 
K_KATP = pars(2); 
K_MATP = pars(3);
K_HADP = pars(4); 
K_KADP = pars(5); 
K_MADP = pars(6);
K_HPi  = pars(7); 
K_KPi  = pars(8); 
K_MPi  = pars(9);

% J_F1F0
X_F1F0 = pars(10); 
n_H  = pars(11); 
DGr0_F1F0 = pars(12);

% J_ANT 
E_ANT = pars(13); 
k2_ANT_o = pars(14); 
k3_ANT_o = pars(15); 
Ko_Do = pars(16); 
Ko_To = pars(17); 
alpha_1 = pars(18); 
alpha_2 = pars(19); 
alpha_3 = pars(20); 
delta_T = pars(21); 
delta_D = pars(22);

% J_ATPase 
X_ATPase = pars(23); 
K_iADP   = pars(24); 
DGr0_ATPase = pars(25); 

% J_PiC
X_PiC = pars(26); 
k_PiC = pars(27);

% J_CK
DGr0_CK = pars(28); 
X_CK    = pars(29); 
kref    = pars(30);

%% States 

ATP_x = x(1); 
ADP_x = x(2); 
Pi_x  = x(3); 
ATP_e = x(4); 
ADP_e = x(5); 
Pi_e  = x(6); 
Cr    = x(7); 

CrP = Cr_tot - Cr; 

%% pH 

H_x = 10^(-pH_inner); 
H_e = 10^(-pH_outer); 

%% Binding polynominals

P_ATP = 1 + H_x/K_HATP + Mg/K_MATP + K/K_KATP;
P_ADP = 1 + H_x/K_HADP + Mg/K_MADP + K/K_KADP;
P_Pi  = 1 + H_x/K_HPi  + Mg/K_MPi  + K/K_KPi; 

%% J_F1F0

Keq_prime_F1F0 = exp(-(DGr0_F1F0 - F * psi * n_H)/(R * T)) ... 
    * (P_ATP / (P_ADP * P_Pi)) ...
    * (H_e^n_H / H_x^(n_H - 1)); 

J_F1F0 = X_F1F0 * (Keq_prime_F1F0 * Pi_x * ADP_x - ATP_x); 

%% J_ANT 

k2_ANT = k2_ANT_o*exp((alpha_1*(-3)+alpha_2*(-4)+alpha_3)*F*psi/(R * T));
k3_ANT = k3_ANT_o*exp((alpha_1*(-4)+alpha_2*(-3)+alpha_3)*F*psi/(R * T));

Ko_D = Ko_Do * exp(3*delta_D*F*psi/(R * T));
Ko_T = Ko_To * exp(4*delta_T*F*psi/(R * T));
q = k3_ANT * Ko_D * exp(F*psi/(R * T)) / (k2_ANT * Ko_T);

J_ANT = E_ANT * ...
    (k2_ANT * q * ATP_x * ADP_e / Ko_D - k3_ANT * ATP_e * ADP_x / Ko_T ) ... 
    / ((1 + ATP_e/Ko_T + ADP_e/Ko_D) * (ADP_x + q * ATP_x));

%% J_ATPase 

Keq_ATPase = exp(-DGr0_ATPase / (R * T)) / H_e * (P_ADP * P_Pi) / P_ATP; 

J_ATPase = X_ATPase;% * (1 - (ADP_e * Pi_e)/(ATP_e * Keq_ATPase)) ...
    %/ (1 + ADP_e / K_iADP); 

%% J_PiC

J_PiC = X_PiC * (Pi_e * H_e - Pi_x * H_x) / (Pi_e + k_PiC); 

%% J_CK

Keq_CK = kref * H_e * P_ATP / P_ADP;
%Keq_CK = exp(-DGr0_CK / (R * T)) * H_e * P_ATP / P_ADP; 
%{
Changed from equation in Bazil et al. to include kref from Wu et al 2008.
The original equation has a rate that is much too fast for the current
model
%} 

J_CK = X_CK * (Keq_CK * ADP_e * CrP - ATP_e * Cr); 

%% Differential equations 

dATP_x = (J_F1F0 - J_ANT) / V_x; 
dADP_x = (-J_F1F0 + J_ANT) / V_x; 
dPi_x  = (-J_F1F0 + J_PiC) / V_x; 
dATP_e = J_ANT / V_e - J_ATPase + J_CK; 
dADP_e = -J_ANT / V_e + J_ATPase - J_CK; 
dPi_e  = -J_PiC / V_e + J_ATPase; 
dCr    = J_CK; 

dxdt = [dATP_x; dADP_x; dPi_x; 
    dATP_e; dADP_e; dPi_e;
    dCr; 
    ]; 

end 




