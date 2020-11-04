function dxdt = model(~,x,X_ATPase)


%% Constants 

V_x = 0.650;  % L maxtrix water / L mito 
V_i = V_x / .9 * .1; % L intermembrane space water / L mito
V_e = 0.8425; % L cyto water / L cyto 
Vmito = 0.2882; % L / L tissue
Vcyto = 0.6801; % L / L tissue

R      = 8.3145;
T      = 310.15;
F      = 96485.33; %C/mole
pH_e   = 7.2; 
pH_x   = 7.4; 
Mg     = 1e-3; %M 
K      = 100e-3;  %M 
DPsi   = 0.175; % V 
O2     = 2.6e-5; 

% Total concentration pools
Cr_tot   = 54e-3; % total creatine (M)
Q_tot    = 1.35e-3; 
NAD_tot  = 2.97e-3; 
cytC_tot = 2.7e-3; 

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

%% States 

ATP_x   = x(1); 
ADP_x   = x(2); 
Pi_x    = x(3); 
ATP_e   = x(4); 
ADP_e   = x(5); 
Pi_e    = x(6); 
Cr      = x(7); 
QH2     = x(8); 
NADH_x  = x(9); 
cytCred = x(10); 

CrP    = Cr_tot - Cr; 
Q      = Q_tot - QH2; 
NAD_x  = NAD_tot - NADH_x; 
cytCox = cytC_tot - cytCred; 

%% Hydrogen ion concentration 

H_x = 10^(-pH_x);
H_e = 10^(-pH_e);
H_w = 10^(-7); 

%% Binding polynominals

P_ATPx = 1 + H_x/K_HATP + Mg/K_MATP + K/K_KATP;
P_ADPx = 1 + H_x/K_HADP + Mg/K_MADP + K/K_KADP;
P_PIx  = 1 + H_x/K_HPi  + Mg/K_MPi  + K/K_KPi; 

P_ATPe = 1 + H_e/K_HATP + Mg/K_MATP + K/K_KATP;
P_ADPe = 1 + H_e/K_HADP + Mg/K_MADP + K/K_KADP;
P_PIe  = 1 + H_e/K_HPi  + Mg/K_MPi  + K/K_KPi; 

%% J_F1F0

X_F1F0    = 812; 
n_H       = 8/3;
DGr0_F1F0 = -4510;

Keq_prime_F1F0 = exp(-(DGr0_F1F0 - F * DPsi * n_H)/(R * T)) ... 
    * (P_ATPx / (P_ADPx * P_PIx)) ...
    * (H_e^n_H / H_x^(n_H - 1)); 

J_F1F0 = X_F1F0 * (Keq_prime_F1F0 * Pi_x * ADP_x - ATP_x); 

%% J_ANT 

X_ANT     = 0.141;           %mol (L mito)^(-1)
del_D     = 0.0167;
del_T     = 0.0699;
k2_ANT    = 9.54/60;
k3_ANT    = 30.05/60;
K_D_o_ANT = 38.89e-6;
K_T_o_ANT = 56.05e-6;
A         = +0.2829;
B         = -0.2086;
C         = +0.2372;

fi = F*DPsi/(R * T);

ADP_i1 = ADP_e/P_ADPe; % [ADP^3-]_e;
ATP_i1 = ATP_e/P_ATPe; % [ATP^4-]_e;
ADP_x1 = ADP_x/P_ADPx; % [ADP^3-]_x;
ATP_x1 = ATP_x/P_ATPx; % [ATP^4-]_x;

k2_ANT_fi   = k2_ANT*exp((A*(-3)+B*(-4)+C)*fi);
k3_ANT_fi   = k3_ANT*exp((A*(-4)+B*(-3)+C)*fi);

K_D_o_ANT_fi = K_D_o_ANT*exp(3*del_D*fi);
K_T_o_ANT_fi = K_T_o_ANT*exp(4*del_T*fi);

q     = k3_ANT_fi*K_D_o_ANT_fi*exp(fi)/(k2_ANT_fi*K_T_o_ANT_fi);
term2 = k2_ANT_fi*ATP_x1*ADP_i1*q/K_D_o_ANT_fi ;
term3 = k3_ANT_fi.*ADP_x1*ATP_i1/K_T_o_ANT_fi;
num   = term2 - term3;
den   = (1 + ATP_i1/K_T_o_ANT_fi + ADP_i1/K_D_o_ANT_fi)*(ADP_x1 + ATP_x1*q);

J_ANT = X_ANT*num/den;

%% J_ATPase 

J_ATPase = X_ATPase;

%% J_PiC

x_PIH = 3.34e7; %mol/s/M/lmito
k_PIH = 1.61e-3;

a = Pi_e*(H_e/K_HPi)/P_PIe;
p = Pi_x*(H_x/K_HPi)/P_PIx;

J_PiC = (x_PIH/k_PIH)*(H_e*a - H_x*p)/(1+a/k_PIH)/(1+k_PIH);

%% J_CK

X_CK = 1e6; %M/s
kref = 7.14e8; % Wu et al 2008=
Keq_CK = kref * H_e * P_ATPe / P_ADPe;
%Keq_CK = exp(-DGr0_CK / (R * T)) * H_e * P_ATP / P_ADP; 
%{
Changed from equation in Bazil et al. to include kref from Wu et al 2008.
The original equation has a rate that is much too fast for the current
model
%} 

J_CK = X_CK * (Keq_CK * ADP_e * CrP - ATP_e * Cr); 

%% J_C1 

DGo_C1 = -69.37 * 1000; %J mol^(-1)
X_C1 = 4405; 

DG_H = F * DPsi + R * T * log(H_e / H_x); 
Keq_C1 = exp(-(DGo_C1 + 4 * DG_H - R * T * log(H_x / H_w)) / (R * T)); 

J_C1 = X_C1 * (Keq_C1 * NADH_x * Q - NAD_x * QH2); 

%% J_C3 

DGo_C3 = -32.53 * 1000; %J mol^(-1)
X_C3 = 4.887; 

Keq_C3 = exp(-(DGo_C3 + 4 * DG_H + 2 * R * T * log(H_x / H_w) ...
    - 2 * F * DPsi) / (2 * R * T)); 

J_C3 = X_C3 * (Keq_C3 * cytCox * sqrt(QH2) - cytCred * sqrt(Q)); 

%% J_C4 

DGo_C4 = -122.94 * 1000; %J mol^(-1)
X_C4 = 6.766e-5; 
k_O2   = 120e-6;

Keq_C4 = exp(-(DGo_C4 + 2 * DG_H - 2 * R * T * log(H_x / H_w)) / (2 * R * T)); 
f_C4   = (1 / (1 + k_O2 / O2)) * (cytCred / cytC_tot); 

J_C4 = X_C4 * f_C4 * ... 
    (Keq_C4 * cytCred * O2^(1/4) - cytCox * exp(F * DPsi / (R * T)));

%% Differential equations 

dATP_x   = (+J_F1F0 - J_ANT) / V_x; 
dADP_x   = (-J_F1F0 + J_ANT) / V_x; 
dPi_x    = (-J_F1F0 + J_PiC) / V_x; 
dATP_e   = (+J_ANT*Vmito/Vcyto - J_ATPase + J_CK) / V_e; 
dADP_e   = (-J_ANT*Vmito/Vcyto + J_ATPase - J_CK) / V_e; 
dPi_e    = (-J_PiC*Vmito/Vcyto + J_ATPase) / V_e; 
dCr      = J_CK / V_e;
dQH2     = (+J_C1 - J_C3) / V_x;
dNADH_x  = -J_C1 / V_x; 
dcytCred = 2 * (J_C3 - J_C4) / V_i; 



dxdt = [dATP_x; dADP_x; dPi_x; 
    dATP_e; dADP_e; dPi_e;
    dCr; dQH2; dNADH_x; dcytCred; 
    ]; 

end 




