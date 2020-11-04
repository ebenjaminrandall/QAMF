
clear all
close all

% Constants
R = 8.314; % J / mol / K
T = 273.15 + 37; % K

% Dissociation constants (M)
K_HATP = 2.757e-7; 
K_HADP = 4.106e-7; 
K_HPi  = 2.308e-7; 

K_KATP = 9.809e-2; 
K_KADP = 1.319e-1; 
K_KPi  = 3.803e-1; 

K_MgATP = 8.43e-5; 
K_MgADP = 7.149e-4; 
K_MgPi  = 2.815e-2; 

% Equilibrium constant
K_eq = 6; 

% Concentrations (M)
Mg     = 1e-3;
H      = 10^(-7.2); 
SigATP = 0.5e-3; 
SigADP = 9.5e-3; 
SigPi  = 1e-3; 
K      = 150e-3; 

% Polynomials
P_ATP = 1 + Mg/K_MgATP + H/K_HATP + K/K_KATP; 
P_ADP = 1 + Mg/K_MgADP + H/K_HADP + K/K_KADP; 
P_Pi = 1 + Mg/K_MgPi + H/K_HPi + K/K_KPi; 


K_apparent = K_eq * H * P_ATP / (P_ADP * P_Pi); 
DeltaG0_apparent = -R * T * log(K_apparent); 
Q_r = SigATP / (SigADP * SigPi); 
DeltaG = DeltaG0_apparent + R * T * log(Q_r)

return 





% Method 1
MgATP = (Mg / K_MgATP) * (SigATP / P_ATP); 
MgADP = (Mg / K_MgADP) * (SigADP / P_ADP); 
Pi    = SigPi / P_Pi; 
DeltaG0 = -R * T * log(K_eq);
DeltaG  = DeltaG0 - R * T * log(MgADP * Pi * H / MgATP)

% Method 2
K_apparent = K_eq * (P_ATP / (P_ADP * P_Pi)) * H * (K_MgATP / K_MgADP)
DeltaG0app = -R * T * log(K_apparent) 
DeltaG  = DeltaG0app + R * T * log(SigATP / (SigADP * SigPi))







