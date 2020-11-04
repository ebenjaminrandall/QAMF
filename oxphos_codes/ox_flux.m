function [J_C4] = ox_Flux(c,Mg_tot,Pi_e,ADP_e);

% Fixed parameters
RT = 2.5775;                    % kJ  mol^{-1}
F = 0.096484;                   % kJ mol^{-1} mV^{-1}
dG_C4o = -122.94;               % kJ mol^{-1}
pH_e = 7.1; 
H_e = 10^(-pH_e);               % Molar
Ctot   = 2.70e-3;

% Adjustable parameters:
x_C4   = 6.766e-5;

k_O2   = 120e-6;

% Concentration state variables:
H_x    = c(1);
Cred   = c(6);
O2     = c(7);
dPsi   = c(19);
Cox    = Ctot - Cred;

% Mg, K, H in IM space
H_i  = H_e;

% Membrane proton motive force:
dG_H = F*dPsi + (RT)*log(H_i/H_x);

dG_C4op = dG_C4o - 2*RT*log(H_x/1e-7);

% Complex IV flux
J_C4 = x_C4*(O2/(O2+k_O2))*(Cred/Ctot)*( exp(-(dG_C4op+2*dG_H)/(2*RT))*Cred*(O2^0.25) - Cox*exp(F*dPsi/RT) );
