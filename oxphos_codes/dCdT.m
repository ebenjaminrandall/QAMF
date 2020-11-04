function [f] = dCdT(t,c,Mg_tot,Pi_e,ADP_e);

% FUNCTION DCDT1: computes the derivatives for an ODE
%   model of oxidative phosphorylation.
%
% INPUTS: 
%   t - time (not used)
%   c - vector of state variables
%   x - vector of free parameters
%   Mg - buffer (free) magnesium concentration
%   Pi_e - buffer Pi concentration
%   ADP_e - buffer ADP concentration
%
% OUTPUTS:
%   f - vector of time derivatives of concentrations
%%

% Fixed parameters
RT = 2.5775;                    % kJ  mol^{-1}
F = 0.096484;                   % kJ mol^{-1} mV^{-1}
dG_C1o = -69.37;                % kJ mol^{-1}
dG_C3o = -32.53;                % kJ mol^{-1}
dG_C4o = -122.94;               % kJ mol^{-1}
dG_F1o = 36.03;                 % kJ mol^{-1}
n_A = 3.0;
pH_e = 7.1; 
H_e = 10^(-pH_e);               % Molar
K_e = 150.0e-3;                 % Molar
ATP_e = 0;                      % Molar
AMP_e = 0;                      % Molar
k_dHPi  = 10^(-6.75);           % H/Pi binding constant Molar
k_dHatp = 10^(-6.48);           % H/ATP binding constant Molar
k_dHadp = 10^(-6.29);           % H/ADP binding constant Molar
K_DT = 24e-6;                   % Mg/ATP binding constant Molar
K_DD = 347e-6;                  % Mg/ADP binding constant Molar
K_AK = 0.4331;                  % nondimensional
W_m = 0.143/0.20;               % mitochondrial water space (ml water per ml mito)
W_x = 0.9*W_m;                  % Matrix water space
W_i = 0.1*W_m;                  % IM water space
gamma = 5.99;                   % mito outer membrane area per unit volume micron^{-1}
Ctot   = 2.70e-3;
Qtot   = 1.35e-3;               
NADtot = 2.97e-3;               % M; total NAD/NADH

% Adjustable parameters:
k_Pi1  = 0.1553e-3;
k_Pi2  = 0.8222e-3;
k_Pi3  = 0.3601e-3;
k_Pi4  = 5.9242e-3;
k_PiH  = 0.2542e-3;
r      = 4.559;
x_DH   = 0.0866;
x_C1   = 4405; 
x_C3   = 4.887;
x_C4   = 6.766e-5;
x_F1   = 1000.0;
x_ANT  = 0.008123;
x_Pi1  = 3.850e5; 
x_KH   = 5.651e7;
x_Hle  = 200.0;
x_K    = 0;

k_O2   = 120e-6;
k_mADP = 3.50e-6;
x_buff = 100;
x_AK   = 0e6;
x_MgA  = 1e6;
x_A    = 85.0;                  % micron sec^{-1}
x_Pi2  = 327;                   % micron sec^{-1}

% Concentration state variables:
H_x    = c(1);
K_x    = c(2);
Mg_x   = c(3);
NADH_x = c(4);
QH2    = c(5);
Cred   = c(6);
O2     = c(7);
ATP_x  = c(8);
ADP_x  = c(9);
ATP_mx = c(10);
ADP_mx = c(11);
Pi_x   = c(12);
ATP_i  = c(13);
ADP_i  = c(14);
AMP_i  = c(15);
ATP_mi = c(16);
ADP_mi = c(17);
Pi_i   = c(18);
dPsi   = c(19);

% Other concentrations computed from the state variables:
NAD_x  = NADtot - NADH_x;
Q      = Qtot - QH2;
Cox    = Ctot - Cred;
ATP_fx = ATP_x - ATP_mx;
ADP_fx = ADP_x - ADP_mx;
ATP_fi = ATP_i - ATP_mi;
ADP_fi = ADP_i - ADP_mi;

% ADP/Mg binding in E space
ADP_me = ( (K_DD+ADP_e+Mg_tot) - sqrt((K_DD+ADP_e+Mg_tot)^2-4*(Mg_tot*ADP_e)) )/2;
ADP_fe = ADP_e - ADP_me;
Mg_e   = Mg_tot - ADP_me;

% Mg, K, H in IM space
H_i  = H_e;
Mg_i = Mg_e;
K_i  = K_e;

% Membrane proton motive force:
dG_H = F*dPsi + (RT)*log(H_i/H_x);

dG_C1op = dG_C1o - 1*RT*log(H_x/1e-7);
dG_C3op = dG_C3o + 2*RT*log(H_x/1e-7);
dG_C4op = dG_C4o - 2*RT*log(H_x/1e-7);

% Fluxes:
J_DH = x_DH*(r*NAD_x-NADH_x)*((1+Pi_x/k_Pi1)/(1+Pi_x/k_Pi2));
J_C1 = x_C1*( exp(-(dG_C1op+4*dG_H)/RT)*NADH_x*Q - NAD_x*QH2 );
J_C3 = x_C3*((1+Pi_x/k_Pi3)/(1+Pi_x/k_Pi4))*(exp(-(dG_C3op+4*dG_H-2*F*dPsi)/(2*RT))*Cox*sqrt(QH2) - Cred*sqrt(Q) );
J_C4 = x_C4*(O2/(O2+k_O2))*(Cred/Ctot)*( exp(-(dG_C4op+2*dG_H)/(2*RT))*Cred*(O2^0.25) - Cox*exp(F*dPsi/RT) );
J_F1   = x_F1*( exp(-(dG_F1o-n_A*dG_H)/RT)*(K_DD/K_DT)*ADP_mx*Pi_x - ATP_mx );
%J_F1   = x_F1*( exp(-(dG_F1o-n_A*dG_H)/RT)*ADP_x*Pi_x - ATP_x );
Psi_x = -0.65*dPsi;
Psi_i = +0.35*dPsi;
if (ADP_fi > 1e-12) || (ATP_fi > 1e-12)
  J_ANT  = x_ANT*( ADP_fi/(ADP_fi+ATP_fi*exp(-F*Psi_i/RT)) - ADP_fx/(ADP_fx+ATP_fx*exp(-F*Psi_x/RT)) )*(ADP_fi/(ADP_fi+k_mADP));
else
  J_ANT  = 0;
end  
H2PIi = Pi_i*H_i/(H_i+k_dHPi); H2PIx = Pi_x*H_x/(H_x+k_dHPi);
J_Pi1 = x_Pi1*(H_i*H2PIi - H_x*H2PIx)/(H2PIi+k_PiH);
J_Hle  = x_Hle*dPsi*(H_i*exp(F*dPsi/RT)-H_x)/(exp(F*dPsi/RT)-1);
J_KH   = x_KH*( K_i*H_x - K_x*H_i );
J_K    = x_K*dPsi*(K_i*exp(F*dPsi/RT)-K_x)/(exp(F*dPsi/RT)-1);
J_AKi  = x_AK*( K_AK*ADP_i*ADP_i - AMP_i*ATP_i );
J_AMP  = gamma*x_A*(AMP_e-AMP_i);
J_ADP  = gamma*x_A*(ADP_e-ADP_i);
J_ATP  = gamma*x_A*(ATP_e-ATP_i);
J_Pi2  = gamma*x_Pi2*(Pi_e-Pi_i);    
J_MgATPx = x_MgA*(ATP_fx*Mg_x-K_DT*ATP_mx);
J_MgADPx = x_MgA*(ADP_fx*Mg_x-K_DD*ADP_mx);
J_MgATPi = x_MgA*(ATP_fi*Mg_i-K_DT*ATP_mi);
J_MgADPi = x_MgA*(ADP_fi*Mg_i-K_DD*ADP_mi);

% Derivatives:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Model for ATP, ADP, and Pi proton buffering
% dATPdT = (+J_F1 - J_ANT - J_MgATPx)/W_x; % Rate of production of Mg-unbound ATP in matrix
% dADPdT = (-J_F1 + J_ANT - J_MgADPx)/W_x; % Rate of production of Mg-unbound ADP in matrix
% dPIdT1 = (-J_F1 + J_Pi1)/W_x; % Rate of production of Pi in matrix
% J_H   = (-4*J_C1 - 2*J_C3 - 4*J_C4 + (n_A-1)*J_F1 + 2*J_Pi1 + J_Hle - J_KH)/W_x;
% num = -dATPdT/(1+k_dHatp/H_x) - dADPdT/(1+k_dHadp/H_x) - dPIdT1/(1+k_dHPi/H_x) + J_H;
% den = 1 + (ATP_fx*k_dHatp/H_x^2)/((1+k_dHatp/H_x)^2) + (ADP_fx*k_dHadp/H_x^2)/((1+k_dHadp/H_x)^2) + (Pi_x*k_dHPi/H_x^2)/((1+k_dHPi/H_x)^2);
% f(1) = num/den;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f(1)  = x_buff*H_x*( -4*J_C1 - 2*J_C3 - 4*J_C4 + (n_A-1)*J_F1 + 2*J_Pi1 + J_Hle - J_KH )/W_x; % H_x
f(2)  = (J_KH + J_K)/W_x; % K_x
f(3)  = (-J_MgATPx - J_MgADPx)/W_x; % Mg_x
f(4)  = (+J_DH - J_C1)/W_x; % NADH
f(5)  = (+J_C1 - J_C3)/W_x; % QH2
f(6)  = (+2*J_C3 - 2*J_C4)/W_i; % Cred
f(7)  = 0; % O2
f(8)  = (+J_F1 - J_ANT)/W_x; % ATP_x
f(9)  = (-J_F1 + J_ANT)/W_x; % ADP_x
f(10) = (J_MgATPx)/W_x; % ATP_mx
f(11) = (J_MgADPx)/W_x; % ADP_mx
f(12) = (-J_F1 + J_Pi1 )/W_x;  % Pi_x
f(13) = (+J_ATP + J_ANT + J_AKi )/W_i; % ATP_i
f(14) = (+J_ADP - J_ANT - 2*J_AKi )/W_i; % ADP_i
f(15) = (+J_AMP + J_AKi)/W_i; % AMP_i
f(16) = (J_MgATPi)/W_i; % ATP_mi
f(17) = (J_MgADPi)/W_i; % ADP_mi
f(18) = (-J_Pi1 + J_Pi2 )/W_i; % Pi_i
f(19) = (1.48e5)*( 4*J_C1 + 2*J_C3 + 4*J_C4 - n_A*J_F1 - J_ANT - J_Hle - J_K );
f = f';
