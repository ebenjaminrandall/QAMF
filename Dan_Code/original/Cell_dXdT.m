function [f] = Cell_dXdT(~,x,x_AtC)
% This function is used to calculate time derivatives of state
% variables of the cell-level cardaic energetics.

% Copyright: 
% Fan Wu and Daniel. A Beard, 2008
% Computational Engineering Group
% Biotechnology and Bioengineering Center
% Medical College of Wisconsin


%% Setting adjustable parameter values

% Parameters for oxidative phosphorylation
x_C1 = 32365;
x_C3 = 0.79081;
x_C4 = 0.00010761;
x_F1 = 8120.4;
x_ANT = 0.006762;
x_PI1 = 3.3356e+07;
k_PIH = 0.0016125;
x_Hle = 50;
k_PI1 = 3.2489e-05;
k_PI2 = 0.0034745;


%% define indices for all state variables
% (i) Matrix species and dPsi
idPsi        = 1;
iATP_x       = 2;
iADP_x       = 3;
iPI_x        = 4;
iNADH_x      = 5;
iQH2_x       = 6;
 
%  (ii) IM space species
iCred_i      = 7;
iATP_i       = 8;
iADP_i       = 9;
iAMP_i       = 10;
iPI_i        = 11;
 
%  (iii) Cytoplasmic species
iATP_c       = 12;
iADP_c       = 13;
iPI_c        = 14;
iPCr_c       = 15;
iAMP_c       = 16;
iCr_c        = 17;

   
%% Defining indices of reactants
iATP    =   1;
iADP    =   2;
iAMP    =   3;
iPI     =   4;
iNADH   =   5;
iNAD    =   6;
iQH2    =   7;
iCOQ    =   8;
iCox    =   9;
iCred   =   10;
iO2     =   11;
iH2O    =   12;
iCO2tot =   13;
iPCr    =   14;
iCr     =   15;

N_reactant = 15;

%%  Listing fixed model parameters
%  (i) Thermochemical constants
RT = 8.314*(37+273.15)/1e3;          % kJ  mol^{-1}
F = 0.096484;                   % kJ mol^{-1} mV^{-1}

%  (ii) Subcelular volumes and water spaces
Vmito = 0.2882;                 % (ml mito / ml cell)
Vcyto = 0.6601;                 % (ml cyto / ml cell)
Rm_cyto = Vmito / Vcyto;        % Volume ratio mito volume / cytoplasm volume
% Rc_cell = Vcyto;                % Volume ratio of cytosol / cell volume
W_c = 0.807*1.044;              % cytosol water space (ml water per ml cytosol) [VB2002]
W_m = 0.664*1.09;               % mitochondrial water space (ml water per ml mito) [VB2002]
W_x = 0.9*W_m;                  % Matrix water space
W_i = 0.1*W_m;                  % IM water space

%  (iii) Pooled concentrations
Ctot   = 2.70e-3;               % M; total cytoC
Qtot   = 1.35e-3;               % M; total Q + QH2
NADtot = 2.97e-3;               % M; total NAD+NADH

%  (iv) Ox Phos Model Parameters
n_A     = 3.0;  % unitless
k_O2   = 1.2e-4;

%  (v) Outer membrane transport parameters
x_A    = 85;        % micron sec^{-1}
x_PI2  = 327;       % micron sec^{-1}
gamma  = 5.99;      % mito membrane area per cell volume micron^{-1}

% (vi) Oxygen solubility
a_3 = 1.74e-6;                  % oxygen solubility in cell

%% Loading values of state variables
MinCon = 1e-32;

% (i) Matrix species and dPsi
dPsi        = x(idPsi);
ATP_x       = x(iATP_x);
ADP_x       = x(iADP_x);
PI_x        = x(iPI_x);
NADH_x      = x(iNADH_x);
QH2_x       = x(iQH2_x);

% (ii) IM space species
Cred_i      = x(iCred_i);
ATP_i       = x(iATP_i);
ADP_i       = x(iADP_i);
AMP_i       = x(iAMP_i);
PI_i        = x(iPI_i);

% (iii) Cytoplasmic species
ATP_c       = x(iATP_c);
ADP_c       = x(iADP_c);
PI_c        = x(iPI_c);
PCr_c       = x(iPCr_c);
AMP_c       = x(iAMP_c);
Cr_c        = x(iCr_c);

% (iv) Other concentrations computed from the state variables:
NAD_x  = NADtot - NADH_x;
COQ_x  = Qtot - QH2_x;
Cox_i  = Ctot - Cred_i;

% (v) set the H+, Mg2+, and K+ concentrations
H_x  = 10^(-7.4);
K_x  = 100e-3;
Mg_x = 1e-3;
H_c  = 10^(-7.15);
K_c  = 140e-3;
Mg_c  = 1e-3;
H_i  = H_c;
Mg_i = Mg_c;
K_i  = K_c;

% Oxygen concentrations:
PO2 = 25; % mmHg
CfO2 = a_3*PO2;

%% Loading thermodynamic data (deltG, pK, etc.)
% T = 298.15K (25 C)    I = 0.17 M
% standard Gibbs free energy of formation of reference species, (kJ/mol)

% without temperature correction on dGf
dGf1(1:N_reactant) = 0;
dGf1(iH2O) = -235.74;        % H2O
dGf1(iO2) = 16.40;           % O2(aq)
dGf1(iNADH) = 39.31;         % NADH
dGf1(iNAD) = 18.10;          % NAD+
dGf1(iQH2) = -23.30;         % QH2
dGf1(iCOQ) = 65.17;          % Q
dGf1(iATP) = -2771.00;       % ATP4-
dGf1(iADP) = -1903.96;       % ADP3-
dGf1(iAMP) = -1034.66;       % AMP2-
dGf1(iCred) = -27.41;        % CytoC(red)2+
dGf1(iCox) = -6.52;          % CytoC(ox)3+
dGf1(iPI) = -1098.27;        % HPO42-
dGf1(iPCr) = 0;              % PCr2-
dGf1(iCr) = -252.68;         % HCr
dGf1(iCO2tot) = -530.71;     % CO2tot

% K values for reference species 
% pK_KATP is corrected to be 1.013, 08/26/08
% pK_KADP is corrected to be 0.882, 08/26/08
% pK_KAMP is corrected to be 0.6215, 08/26/08
% pK_MgOAA is corrected to be 0.8629, 08/26/08
% pK_KSUC is corrected to be 0.3525, 08/26/08
Kh(1:N_reactant) = inf; Km(1:N_reactant) = inf; Kk(1:N_reactant) = inf;
Kh(iATP) = 10^(-6.59); Km(iATP) = 10^(-3.82); Kk(iATP) = 10^(-1.013); 
Kh(iADP) = 10^(-6.42); Km(iADP) = 10^(-2.79); Kk(iADP) = 10^(-0.882); 
Kh(iAMP) = 10^(-6.22); Km(iAMP) = 10^(-1.86); Kk(iAMP) = 10^(-0.6215); 
Kh(iPI) = 10^(-6.71); Km(iPI) = 10^(-1.69); Kk(iPI) = 10^(+0.0074);
Kh(iCO2tot) = 10^(-9.82);

% compute binding polynomials for reactants
P_x(1:N_reactant) = 1; P_c(1:N_reactant) = 1; P_i(1:N_reactant) =1;
P_x(iATP) = 1 + H_x/Kh(iATP) + Mg_x/Km(iATP) + K_x/Kk(iATP);
P_c(iATP) = 1 + H_c/Kh(iATP) + Mg_c/Km(iATP) + K_c/Kk(iATP);
P_i(iATP) = 1 + H_i/Kh(iATP) + Mg_i/Km(iATP) + K_i/Kk(iATP);
P_x(iADP) = 1 + H_x/Kh(iADP) + Mg_x/Km(iADP) + K_x/Kk(iADP);
P_c(iADP) = 1 + H_c/Kh(iADP) + Mg_c/Km(iADP) + K_c/Kk(iADP);
P_i(iADP) = 1 + H_i/Kh(iADP) + Mg_i/Km(iADP) + K_i/Kk(iADP);
P_x(iAMP) = 1 + H_x/Kh(iAMP) + Mg_x/Km(iAMP) + K_x/Kk(iAMP);
P_c(iAMP) = 1 + H_c/Kh(iAMP) + Mg_c/Km(iAMP) + K_c/Kk(iAMP);
P_i(iAMP) = 1 + H_i/Kh(iAMP) + Mg_i/Km(iAMP) + K_i/Kk(iAMP);
P_x(iPI) = 1 + H_x/Kh(iPI) + Mg_x/Km(iPI) + K_x/Kk(iPI); % add K-bound item, 06/10/08
P_c(iPI) = 1 + H_c/Kh(iPI) + Mg_c/Km(iPI) + K_c/Kk(iPI); % add K-bound item, 06/10/08
P_i(iPI) = 1 + H_i/Kh(iPI) + Mg_i/Km(iPI) + K_i/Kk(iPI); % add K-bound item, 06/10/08
P_x(iCO2tot) = 1 + H_x/Kh(iCO2tot);
P_i(iCO2tot) = 1 + H_i/Kh(iCO2tot);
P_c(iCO2tot) = 1 + H_c/Kh(iCO2tot);


%% I. Flux expresssions in the oxidative phosphorylation 

%% ------------------------------- 
% Substrate/ion transport
% Transport of ADP, ATP, and PI across outer membrane:
J_ADP   = gamma*x_A*(ADP_c-ADP_i);
J_ATP   = gamma*x_A*(ATP_c-ATP_i);
J_AMP   = gamma*x_A*(AMP_c-AMP_i);
J_PI2   = gamma*x_PI2*(PI_c-PI_i);    


%% -------------------------------

% 0. Dehydrogenase
r      = 4.559;
x_DH   = 0.0866;
k_Pi1  = 0.1553e-3;
k_Pi2  = 0.8222e-3;
J_DH = x_DH*(r*NAD_x-NADH_x)*((1+PI_x/k_Pi1)/(1+PI_x/k_Pi2));

% 1. Complex I
% NADH_x + Q_x + 5H+_x <-> NAD+_x + QH2_x + 4H+_i + 4dPsi

% compute free energy of reaction from free energe of formation
dGr_C1o = dGf1(iNAD) + dGf1(iQH2) - dGf1(iNADH) - dGf1(iCOQ);
% compute Keq from combined free energy of reactions (including potential
% change due to charge translocation)
Keq_C1 = exp(-(dGr_C1o + 4*F*dPsi)/RT);
% compute Kapp from Kapp = Keq*H_x^n/H_i^m*Product(Poly)  (n,m are
% stochi.coefficients)
Kapp_C1 = Keq_C1*H_x^5/H_i^4;
J_C1 = x_C1*( Kapp_C1*NADH_x*COQ_x - NAD_x*QH2_x );

%% -------------------------------
% 2. Complex III
% QH2_x + 2cytoC(ox)3+_i + 2H+_x <-> Q_x + 2cytoC(red)2+_i + 4H+_i + 2dPsi

dGr_C3o = dGf1(iCOQ) + 2*dGf1(iCred) - dGf1(iQH2) - 2*dGf1(iCox);
Keq_C3 = exp(-(dGr_C3o + 2*F*dPsi)/RT);
Kapp_C3 = Keq_C3*H_x^2/H_i^4;
QH2_x = max(MinCon,QH2_x); COQ_x = max(MinCon,COQ_x);
J_C3 = x_C3*((1+PI_x/k_PI1)/(1+PI_x/k_PI2))*...
          (Kapp_C3^0.5*Cox_i*sqrt(QH2_x) - Cred_i*sqrt(COQ_x) );

%% -------------------------------
% 3. Complex IV
% 2cytoC(red)2+_i + 0.5O2_x + 4H+_x <-> 2cytoC(ox)3+_x + H2O_x + 2H+_i +
% 2dPsi

dGr_C4o = 2*dGf1(iCox) + dGf1(iH2O) - 2*dGf1(iCred) - 0.5*dGf1(iO2);
% 2 charges from translocation of proton, and the other 2 from cytoC
Keq_C4 = exp(-(dGr_C4o + 4*F*dPsi)/RT);
Kapp_C4 = Keq_C4*H_x^4/H_i^2;
J_C4 = x_C4*(CfO2/(CfO2+k_O2))*exp(F*dPsi/RT)*(Cred_i/Ctot)*( Kapp_C4^0.5*Cred_i*(CfO2^0.25) - Cox_i );


%% -------------------------------
% 4. F1Fo-ATPase
% ADP3-_x + HPO42-_x + H+_x + n_A*H+_i <-> ATP4- + H2O + n_A*H+_x

dGr_F1o = dGf1(iATP) + dGf1(iH2O) - dGf1(iADP) - dGf1(iPI);
Keq_F1 = exp(-(dGr_F1o-n_A*F*dPsi)/RT);
Kapp_F1 = Keq_F1*H_i^n_A/H_x^(n_A-1)*P_x(iATP)/(P_x(iADP)*P_x(iPI));
J_F1 = x_F1*(Kapp_F1*ADP_x*PI_x - ATP_x);


%% -------------------------------
% 5. ANT
% ATP4-_x + ADP3-_i <-> ATP4-_i + ADP3-_x

ADP_i1 = ADP_i/P_i(iADP); % ADP^3-
ATP_i1 = ATP_i/P_i(iATP); % ATP^4-
ADP_x1 = ADP_x/P_x(iADP); % ADP^3-
ATP_x1 = ATP_x/P_x(iATP); % ATP^4-

del_D = 0.0167;
del_T = 0.0699;
k2_ANT = 9.54/60; % = 1.59e-1
k3_ANT = 30.05/60; % = 5.01e-1
K_D_o_ANT = 38.89e-6; 
K_T_o_ANT = 56.05e-6;
A = +0.2829;
B = -0.2086;
C = +0.2372;
fi = F*dPsi/RT;
k2_ANT_fi = k2_ANT*exp((A*(-3)+B*(-4)+C)*fi);
k3_ANT_fi = k3_ANT*exp((A*(-4)+B*(-3)+C)*fi);

K_D_o_ANT_fi = K_D_o_ANT*exp(3*del_D*fi);
K_T_o_ANT_fi = K_T_o_ANT*exp(4*del_T*fi);
q = k3_ANT_fi*K_D_o_ANT_fi*exp(fi)/(k2_ANT_fi*K_T_o_ANT_fi);
term2 = k2_ANT_fi*ATP_x1.*ADP_i1*q/K_D_o_ANT_fi ;
term3 = k3_ANT_fi.*ADP_x1.*ATP_i1/K_T_o_ANT_fi;
num = term2 - term3;
den = (1 + ATP_i1/K_T_o_ANT_fi + ADP_i1/K_D_o_ANT_fi)*(ADP_x1 + ATP_x1*q);
J_ANT = x_ANT/7.2679e-003*(0.70e-1)*num/den; % x_ANT'(in the paper) = x_ANT/7.2679e-003*(0.70e-1)

%% -------------------------------
% 6. H+-PI2 cotransporter

H2PIi1 = PI_i*(H_i/Kh(iPI))/P_i(iPI);
H2PIx1 = PI_x*(H_x/Kh(iPI))/P_x(iPI);
J_PI1 = x_PI1/k_PIH*(H_i*H2PIi1 - H_x*H2PIx1)/(1+H2PIi1/k_PIH)/(1+H2PIx1/k_PIH);


%% -------------------------------
% 7. H+ leak
if abs(dPsi) > 1e-9
  J_Hle = x_Hle*dPsi*(H_i*exp(F*dPsi/RT)-H_x)/(exp(F*dPsi/RT)-1);
else
  J_Hle = x_Hle*RT*( H_i - H_x )/F;
end


% %% -------------------------------
% % 8. K+/H+ anti-porter
% J_KH   = x_KH*( K_i*H_x - K_x*H_i);
% 


%% ---------------------------------------
% 10. Creatine kinase reaction
% ADP3- + PCr2- + H+ = ATP4- + Cr0
% set Cr and PCr concentrations according to experiments

x_CK = 1e7;
% K_CK = exp(50.78/RT); 
K_CK = 7.408e8; 
Kapp_CK = K_CK*H_c*P_c(iATP)*P_c(iCr)/P_c(iADP)/P_c(iPCr);
J_CKe  = x_CK * (Kapp_CK*ADP_c*PCr_c - ATP_c*Cr_c );


%% ------------------------------------------
% 11. Adenylate kinase reaction
% 2ADP3- = ATP4- + AMP2-
% dGr_AKo = dGf1(iATP) + dGf1(iAMP) - 2*dGf1(iADP);
% Keq_AK = exp(-dGr_AKo/RT);
Keq_AK = 3.97e-1;
Kapp_AK = Keq_AK*P_c(iATP)*P_c(iAMP)/(P_c(iADP)*P_c(iADP));
x_AK = 1e7;
J_AKi  = 0*x_AK*( Kapp_AK*ADP_i*ADP_i - AMP_i*ATP_i );
J_AKe  = x_AK*( Kapp_AK*ADP_c*ADP_c - AMP_c*ATP_c );

%% -------------------------------------------
% 12. ATPase flux
% ATP4- + H2O = ADP3- + PI2- + H+
J_AtC  = x_AtC; 


%% ---------------------------------------------------
%% Computing time derivatives of state variables

% Time derivatives:
  
%  (i) Matrix species and dPsi
f(idPsi)   = (1.48e5)*( 4*J_C1 + 2*J_C3 + 4*J_C4 - n_A*J_F1 - J_ANT - J_Hle );
f(iATP_x)  = (+ J_F1 - J_ANT)/W_x; 
f(iADP_x)  = (- J_F1 + J_ANT)/W_x;
f(iPI_x)   = (- J_F1 + J_PI1)/W_x;
f(iNADH_x) = (+J_DH - J_C1)/W_x; % NADH
f(iQH2_x)  = ( + J_C1 - J_C3)/W_x; 

%  (ii) IM space species
f(iCred_i) = (+2*J_C3 - 2*J_C4)/W_i; 
f(iATP_i) = (J_ATP + J_ANT + J_AKi)/W_i; 
f(iADP_i) = (J_ADP - J_ANT - 2*J_AKi)/W_i; 
f(iAMP_i) = (J_AMP + J_AKi)/W_i;   
f(iPI_i)  = (-J_PI1 + J_PI2)/W_i; 

%  (iii) Buffer species
f(iATP_c) = (-Rm_cyto*J_ATP - J_AtC + J_CKe + J_AKe)/W_c;  
f(iADP_c) = (-Rm_cyto*J_ADP + J_AtC - J_CKe - 2*J_AKe)/W_c; 
f(iPI_c)  = (-Rm_cyto*J_PI2 + J_AtC)/W_c; 
f(iPCr_c) = (-J_CKe)/W_c;
f(iAMP_c) = (-Rm_cyto*J_AMP + J_AKe)/W_c;
f(iCr_c)  = (+J_CKe)/W_c;

f = f';
