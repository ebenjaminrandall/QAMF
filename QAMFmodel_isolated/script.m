clear
clc;
close all;

flags = 1;

%% Parameters defining metabolite pools

% Volume fractions and water space fractions
V_c = 1.0;       % buffer volume fraction       % L buffer (L cuvette)^(-1)
V_m = 0.0005;       % mitochondrial volume fraction % L mito (L cuvette)^(-1)
R_m2c = V_m / V_c;  % mito to cyto volume ratio     % L mito (L cuvette)^(-1)
W_c = 1.0;       % buffer water space           % L buffer water (L buffer)^(-1)
W_m = 0.7238;       % mitochondrial water space     % L mito water (L mito)^(-1)
W_x = 0.9*W_m;      % matrix water space            % L matrix water (L mito)^(-1)
W_i = 0.1*W_m;      % intermembrane water space     % L IM water (L mito)^(-1)

V_cuvette = 2e-3; % volume of experimental cuvette % mL
V_mito_total = V_m * V_cuvette; % total volume of mito in cuvette % L
X_CS = 0.67; % activity of citrate synthase per 2 mL cuvette

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

fixedpars.pools.NAD_tot = NAD_tot;
fixedpars.pools.Q_tot = Q_tot;
fixedpars.pools.c_tot = c_tot;

% Membrane capacitance
Cm = 3e-6;

fixedpars.capacitance.Cm = Cm;

%% Set fixed pH, cation concentrations, and O2 partial pressure

% pH
pH_x = 7.4;
pH_c = 7.2;

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
X_AtC  = 0;

% NEW values for OXPHOS
X_C1   = 3e5;
X_C3   = 5e7;
X_C4   = 0.20;
X_F1F0 = 100;
E_ANT  = 0.065;
E_PiC  = 1.0e6;
E_H = 2e3; % 1e3

load('params.mat')
% NEW values for OXPHOS
X_C1 = X_fitted(1);
X_C3 = X_fitted(2);
X_C4 = X_fitted(3);
X_F1F0 = X_fitted(4);
% E_ANT = X_fitted(5);
% E_PiC = X_fitted(6);
%% Initial conditions

% Membrane potential
DPsi_0    = 175*1e-3;

% Matrix species
ATP_x_0   = 0.5e-3;
ADP_x_0   = 9.5e-3;
Pi_x_0    = 0.3e-3;
NADH_x_0  = 0.1 * NAD_tot;
QH2_x_0   = 0.1 * Q_tot;

% IM species
cred_i_0  = 0.1 * c_tot;

% Cytosol species
ATP_c_0   = 0e-3;
ADP_c_0   = 0e-3;
Pi_c_0    = 5.0e-3;

x0 = [DPsi_0;
    ATP_x_0; ADP_x_0; Pi_x_0; NADH_x_0; QH2_x_0;
    cred_i_0;
    ATP_c_0; ADP_c_0; Pi_c_0;
    ];

%% Run state-1, -2, -3, -4 experiment

opts = odeset('NonNeg',2:10);
% State 1:
X_DH   = 0; % (no substrate/no dehydrogenase)
pars = [X_DH; X_C1; X_C3; X_C4; X_F1F0; E_ANT; E_PiC; E_H; X_AtC];
[t1,x1] = ode15s(@Mito_dXdT,[0 30],x0,opts,pars,conc,fixedpars,1);
% State 2:
X_DH   = 0.0866; % (turn on dehydrogenase)
pars = [X_DH; X_C1; X_C3; X_C4; X_F1F0; E_ANT; E_PiC; E_H; X_AtC];
[t2,x2] = ode15s(@Mito_dXdT,[30 90],x1(end,:),opts,pars,conc,fixedpars,1);
% State 3 and 4:
x0 = x2(end,:);
ADP_injection = 0.375e-3; % Molar
x0(9) = ADP_injection; % Molar
[t3,x3] = ode15s(@Mito_dXdT,[90 300],x0,opts,pars,conc,fixedpars,1);

t = [t1(1:end-1); t2(1:end-1); t3];
x = [x1(1:end-1,:); x2(1:end-1,:); x3];

NADrel = x(:,5)./NAD_tot;

% obtain Oxygen Flux
for i = 1:length(t)
   [f,J] = Mito_dXdT(0,x(i,:),pars,conc,fixedpars,2);
   J_C4 = J(10);
   % Units to match with Bazil
   % JO2(i) = J_C4 / 2 * V_mito_total / X_CS * 60 * 1e9;
   JO2(i) = J_C4 / 2 * V_mito_total; % mol/s
end

figure(1); plot(t,NADrel)
ylabel('NADH')

figure(2); plot(t,JO2)
ylabel('J_{O2} (mol/s)')

t_stop_2 = 90;
t_stop_3 = 180;
t_int = t((t_stop_2<= t) & (t<=t_stop_3));
JO2_int = JO2((t_stop_2<= t) & (t<=t_stop_3));

O2_integration = trapz(t_int, JO2_int);



%ATP_total = x(end,2) * V_mito_total * W_x + x(end,8) * V_cuvette;
ATP_total = ADP_injection*V_cuvette;
PO_ratio = ATP_total / (O2_integration*2);
disp('PO Ratio:');
disp(PO_ratio);


%% Electrode Smoothing

t_e = linspace(5,250,500);
delta_t = t_e(2) - t_e(1);
JO2_e(1) = 0;

for i = 2:length(t_e)
    % Manually integrating because MATLAB's dynamic step size fails here
    JO2_e(i) = dXdT_electrode(t_e(i), JO2_e(i-1), t, JO2) * delta_t+JO2_e(i-1);
end

hold on; plot(t_e,JO2_e);
legend('Numerical Simuation','Simulated Electrode Reading');
xlabel('Time (s)')
stage_2 = max(JO2_e((t_e>=60) & (t_e<=70)));
state_3 = max(JO2_e);
RCR = state_3/stage_2;
disp('RCR:');
disp(RCR);

%% Steady-state experiments

% Resetting steady state with Pi and ADP present in buffer. I am starting
% the amount of [Pi] that will result in roughly 5 mM [Pi] in the buffer
% in the resulting steady states.

%% Initial conditions

% Membrane potential
DPsi_0    = 175*1e-3;

% Matrix species
ATP_x_0   = 0.5e-3;
ADP_x_0   = 9.5e-3;
Pi_x_0    = 0.3e-3;
NADH_x_0  = 0.1 * NAD_tot;
QH2_x_0   = 0.1 * Q_tot;

% IM species
cred_i_0  = 0.1 * c_tot;

% Cytosol species
ATP_c_0   = 0e-3;
ADP_c_0 = 0.5e-3;
Pi_c_0 = 5.5e-3;

x0 = [DPsi_0;
    ATP_x_0; ADP_x_0; Pi_x_0; NADH_x_0; QH2_x_0;
    cred_i_0;
    ATP_c_0; ADP_c_0; Pi_c_0;
    ];

% range of ATP consumption rates
X_AtC = (0:0.1:6)*1e-6 / V_c;  % Increase max hydrolysis to find apparent Km.
%X_AtC = (0:0.6:50)*1e-6;

clear JO2
clear ADP
clear DPsi

opts = odeset('NonNeg',2:10);

% looping through different ATP consumptions states
for i = 1:length(X_AtC)
   pars = [X_DH; X_C1; X_C3; X_C4; X_F1F0; E_ANT; E_PiC; E_H; X_AtC(i)];
   % run for long time to acheive steady-state
   [t,x] = ode15s(@Mito_dXdT,[0 1600],x0,opts,pars,conc,fixedpars,1);
   [f,J] = Mito_dXdT(0,x(end,:),pars,conc,fixedpars,2);
   J_C4 = J(10);

   % Units to match with Bazil
   JO2(i) = J_C4 / 2 * V_mito_total / X_CS * 60 * 1e9;

   ADP(i) = x(end,9);
   DPsi(i) = x(end,1);
   NADH(i) = x(end,5);
   cred(i) = x(end,7);
   Pi_c(i) = x(end,10); 
end

figure(8);
plot(X_AtC * V_c*1e3,Pi_c*1e3); 
xlabel('X_AtC'); ylabel('Pi_c');


figure(3); plot(JO2,ADP*1e3)
xlabel('J_{O2}');
ylabel('ADP (mM)');

figure(4); plot(JO2,DPsi*1e3)
ylabel('DPsi (mV)');
xlabel('J_{O2}');

figure(5); plot(JO2, cred/c_tot)
xlabel('J_{O2}');
ylabel('Cytochrome C_{red} (normalized)');

figure(6); plot(JO2, NADH/NAD_tot)
xlabel('J_{O2}');
ylabel('NADH (normalized)');

subplot(2,2,1); figure(6); plot(JO2, NADH/NAD_tot)
xlabel('J_{O2}');
ylabel('NADH (normalized)');

subplot(2,2,2);plot(JO2,DPsi*1e3)
ylabel('DPsi (mV)');
xlabel('J_{O2}');

subplot(2,2,3);plot(JO2, cred/c_tot)
xlabel('J_{O2}');
ylabel('Cytochrome C_{red} (normalized)');

subplot(2,2,4); plot(ADP*1e3, JO2)
ylabel('J_{O2}');
xlabel('ADP (mM)');

% Calculate apparent Km (Note: need to increase max hydrolysis)
[C I] = min(abs(JO2 - max(JO2)/2));
disp('Apparent ADP Km (uM):');
disp(ADP(I)*1e6);


%% Data
MVO2_high_data = [14.509
    53.625
    77.07
    108.57
    144.36];

NADH_norm_high_data = [0.6642
    0.4821
    0.3901
    0.3688
    0.3695];

ADP_high_data = [0
    0.0049
    0.0578
    0.1648
    0.4716];

figure();
plot(ADP*1e3, JO2);
hold on;
plot(ADP_high_data,MVO2_high_data, 'o');
xlabel('ADP (mM)');
ylabel('VO2 nmoll O2 / min / U CS')