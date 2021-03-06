function [SSE_total] = calculate_SSE(X_activities, plotting)

%% Adjustable parameters


X_AtC  = 0;

% NEW values for OXPHOS
X_C1 = X_activities(1);
X_C3 = X_activities(2);
X_C4 = X_activities(3);
X_F1F0 = X_activities(4);
E_ANT = X_activities(5);
E_PiC = X_activities(6);

% E_ANT  = 0.065;
% E_PiC  = 1.0e6;

X_DH   = 0.0866; 

E_H = 2e3; % 1e3

flags = 1;

% Load Data 
load vinnakota.mat;

%% Parameters defining metabolite pools

% Volume fractions and water space fractions
V_c = 1.0;       % buffer volume fraction       % L buffer (L cuvette)^(-1)
V_m = 0.455e-3; %0.0005;       % mitochondrial volume fraction % L mito (L cuvette)^(-1)
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




%% Steady-state experiments

% Resetting steady state with Pi and ADP present in buffer. I am starting
% the amount of [Pi] that will result in roughly 5 mM [Pi] in the buffer
% in the resulting steady states.

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
X_AtC = (0:0.1:8)*1e-6;  % Increase max hydrolysis to find apparent Km.

opts = odeset('NonNeg',2:10);

% looping through different ATP consumptions states
for i = 1:length(X_AtC)
   pars = [X_DH; X_C1; X_C3; X_C4; X_F1F0; E_ANT; E_PiC; E_H; X_AtC(i)];
   % run for long time to acheive steady-state
   [t,x] = ode15s(@Mito_dXdT,[0 800],x0,opts,pars,conc,fixedpars,1);
   [f,J] = Mito_dXdT(0,x(end,:),pars,conc,fixedpars,2);
   J_C4 = J(10);

   % Units to match with Bazil
   O2Conv = 1/ 2 * V_mito_total / X_CS * 60 * 1e9;
   JO2_hi(i) = J_C4 / 2 * V_mito_total / X_CS * 60 * 1e9;

   ADP_hi(i) = x(end,9);
   DPsi_hi(i) = x(end,1);
   NADH_hi(i) = x(end,5);
   cred_hi(i) = x(end,7);
   Pi_c_hi(i) = x(end,10); 
end

%% Low Pi Concentration
Pi_c_0 = 1.5e-3;


x0 = [DPsi_0;
    ATP_x_0; ADP_x_0; Pi_x_0; NADH_x_0; QH2_x_0;
    cred_i_0;
    ATP_c_0; ADP_c_0; Pi_c_0;
    ];

% looping through different ATP consumptions states
for i = 1:length(X_AtC)
   pars = [X_DH; X_C1; X_C3; X_C4; X_F1F0; E_ANT; E_PiC; E_H; X_AtC(i)];
   % run for long time to acheive steady-state
   [t,x] = ode15s(@Mito_dXdT,[0 800],x0,opts,pars,conc,fixedpars,1);
   [f,J] = Mito_dXdT(0,x(end,:),pars,conc,fixedpars,2);
   J_C4 = J(10);

   % Units to match with Bazil
   O2Conv = 1/ 2 * V_mito_total / X_CS * 60 * 1e9;
   JO2_low(i) = J_C4 / 2 * V_mito_total / X_CS * 60 * 1e9;

   ADP_low(i) = x(end,9);
   DPsi_low(i) = x(end,1);
   NADH_low(i) = x(end,5);
   cred_low(i) = x(end,7);
   Pi_c_low(i) = x(end,10); 
end

%% Residuals and Fitting
% Calculate residuals for ADP
f=fit(ADP_hi_pi, JO2_hi_pi,'linearinterp');
ADP_range = ADP_hi((ADP_hi*1e3 <= max(ADP_hi_pi)));
linear_int_JO2 = f(ADP_range.*1e3)';
resid_ADP = JO2_hi(ADP_hi*1e3 <= max(ADP_hi_pi)) - linear_int_JO2;
%resid_ADP = ADP_hi_pi- interp1(JO2, ADP*1e3,JO2_hi_pi );

% % Remove NaN;
% resid_ADP = resid_ADP(~isnan(resid_ADP));
SSE_adp = sum((resid_ADP).^2);


% Calculate Residuals for NADH_rel
% resid_NADH_nonrm = JO2_hi_pi - interp1(NADH/NAD_tot,JO2,NADH_hi_pi);

% resid_NADH_nonrm = resid_NADH_nonrm(~isnan(resid_NADH_nonrm)); % Remove NaN
% SSE_NADH_norm = sum((resid_NADH_nonrm).^2);

SSE_total = SSE_adp ;% + SSE_NADH_norm;


if plotting == 1
    figure();
    subplot(2,2,1); plot(JO2_hi, NADH_hi/NAD_tot); hold on;
    plot(JO2_hi_pi, NADH_hi_pi, '+');
    plot(JO2_low, NADH_low/NAD_tot); 
    plot(JO2_low_pi, NADH_low_pi, '+');
    xlabel('J_{O2}');
    ylabel('NADH (normalized)');
    legend('5 mM Pi','5 mM Pi Data','1 mM Pi','1 mM Pi Data');

    subplot(2,2,2);plot(JO2_hi,DPsi_hi*1e3); hold on;
    plot(JO2_hi_pi_few, deltaPsi_hi_pi,'+');
    plot(JO2_low,DPsi_low*1e3);
    plot(JO2_low_pi_few, deltaPsi_low_pi,'+');
    ylabel('DPsi (mV)');
    xlabel('J_{O2}');

    subplot(2,2,3);plot(JO2_hi, cred_hi/c_tot); hold on;
    plot(JO2_hi_pi_few, cytC_hi_pi,'+');
    plot(JO2_low, cred_low/c_tot);
    plot(JO2_low_pi_few, cytC_low_pi,'+');
    xlabel('J_{O2}');
    ylabel('Cytochrome C_{red} (normalized)');

    subplot(2,2,4); plot(ADP_hi*1e3, JO2_hi); hold on;
    plot(ADP_hi_pi,JO2_hi_pi, '+');
    plot(ADP_low*1e3, JO2_low);
    plot(ADP_low_pi,JO2_low_pi, '+');
    ylabel('J_{O2}');
    xlabel('ADP (mM)');
end


end

