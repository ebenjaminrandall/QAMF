 clear 

%% Flags 
%{
This section sets certain flags within the code: 
    flags       determines whether the model solves for dXdT or solves the
                auxiliary equations 
    printon     prints the figures at the end of the the code
%}
flags   = 1;    % 1 - solves model for dXdT;    2 - solves for auxiliary equations
printon = 0;    % 0 - does not print figues;    1 - prints figures 

%% Fixed parameters and specifications 
%{
This section creates a structure containing the various volume and water
fractions, total metabolite pool concentrations, and the membrane
capacitance.
%}

% Volume fractions and water space fractions
V_c = 0.6601;       % cytosol volume fraction       % L cyto (L cell)^(-1)
V_m = 0.2882;       % mitochondrial volume fraction % L mito (L cell)^(-1)
V_m2c = V_m / V_c;  % mito to cyto volume ratio     % L mito (L cyto)^(-1)
W_c = 0.8425;       % cytosol water space           % L cyto water (L cyto)^(-1) 
W_m = 0.7238;       % mitochondrial water space     % L mito water (L mito)^(-1)
W_x = 0.9*W_m;      % matrix water space            % L matrix water (L mito)^(-1)
W_i = 0.1*W_m;      % intermembrane water space     % L IM water (L mito)^(-1)

fixedpars.fractions.V_c = V_c; 
fixedpars.fractions.V_m = V_m; 
fixedpars.fractions.V_m2c = V_m2c; 
fixedpars.fractions.W_c = W_c; 
fixedpars.fractions.W_m = W_m; 
fixedpars.fractions.W_x = W_x; 
fixedpars.fractions.W_i = W_i; 

% Total pool concentrations 
NAD_tot = 2.97e-3;  % NAD+ and NADH conc            % mol (L matrix water)^(-1)
Q_tot   = 1.35e-3;  % Q and QH2 conc                % mol (L matrix water)^(-1) 
c_tot   = 2.7e-3;   % cytochrome c ox and red conc  % mol (L IM water)^(-1)

Cr_tot  = 35e-3;    % creatine and creatine phosphate conc % mol (L cell)^(-1)
Cr_tot_c = Cr_tot / (V_c * W_c); % convert to mol (L cyto water)^(-1)

fixedpars.pools.NAD_tot = NAD_tot; 
fixedpars.pools.Q_tot = Q_tot; 
fixedpars.pools.c_tot = c_tot; 
fixedpars.pools.Cr_tot_c = Cr_tot_c; 

% Membrane capacitance 
Cm = 3e-3;       % mol (V * L mito)^(-1)

fixedpars.capacitance.Cm = Cm; 

%% Set pH, metabolite concentrations, and O2 partial pressure 
%{
This section creates a vector with the pH, K+ and Mg2+ concentrations, and
the O2 partial pressure.
%} 

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
PO2 = 25;           % mmHg

conc = [pH_x; pH_c; K_x; K_c; Mg_x; Mg_c; PO2]; 

%% Adjustable parameters 
%{
This section creates a vector of adjustable parameters, namely the flux
rates. 
%} 

X_DH  = 0.0866;       % mol (s * L mito)^(-1)             %0.0866;
X_C1  = 3e5;       % mol (s * L mito)^(-1)             %32365;
X_C3  = 5e7;       % mol (s * L mito)^(-1)             %0.79081* 40 * 1e5;
X_C4  = 0.2;       % mol (s * L mito)^(-1)             %0.00010761 * 7e2;
X_F   = 1e2;       % mol (s * L mito)^(-1)             %100; 
E_ANT = 0.065;     % mol (L mito)^(-1)    %1.5*0.006762/7.2679e-003*(0.70e-1);
E_PiC = 1e6;       % (L cell) (s * L mito)^(-1)        %3.3356e+07 * 2e-2; 
E_H   = 3e2;       % mol (s * V * L mito)^(-1)
X_CK  = 1e7;       % mol (s * L cyto)^(-1) 
X_AtC = 0.5e-3;    % mol (s * L cyto)^(-1)

pars = [X_DH; X_C1; X_C3; X_C4; X_F; 
    E_ANT; E_PiC; E_H; 
    X_CK; X_AtC; 
    ]; 

%% Initial conditions 
%{ 
This section creates a vector of initial conditions for the state
variables, the membrane potential and the concentrations of various
metabolites in the mitochondrial matrix, intermembane (IM) space, and
cytoplasm.
%}

% Membrane potential
DPsi_0    = 175*1e-3;       % V

% Matrix species 
ATP_x_0  = 0.5e-3;         % mol (L matrix water)^(-1)
ADP_x_0  = 9.5e-3;         % mol (L matrix water)^(-1)
Pi_x_0   = 0.3e-3;         % mol (L matrix water)^(-1)
NADH_x_0 = 2/3 * NAD_tot;  % mol (L matrix water)^(-1)
QH2_x_0  = Q_tot/2;        % mol (L matrix water)^(-1)

% IM species
cred_i_0  = c_tot/2;        % mol (L IM water)^(-1)

% Cytoplasmic species
ATP_c_0 = 9.95e-3;        % mol (L cyto water)^(-1)
ADP_c_0 = 0.05e-3;        % mol (L cyto water)^(-1)
Pi_c_0  = 5.0e-3;         % mol (L cyto water)^(-1)
CrP_c_0 = .3 * Cr_tot_c;  % mol (L cyto water)^(-1)

x0 = [DPsi_0; 
    ATP_x_0; ADP_x_0; Pi_x_0; NADH_x_0; QH2_x_0; 
    cred_i_0; 
    ATP_c_0; ADP_c_0; Pi_c_0; CrP_c_0; 
    ]; 

%% Solve model and plot results
%{ 
This section solves the model using stiff ODE integrator ode15s. All
outputs are converted to the appropriate units for plotting. Plots can be
printed by setting printon = 1 in the Flags section. 
%} 

[t,x] = ode15s(@model,[0 60],x0,[],pars,conc,fixedpars,flags);

% Find steady-state fluxes and polynomials 
f = zeros(14,length(t)); 
for i = 1:length(t)
    [~,f1] = model([],x(i,:),pars,conc,fixedpars,2); 
    f(:,i) = f1; 
end 

PATP_c = f(4,:)'; 
J_C4   = f(10,:)'; 

% convert to mV 
DPsi = x(:,1) * 1e3; 

% convert to mmol (L matrix water)^(-1)
sumATP_x = x(:,2) * 1e3;
sumADP_x = x(:,3) * 1e3; 
sumPi_x  = x(:,4) * 1e3; 
NADH_x   = x(:,5) * 1e3;

% convert to mmol (L IM water)^(-1)
cred_i = x(:,7) * 1e3; 

% convert to mmol (L cyto water)^(-1)
sumATP_c = x(:,8)  * 1e3; 
sumADP_c = x(:,9)  * 1e3; 
sumPi_c  = x(:,10) * 1e3; 
CrP_c    = x(:,11) * 1e3; 
Cr_c     = Cr_tot_c * 1e3 - CrP_c;

% convert to mmol (L cell)^(-1)
sumATP_x_cell = sumATP_x * (V_m * W_x); 
sumADP_x_cell = sumADP_x * (V_m * W_x); 
sumPi_x_cell  = sumPi_x  * (V_m * W_x); 

sumATP_c_cell = sumATP_c * (V_c * W_c + V_m * W_i); 
sumADP_c_cell = sumADP_c * (V_c * W_c + V_m * W_i);
sumPi_c_cell  = sumPi_x  * (V_c * W_c + V_m * W_i); 
CrP_c_cell    = CrP_c * (V_c * W_c);
Cr_c_cell     = Cr_c  * (V_c * W_c);  

% Oxygen consumption rate
R_c2t = 0.73; % cell volume (tissue volume)^(-1) 
rho   = 1.05; % tissue mass density 
MVO2  = J_C4 / 2 * V_m * R_c2t / 1000 / rho * 60 * 1e6; % umol O2 (min * g tissue)^(-1)

% [CrP]/[ATP] ratio 
ATP_c_cell = sumATP_x_cell + sumATP_c_cell; % calculate total ATP within cell
CrP_ATP = CrP_c_cell ./ ATP_c_cell; 

%% Plot figures 

linethickness = 3; 
fontsize = 20; 

% Membrane potential 
hfig1 = figure(1);
plot(t,DPsi,'b','linewidth',linethickness)
ylabel('\Delta \Psi (mV)')
xlabel('Time (s)')
set(gca,'FontSize',fontsize)

% [NADH]_x/[NAD]_tot ratio 
hfig2 = figure(2);
plot(t,NADH_x/(NAD_tot * 1e3),'b','linewidth',linethickness)
xlabel('Time (s)')
ylabel('NADH/NAD_{tot} ratio');
set(gca,'FontSize',fontsize)

% Matrix adenine species
hfig3 = figure(3);
clf 
plot(t,sumATP_x,'b',t,sumADP_x,'r',t,sumPi_x,'c','linewidth',linethickness)
title('Matrix species')
legend('ATP','ADP','Pi')
set(gca,'FontSize',fontsize)
xlabel('Time (s)')
ylabel('Concentration (mmol (L matrix water)^{-1})')

% Cytosplasm adenine species
hfig4 = figure(4);
clf
plot(t,sumATP_c,'b',t,sumADP_c,'r',t,sumPi_c,'c','linewidth',linethickness)
title('Cytoplasm species')
legend('ATP','ADP','Pi')
xlabel('Time (s)')
ylabel('Concentration (mmol (L cyto water)^{-1})')
set(gca,'FontSize',fontsize)

% Creatine species
hfig5 = figure(5);
clf
plot(t,CrP_c_cell,'b',t,Cr_c_cell,'r','linewidth',linethickness)
legend('[CrP]','[Cr]')
title('Creatine')
xlabel('Time (s)')
ylabel('Concentration (mmol (L cell)^{-1})')
set(gca,'FontSize',fontsize)

% cytochrome c
hfig6 = figure(6);
clf
plot(t,cred_i,'b','linewidth',linethickness)
title('Cytochrome c')
xlabel('Time (s)')
ylabel('Concentration (mmol (L IM water)^{-1})')
set(gca,'FontSize',fontsize)

if printon == 1
    print(hfig1,'-dpng','DPsi.png')
    print(hfig2,'-dpng','NADHratio.png')
    print(hfig3,'-dpng','Matrix_conc.png')
    print(hfig4,'-dpng','Cyto_conc.png')
    print(hfig5,'-dpng','Creatine.png')
    print(hfig6,'-dpng','cytochromec.png')
end 


%% CrP/ATP vs ATP consumption rate 
%{ 
This section calculates [CrP]/[ATP] ratio and [Pi]_c as ATP consumption
rate increases.
%} 

clear CrP_ATP

% Range of ATP hydrolysis from 0.36 to 1.2e-3 mmol (s * L cell)^(-1)
X_AtC = (0:0.1:1.52) * 1e-3 / V_c; % convert to mol (s * L cyto)^(-1)
p = pars; 

CrP_ATP = zeros(size(X_AtC)); 
Pi_c = zeros(size(X_AtC)); 
cred_i = zeros(size(X_AtC)); 
for i = 1:length(X_AtC)
    p(10) = X_AtC(i); % reassign value for X_AtC
    [t,x] = ode15s(@model,[0 60],x0,[],p,conc,fixedpars,1);
    
    sumATP_x_cell = x(end,2) * (V_m * W_x); 
    sumATP_c_cell = x(end,8) * (V_c * W_c + V_m * W_i); 
    CrP_c_cell    = x(end,11) * (V_c * W_c); 
    CrP_ATP(i)    = CrP_c_cell/(sumATP_c_cell + sumATP_x_cell); 
    
    Pi_c(i) = x(end,10); 
    cred_i(i) = x(end,7) / c_tot; 
end 

% Convert mol (s * (L cyto))^(-1) to mmol (s * (L cell))^(-1)
X_AtC = X_AtC * 1e3 * V_c; 

hfig7 = figure(7);
clf
plot(X_AtC,CrP_ATP,'b','linewidth',linethickness)
xlabel('ATP consumption rate (mmol (s * L cell)^{-1})')
ylabel('[CrP]/[ATP]')
set(gca,'FontSize',fontsize)
ylim([0 3])
xlim([0 1.5])

hfig8 = figure(8);
clf
plot(X_AtC,Pi_c*1e3,'b','linewidth',linethickness)
xlabel('ATP consumption rate (mmol (s * L cell)^{-1})')
ylabel('[Pi]_c (mol (L matrix water)^{-1})')
set(gca,'FontSize',fontsize)
xlim([0 1.5])
ylim([0 ceil(max(Pi_c*1e3))])

hfig9 = figure(9);
clf
plot(X_AtC,cred_i,'b','linewidth',linethickness)
xlabel('ATP consumption rate (mmol (s * L cell)^{-1})')
ylabel('[c^{2+}]_i / [c]_{tot}')
set(gca,'FontSize',fontsize)
xlim([0 1.5])

if printon == 1
    savefig(hfig7,'CrP_ATP.fig')
    savefig(hfig8,'Pi_c.fig')
end 