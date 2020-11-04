% Runs metabolite code
clear

Cr_tot   = 40   * 1e-3;
Q_tot    = 1.35 * 1e-3; 
NAD_tot  = 2.97 * 1e-3; 
cytc_tot = 2.7  * 1e-3; 

%% Initial conditions

% Concentrations (M)
ATP_x_0    = 1e-3; 
ADP_x_0    = 10 * 1e-3; 
Pi_x_0     = 1e-3; 
ATP_e_0    = 7.5 * 1e-3; 
ADP_e_0    = 0.5 * 1e-3; 
Pi_e_0     = 0.3 * 1e-3; 
Cr_0       = 0.65 * Cr_tot; 
QH2_0      = 0.5 * Q_tot; 
NADH_x_0   = 0.5 * NAD_tot; 
cytC_red_0 = 0.5 * cytc_tot; 

% Membrance potential (V) 
DPsi_0 = 175 * 1e-3; 

x0 = [ATP_x_0; ADP_x_0; Pi_x_0; ATP_e_0; ADP_e_0; Pi_e_0; ...
    Cr_0; QH2_0; NADH_x_0; cytC_red_0; DPsi_0]; 

%% Solve model 

x_ATPase = 0.547*1e-3;%5*(1/0.6810)*0.5e-3; % ATP hydrolysis rate: M / s / (liter cytosol)

tend = 60; 

options = odeset('RelTol',1e-8','AbsTol',1e-8); 
[t,x] = ode15s(@model,[0 tend],x0,options,x_ATPase); 

ATP_x   = x(:,1)  * 1e3; 
ADP_x   = x(:,2)  * 1e3; 
Pi_x    = x(:,3)  * 1e3; 
ATP_e   = x(:,4)  * 1e3; 
ADP_e   = x(:,5)  * 1e3; 
Pi_e    = x(:,6)  * 1e3; 
Cr      = x(:,7)  * 1e3; 
QH2     = x(:,8)  * 1e3; 
NADH_x  = x(:,9)  * 1e3; 
cytcred = x(:,10) * 1e3; 
DPsi    = x(:,11) * 1e3; 

CrP    = Cr_tot   * 1e3 - Cr; 
Q      = Q_tot    * 1e3 - QH2; 
NAD_x  = NAD_tot  * 1e3 - NADH_x; 
cytcox = cytc_tot * 1e3 - cytcred; 

%% Plot 

figure(1)
clf
hold on 
plot(t,ATP_x,'linewidth',2)
plot(t,ADP_x,'linewidth',2)
plot(t,Pi_x,'linewidth',2)
legend('ATP','ADP','Pi');
title('Matrix')
xlabel('Time (s)')
ylabel('Concentration (mM)')
set(gca,'FontSize',16)

figure(2)
clf
hold on 
plot(t,ATP_e,'linewidth',2)
plot(t,ADP_e,'linewidth',2)
plot(t,Pi_e,'linewidth',2)
legend('ATP','ADP','Pi');
title('Cytosol')
xlabel('Time (s)')
ylabel('Concentration (mM)')
set(gca,'FontSize',16)

figure(3)
clf
hold on 
plot(t,Cr,'linewidth',2)
plot(t,CrP,'linewidth',2)
legend('Cr','CrP')
xlabel('Time (s)')
ylabel('Concentration (mM)')
set(gca,'FontSize',16)

figure(4)
clf
hold on 
plot(t,QH2,'linewidth',2)
plot(t,Q,'linewidth',2)
legend('QH2','Q')
xlabel('Time (s)')
ylabel('Concentration (mM)')
set(gca,'FontSize',16)

figure(5)
clf
hold on 
plot(t,NADH_x,'linewidth',2)
plot(t,NAD_x,'linewidth',2)
legend('NADH','NAD')
xlabel('Time (s)')
ylabel('Concentration (mM)')
set(gca,'FontSize',16)

figure(6)
clf
hold on 
plot(t,cytcred,'linewidth',2)
plot(t,cytcox,'linewidth',2)
legend('cytcred','cytcox')
xlabel('Time (s)')
ylabel('Concentration (mM)')
set(gca,'FontSize',16)

figure(7)
clf
hold on 
plot(t,DPsi,'linewidth',2)
title('Membrane Potential')
xlabel('Time (s)')
ylabel('Membrane Potential (mV)')
set(gca,'FontSize',16)

return 


%% Solve model 

clear ATP_x ADP_x Pi_x ATP_e ADP_e Pi_e Cr CrP

x_ATPase = (0:1:10).*(1/0.6810)*0.5e-3; % ATP hydrolysis rate: M / s / (liter cytosol)

options = odeset('RelTol',1e-8','AbsTol',1e-8); 

for i = 1:length(x_ATPase)
    
  [t,x] = ode15s(@model,[0 60],x0,options,x_ATPase(i)); 

  ATP_x(i) = x(end,1) * 1e3; 
  ADP_x(i) = x(end,2) * 1e3; 
  Pi_x(i)  = x(end,3) * 1e3; 
  ATP_e(i) = x(end,4) * 1e3; 
  ADP_e(i) = x(end,5) * 1e3; 
  Pi_e(i)  = x(end,6) * 1e3; 
  Cr(i)    = x(end,7) * 1e3; 
  CrP(i)   = Cr_tot*1e3 - Cr(i); 

end

figure(3)
clf
plot(x_ATPase * 1e3, Pi_e)
title('Pi')

figure(4)
clf
plot(x_ATPase * 1e3, CrP./ATP_e)
title('CrP/ATP')
