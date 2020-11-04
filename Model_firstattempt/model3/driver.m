% Runs metabolite code
clear


Cr_tot = 54e-3;
Q_tot    = 1.35e-3; 
NAD_tot  = 2.97e-3; 
cytC_tot = 2.7e-3; 

%% Initial conditions

% Concentrations (M)
ATP_x_0   = 1e-3; 
ADP_x_0   = 9e-3; 
Pi_x_0    = 1e-3; 
ATP_e_0   = 7.5e-3; 
ADP_e_0   = 0.5e-3; 
Pi_e_0    = 0.5e-3; 
Cr_0      = 0.65 * Cr_tot; 
QH2_0     = 0.5 * Q_tot; 
NADH_x_0  = 0.5 * NAD_tot; 
cytCred_0 = 0.5 * cytC_tot;  

x0 = [ATP_x_0; ADP_x_0; Pi_x_0; ATP_e_0; ADP_e_0; Pi_e_0; ...
    Cr_0; QH2_0; NADH_x_0; cytCred_0]; 

%% Solve model 

x_ATPase = 5*(1/0.6810)*0.5e-3; % ATP hydrolysis rate: M / s / (liter cytosol)

options = odeset('RelTol',1e-8','AbsTol',1e-8); 
[t,x] = ode15s(@model,[0 10],x0,options,x_ATPase); 

ATP_x_0 = x(:,1) * 1e3; 
ADP_x_0 = x(:,2) * 1e3; 
Pi_x_0  = x(:,3) * 1e3; 
ATP_e_0 = x(:,4) * 1e3; 
ADP_e_0 = x(:,5) * 1e3; 
Pi_e_0  = x(:,6) * 1e3; 
Cr_0    = x(:,7) * 1e3; 
CrP   = Cr_tot*1e3 - Cr_0; 

%% Plot 

figure(1); clf; hold on 
plot(t,ATP_x_0,'linewidth',2)
plot(t,ADP_x_0,'linewidth',2)
plot(t,Pi_x_0,'linewidth',2)
legend('ATP','ADP','Pi');
title('Matrix')
set(gca,'FontSize',16)

figure(2); clf; hold on 
plot(t,ATP_e_0,'linewidth',2)
plot(t,ADP_e_0,'linewidth',2)
plot(t,Pi_e_0,'linewidth',2)
legend('ATP','ADP','Pi');
title('Cyto')
set(gca,'FontSize',16)


return 


%% Solve model 

clear ATP_x ADP_x Pi_x ATP_e ADP_e Pi_e Cr CrP

x_ATPase = (0:1:10).*(1/0.6810)*0.5e-3; % ATP hydrolysis rate: M / s / (liter cytosol)

options = odeset('RelTol',1e-8','AbsTol',1e-8); 

for i = 1:length(x_ATPase)
    
  [t,x] = ode15s(@model,[0 60],x0,options,x_ATPase(i)); 

  ATP_x_0(i) = x(end,1) * 1e3; 
  ADP_x_0(i) = x(end,2) * 1e3; 
  Pi_x_0(i)  = x(end,3) * 1e3; 
  ATP_e_0(i) = x(end,4) * 1e3; 
  ADP_e_0(i) = x(end,5) * 1e3; 
  Pi_e_0(i)  = x(end,6) * 1e3; 
  Cr_0(i)    = x(end,7) * 1e3; 
  CrP(i)   = Cr_tot*1e3 - Cr_0(i); 

end

figure(3)
clf
plot(x_ATPase * 1e3, Pi_e_0)
title('Pi')

figure(4)
clf
plot(x_ATPase * 1e3, CrP./ATP_e_0)
title('CrP/ATP')
