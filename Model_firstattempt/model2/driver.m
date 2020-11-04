% Runs metabolite code
clear

%% Initial conditions

Cr_tot = 54e-3;
% Concentrations (M)
ATP_x = 1e-3; 
ADP_x = 9e-3; 
Pi_x  = 1e-3; 
ATP_e = 7.5e-3; 
ADP_e = 0.5e-3; 
Pi_e  = 0.5e-3; 
Cr    = 0.65 * Cr_tot; 

x0 = [ATP_x; ADP_x; Pi_x; ATP_e; ADP_e; Pi_e; Cr]; 

%% Solve model 

x_ATPase = 5*(1/0.6810)*0.5e-3; % ATP hydrolysis rate: M / s / (liter cytosol)

options = odeset('RelTol',1e-8','AbsTol',1e-8); 
[t,x] = ode15s(@model,[0 10],x0,options,x_ATPase); 

ATP_x = x(:,1) * 1e3; 
ADP_x = x(:,2) * 1e3; 
Pi_x  = x(:,3) * 1e3; 
ATP_e = x(:,4) * 1e3; 
ADP_e = x(:,5) * 1e3; 
Pi_e  = x(:,6) * 1e3; 
Cr    = x(:,7) * 1e3; 
CrP   = Cr_tot*1e3 - Cr; 

%% Plot 

figure(1); clf; hold on 
plot(t,ATP_x,'linewidth',2)
plot(t,ADP_x,'linewidth',2)
plot(t,Pi_x,'linewidth',2)
legend('ATP','ADP','Pi');
title('Matrix')
set(gca,'FontSize',16)

figure(2); clf; hold on 
plot(t,ATP_e,'linewidth',2)
plot(t,ADP_e,'linewidth',2)
plot(t,Pi_e,'linewidth',2)
legend('ATP','ADP','Pi');
title('Cyto')
set(gca,'FontSize',16)

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
