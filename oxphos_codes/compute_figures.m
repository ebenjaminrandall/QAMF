load co; % initial conditions

%function [] = compute_figures(x,co,Pi_exp_a,Pi_exp_r,...
%  d_Psi_exp_r,NADH_exp_r,d_Psi_exp_a,NADH_exp_a,JO2_exp_a);

load Pi_exp_r
load Pi_exp_a
load d_Psi_exp_r
load d_Psi_exp_a
load NADH_exp_r
load NADH_exp_a
load JO2_exp_a

% Cytochrome C data (percent reduced)
CPR_exp_r(1) = 0.092;
CPR_exp_r(2) = 0.105;
CPR_exp_a(1) = 0.084;
CPR_exp_a(2) = 0.171;

Ctot   = 2.70e-3;
NADtot = 2.97e-3;               % M; total NAD/NADH

% Resting state:
ADP_e = 0;
Mg = 5.0e-3;
for i = 1:10;
  [t,y] = ode15s(@dCdT,[0 1e3],co,[],Mg,Pi_exp_r(i)*1e-3,ADP_e);
  N = length(t); c = y(N,:)';
  J_C4 = ox_flux(c,Mg,Pi_exp_r(i)*1e-3,ADP_e);
  JO2_r(i) = (J_C4/2)*1.8957e5;
  d_Psi_r(i) = c(19);
  NADH_r(i) = c(4) / NADtot;
  pH_r(i) = -log10(c(1));
  CPR_r(i) = c(6) / Ctot;
end

% Active State:
ADP_e = 1.3e-3;
Mg = 5.0e-3;
for i = 1:10
  [t,y] = ode15s(@dCdT,[0 1e3],co,[],Mg,Pi_exp_r(i)*1e-3,ADP_e);
  N = length(t); c = y(N,:)';
  J_C4 = ox_flux(c,Mg,Pi_exp_r(i)*1e-3,ADP_e);
  JO2_a(i) = (J_C4/2)*1.8957e5;
  d_Psi_a(i) = c(19);
  NADH_a(i) = c(4) / NADtot;
  pH_a(i) = -log10(c(1));
  CPR_a(i) = c(6) / Ctot;
end

figure(1); clf; set(gca,'Fontsize',14)
plot(Pi_exp_r,d_Psi_r,'k--',Pi_exp_r,d_Psi_exp_r,'ko','linewidth',1.5,'Markerfacecolor',[1 1 1],'Markersize',8); 
xlabel('Buffer [Pi] (mM)'); ylabel('\Delta\Psi (mV)'); box on; hold on; 
plot(Pi_exp_r,d_Psi_a,'k-',Pi_exp_a,d_Psi_exp_a,'k^','linewidth',1.5,'Markerfacecolor',0.75*[1 1 1],'Markersize',8); 
hold off;

figure(2); clf; set(gca,'Fontsize',14)
plot(Pi_exp_r,NADH_r,'k--',Pi_exp_r,NADH_exp_r,'ko','linewidth',1.5,'Markerfacecolor',[1 1 1],'Markersize',8); 
xlabel('Buffer [Pi] (mM)'); ylabel('NADH (Normalized)'); box on; hold on; 
plot(Pi_exp_r,NADH_a,'k-',Pi_exp_a,NADH_exp_a,'k^','linewidth',1.5,'Markerfacecolor',0.75*[1 1 1],'Markersize',8);
set(gca,'Ytick',0.0:0.2:1.0); axis([0 10 0 1]); hold off;

figure(3); clf; set(gca,'Fontsize',14)
plot(Pi_exp_r,JO2_r,'k--','linewidth',1.5); 
xlabel('Buffer [Pi] (mM)'); ylabel('MV_{O_2} (mol O_2 min^{-1} (mol cyto A)^{-1})'); box on; hold on;
plot(Pi_exp_r,JO2_a,'k-',Pi_exp_a,JO2_exp_a,'k^','linewidth',1.5,'Markerfacecolor',0.75*[1 1 1],'Markersize',8); 
hold off

figure(4); clf; set(gca,'Fontsize',14); hold on;
plot(Pi_exp_r,CPR_r,'k--',[0 3],CPR_exp_r,'ko','linewidth',1.5,'Markerfacecolor',[1 1 1],'Markersize',8); 
xlabel('Buffer [Pi] (mM)'); ylabel('Cytochrome C Reduced Fraction'); box on; 
plot(Pi_exp_r,CPR_a,'k-',[0 3],CPR_exp_a,'k^','linewidth',1.5,'Markerfacecolor',0.75*[1 1 1],'Markersize',8); 
set(gca,'Ytick',0.0:0.05:0.3); axis([0 10 0 0.3]); hold off; 

figure(5); clf; set(gca,'Fontsize',14)
plot(Pi_exp_r,pH_r,'k--',[0 3],[7.14 7.13],'ko','linewidth',1.5,'Markerfacecolor',[1 1 1],'Markersize',8); 
xlabel('Buffer [Pi] (mM)'); ylabel('Matrix pH'); box on; hold on;
plot(Pi_exp_r,pH_a,'k-',[0 3],[7.16 7.10],'k^','linewidth',1.5,'Markerfacecolor',0.75*[1 1 1],'Markersize',8); 
set(gca,'Ytick',6.5:0.25:7.5); axis([0 10 6.5 7.5]); hold off; 

