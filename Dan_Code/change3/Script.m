clear 


%% Parameters defining metabolite pools

% Volume fractions 
V_c = 0.6601;       % cytosol volume fraction       % L cyto (L cell)^(-1)
V_m = 0.2882;       % mitochondrial volume fraction % L mito (L cell)^(-1)
R_m2c = V_m / V_c;  % mito to cyto volume ratio     % L mito (L cyto)^(-1)
W_c = 0.8425;       % cytosol water space           % L cyto water (L cyto)^(-1) 
W_m = 0.7238;       % mitochondrial water space     % L mito water (L mito)^(-1)
W_x = 0.9*W_m;      % matrix water space            % L matrix water (L mito)^(-1)
W_i = 0.1*W_m;      % intermembrane water space     % L IM water (L mito)^(-1)

% Total pool concentrations 
NAD_tot = 2.97e-3;  % NAD+ and NADH conc            % mol (L matrix water)^(-1)
Q_tot   = 1.35e-3;  % Q and QH2 conc                % mol (L matrix water)^(-1) 
c_tot   = 2.7e-3;   % cytochrome c ox and red conc  % mol (L IM water)^(-1)
Cr_tot  = 30e-3;    % creatine and creatine phosphate conc % mol (L cell)^(-1)

Crtot_c = Cr_tot / (V_c * W_c); % convert to mol (L cyto water)^(-1)

% a  = 0.082;
% p  = a*(29.78/8.62); 
% c  = 0.39837; 
% A0 = 8.62  + 20*a;
% P0 = 29.78 + 20*p;
% C0 = 35.04 + 20*c;

%% Adjustable parameters 

x_DH  = 0.0866;
x_C1  = 32365;
x_C3  = 0.79081;
x_C4  = 0.00010761;
x_F1  = 8120.4;
x_ANT = 0.006762;
x_PI1 = 3.3356e+07;
x_CK  = 1e7;
x_AtC = 0.50e-3;

pars = [x_DH; x_C1; x_C3; x_C4; 
    x_F1; x_ANT; x_PI1; 
    x_CK; x_AtC;
    ]; 

%% Initial conditions

% (i) Matrix species and dPsi
idPsi        = 1;
iATP_x       = 2;
iADP_x       = 3;
iPI_x        = 4;
iNADH_x      = 5;
iQH2_x       = 6;
%  (ii) IM space species
iCred_i      = 7;
%iATP_i       = 8;
%iADP_i       = 9;
%iAMP_i       = 10;
%iPI_i        = 11;
%  (iii) Cytoplasmic species
iATP_c       = 8; %12;
iADP_c       = 9; %13;
iPI_c        = 10; %14;
iPCr_c       = 11; %15;
%iAMP_c       = 16;
%iCr_c        = 17;

x0(idPsi)    = 175;
x0(iATP_x)   = 0.5e-3;
x0(iADP_x)   = 9.5e-3;
x0(iPI_x)    = 0.3e-3;
x0(iNADH_x)  = 2*NAD_tot/3;
x0(iQH2_x)   = Q_tot/2;
x0(iCred_i)  = c_tot/2;
%x0(iATP_i)   = 10e-3;
%x0(iADP_i)   = 0;
%x0(iAMP_i)   = 0;
%x0(iPI_i)    = 1.0e-3;
x0(iATP_c)   = 10e-3;
x0(iADP_c)   = 0;
x0(iPI_c)    = 1.0e-3;
x0(iPCr_c)   = .3 * Crtot_c; %22e-3;
%x0(iAMP_c)   = 0;
%x0(iCr_c)    = 18e-3;

%% Solve model 

[t,x] = ode15s(@Cell_dXdT,[0 60],x0,[],pars);

%% Plot figures 

figure(1)
plot(t,x(:,1))
title('\Delta\Psi (mV)');
set(gca,'FontSize',20)

figure(2)
plot(t,x(:,iPCr_c)./x(:,iATP_c) * (V_c*W_c))
title('CrP/ATP');
set(gca,'FontSize',20)

figure(3)
plot(t,x(:,iPI_c)*1e3,t,x(:,iPI_x)*1e3) 
title('[Pi]')
legend('cyto','matrix');
set(gca,'FontSize',20)

figure(4)
plot(t,x(:,iNADH_x)/NAD_tot)
title('NADH');
set(gca,'FontSize',20)

figure(5)
clf 
plot(t,x(:,iATP_c)*1e3,'b',t,x(:,iATP_x)*1e3,'r')
title('[ATP]')
legend('cyto','matrix')
set(gca,'FontSize',20)

figure(6)
clf
plot(t,x(:,iADP_c)*1e3,'b',t,x(:,iADP_x)*1e3,'r')
title('[ADP]')
legend('cyto','matrix')
set(gca,'FontSize',20)

figure(7)
clf
plot(t,x(:,iPCr_c)*1e3 * (V_c*W_c),'b',t,(Crtot_c - x(:,iPCr_c))*1e3 * (V_c*W_c),'r')
legend('[PCr]','[Cr]')
title('Creatine')
set(gca,'FontSize',20)

figure(8)
clf
plot(t,x(:,iCred_i)*1e3,'b')
title('Cytochrome c')
set(gca,'FontSize',20)


%% CrP/ATP vs ATP consumption rate 

% Range of ATP hydr02olysis from 0.36 to 1.2e-3 mmol  / sec / (l cell)
x_AtC = (0.36:0.1:1.52)*1e-3;
p = pars; 

CrP_ATP = zeros(size(x_AtC)); 
for i = 1:length(x_AtC)
    p(end) = x_AtC(i); 
    [t,x] = ode15s(@Cell_dXdT,[0 60],x0,[],p);
    CrP_ATP(i) = x(end,iPCr_c)/x(end,iATP_c) * (V_c*W_c); 
end 

figure(9)
clf
plot(x_AtC,CrP_ATP,'b')
xlabel('ATP consumption rate')
ylabel('[CrP]/[ATP]')
set(gca,'FontSize',20)




    




