clear 


% parameters defining metabolite pools
a  = 0.082;
p  = a*(29.78/8.62); 
c  = 0.39837; 
A0 = 8.62  + 20*a;
P0 = 29.78 + 20*p;
C0 = 35.04 + 20*c;
Ctot   = 2.70e-3;               % M; total cytoC
Qtot   = 1.35e-3;               % M; total Q + QH2
NADtot = 2.97e-3;               % M; total NAD+NADH

% Initial conditions
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

x0(idPsi)    = 175;
x0(iATP_x)   = 0.5e-3;
x0(iADP_x)   = 9.5e-3;
x0(iPI_x)    = 0.3e-3;
x0(iNADH_x)  = 2*NADtot/3;
x0(iQH2_x)   = Qtot/2;
x0(iCred_i)  = Ctot/2;
x0(iATP_i)   = 10e-3;
x0(iADP_i)   = 0;
x0(iAMP_i)   = 0;
x0(iPI_i)    = 1.0e-3;
x0(iATP_c)   = 10e-3;
x0(iADP_c)   = 0;
x0(iPI_c)    = 1.0e-3;
x0(iPCr_c)   = 22e-3;
x0(iAMP_c)   = 0;
x0(iCr_c)    = 18e-3;


% Range of ATP hydrolysis from 0.36 to 1.2e-3 mmol  / sec / (l cell)
% x_AtC = (0.36:0.02:1.52)*1e-3;
x_AtC = 0.50e-3;


[t,x] = ode15s(@Cell_dXdT,[0 60],x0,[],x_AtC);


figure(1)
plot(t,x(:,1))
title('\Delta\Psi (mV)');

figure(2)
plot(t,x(:,iPCr_c)./x(:,iATP_c))
title('CrP/ATP');

figure(3)
plot(t,x(:,iPI_c)*1e3,t,x(:,iPI_x)*1e3) 
title('[Pi]')
legend('cyto','matrix');

figure(4)
plot(t,x(:,iNADH_x)/NADtot)
title('NADH');

figure(5)
clf 
plot(t,x(:,iATP_c)*1e3,'b',t,x(:,iATP_x)*1e3,'r')
title('[ATP]')
legend('cyto','matrix')

figure(6)
clf
plot(t,x(:,iADP_c)*1e3,'b',t,x(:,iADP_x)*1e3,'r')
title('[ADP]')
legend('cyto','matrix')

figure(7)
clf
plot(t,x(:,iPCr_c)*1e3,'b',t,x(:,iCr_c)*1e3,'r')
legend('[PCr]','[Cr]')
title('Creatine')

figure(8)
clf
plot(t,x(:,iCred_i)*1e3,'b')
title('Cytochrome c')




