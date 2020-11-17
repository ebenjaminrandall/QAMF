function f = dXdT_electrode(t,Jel,tsim,Jox)

tau = 4; % seconds

f = (interp1(tsim,Jox,t)- Jel)/tau;

% The Jox and tsim are the time and flux vector from the model simulation.
% 
% The basic model is:
%
% d(J_electrode)/dt = (Jxo - J_electrode)/tau