% Test Safe RRT implementation15
close all; clear all; clc;
% Define Simulation conditions
% % RANDOM CONDITION
% x0 = randn(3,1);
% xd = randn(3,1);
% FIXED conditions
x0 = -[1;1;0];
xd = zeros(3,1);
TFinal = 100;

% Define System
control = @(t,x) unicycle_CLFCBF(x,xd);
closedLoop = @(t,x) unicycle(x, control(t,x));
optsOde=odeset('OutputFcn',@odeplot,'MaxStep',0.01,'Events',@(t,x) systemConverged(t,x,xd));
figure;
[t,x]=ode45(closedLoop,[0 TFinal],x0, optsOde);
u = zeros(length(t),2);
for i = 1:length(t)
    u(i,:) = control(t(i),x(i,:)')';
end
figure;
plot(t,u);
legend('v','w');
title('Control Effort')

% Plot trajectory
animateUnicycle(t,x,'stepsize',10,'target',xd);
% plotUnicycle(t,x,'stepsize',1);
plotCLF(x,t,xd);

% Determine when to stop ode45
function [CLF,isterminal,direction] = systemConverged(t,x,xd)
CLF = compute_CLF(x,xd);
if (CLF < 1e-6)
    CLF = 0;
end
isterminal = 1;  % Halt integration 
direction = 0;   % The zero can be approached from either direction
end