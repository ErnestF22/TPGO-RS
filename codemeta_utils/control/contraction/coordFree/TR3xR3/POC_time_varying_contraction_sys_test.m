% POC_time_varying_contraction_sys_test.m
% Last Edited: Oct. 21, 2020
% Test if a time-varying contracting vector field converges for any two
% curves starting at the same time

close all; clear all; clc;

% Two random initial conditions
x1 = randn(2,1);
x2 = randn(2,1);
TFinal = 10;

% Simluate
sys = @(t,x) time_varying_sys(t,x);

optsOde=odeset('MaxStep',0.01);
[t,x] = ode45(sys,[0 TFinal],[x1;x2], optsOde);

% Plot distance (2-norm) between x1 and x2 for all t
x1_t = x(:,1:2);
x2_t = x(:,3:4);
xref = [cos(t), sin(t)];
dist_x1_x2 = vecnorm(x1_t-x2_t,2,2);
dist_x1_xref = vecnorm(x1_t-xref,2,2);
dist_x2_xref = vecnorm(x2_t-xref,2,2);
figure
plot(t,dist_x1_x2,'LineWidth',5);
hold on
plot(t,dist_x1_xref,'LineWidth',5);
plot(t,dist_x2_xref,'LineWidth',5);
legend('d(x1,x2)','d(x1,xref)','d(x2,xref)');

function dx = time_varying_sys(t,x)
% Define the system f(x) = x_{ref}(t) - x for same systems with different
% initial conditions where
% x_{ref} = [cos(t);sin(t)];
% x = [ [x1(1);x1(2)]; [x2(1);x2(2)] ]
dx = [cos(t) - x(1);...
    sin(t) - x(2);...
    cos(t) - x(3);...
    sin(t) - x(4)];
end