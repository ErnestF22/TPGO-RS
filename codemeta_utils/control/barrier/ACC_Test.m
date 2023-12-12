% Script to test the Adaptive Cruise Control (ACC) from 
% 'Control Barrier Function based Quadratic Programs with Application to
% Adaptive Cruise Control' 2014.

close all; clear all; clc;

x0 = [900; 20; 100; 1000; 13.89];
TFinal = 20;
control = @(t, x) controller(x);
closedLoop = @(t, x) model(x, control(t,x));
optsOde=odeset('OutputFcn',@odeplot,'MaxStep',0.01);
% optsOde=odeset('MaxStep',0.01);
subplot(2,2,1);
[t,x]=ode45(closedLoop,[0 TFinal],x0, optsOde);
title('States Plot');

subplot(2,2,2)
plot(t,x(:,2));
hold on
plot(t,24*ones(size(t)),'r--');
title('Velocity Plot');
subplot(2,2,3)
h = x(:,3)-1.8*x(:,2);
plot(t,h);
title('h Plot (h>0)');


function [dx] = model(x, u)
% INPUTS:
%   x = state vector
%       x(1:3) = [pos, vel, distance of x(5) - x(1)] of controlled car
%       x(4:5) = [pos, vel] of constant car (leader)
%   u = control input into controlled vehicle
% OUTPUTS:
%   dx = new states

% Parameters
m = 1650; % kg
f0 = 0.1; % Newton
f1 = 5; % Newton*s/m
f2 = 0.25; % Newton*s^2/m
vd = 24; % m/s
v0 = 13.89; %m/s
eps = 10; % exp convergence rate for stability
gamma = 1; % max growth rate of Bdot leq gamma
p_sc = 1e-5; % relaxation term compensation
Fr = f0 + f1*x(2) + f2*x(2)^2; % Aerodyn. drag

% State
dx(1) = x(2); % vel.
dx(2) = -1/m*Fr + 1/m*u(1); % acc.
dx(3) = x(5) - x(2); % difference in vel.
dx(4) = v0;
dx(5) = 0;

dx = dx';
end

function [u] = controller(x)
% Parameters
m = 1650; % kg
f0 = 0.1; % Newton
f1 = 5; % Newton*s/m
f2 = 0.25; % Newton*s^2/m
vd = 24; % m/s
v0 = 13.89; %m/s
eps = 10; % exp convergence rate for stability
gamma = 1; % max growth rate of Bdot leq gamma
p_sc = 1e-1; % relaxation term compensation
Fr = f0 + f1*x(2) + f2*x(2)^2; % Aerodyn. drag

% CLF
phi0 = -2*(x(2)-vd)*Fr/m + eps*(x(2)-vd)^2;
phi1 = 2*(x(2)-vd)/m;

% CBF
z = x(4) - x(1);
h = x(3) - 1.8*x(2); % set function
Bf = -log(h/(1+h)); % Barrier value 
denum = m*(1-1.8*x(2)+z)*(-1.8*x(2)+z);
LfB = -(1.8*Fr+m*(x(5)-x(2)))/denum;
LgB = 1.8/denum;

% Ouadratic Problem Solver
A_clf = [phi1, -1];
% A_clf = phi1;
b_clf = -phi0;
A_cbf = [LgB, 0];
b_cbf = -LfB + gamma/Bf;
H_acc = 2*[1/m^2 0;0 p_sc];
F_acc = -2*[Fr/m^2;0];

% CVXR implementation
% cvx_begin quiet
% %     variable u(2) % [u; delta_sc]
%     variable u1
%     u = [u1; 0];
%     minimize( 1/2*u'*H_acc*u + F_acc'*u )
%     subject to
%         A_clf*u <= b_clf;
% %         A_cbf*u <= b_cbf;
%         
% cvx_end

% quadprog Implementation
options = optimoptions('quadprog',...
    'Algorithm','interior-point-convex','Display','off');
[u,fval,exitflag,output,lambda] = ...
   quadprog(H_acc,F_acc,[A_clf; A_cbf], [b_clf; b_cbf] ,...
   [],[],[],2e5*ones(2,1),[],options);
% [u,fval,exitflag,output,lambda] = ...
%    quadprog(H_acc,F_acc,A_clf, b_clf ,...
%    [],[],[],2e5*ones(2,1),[],options);
% % Display u
% u
% Bdot = LfB + LgB*u(1) - gamma/Bf
% Bf
% if size(u,1) == 1
%     u = u';
% end

end