% Test code to compute the gradient of (12) in
% contractionQuadrotorControl.pdf. First we double check that the gradient
% and hessian for a stationary target is correct, then check it for a
% moving target.
close all; clear all; clc;

% Default parameters
ERR_TOL = 1e-6;
t = 10*rand; % Choose a random to eval
x0 = 10*randn(3,1); % Random start point
v0 = 10*randn(3,1); % Random constant velocity
%% Assume xd = 0
% Define a curve(line) starting at x0
x = @(t) x0 + v0*t;
% Define the cost function
f = @(x) 2*(sqrt(1+(norm(x)^2/2))-1);
% Define what we think the analytical grad is
gradf = @(x) x/sqrt(1+norm(x)^2/2);
hessf = @(x) (sqrt(1+norm(x)^2/2)*eye(3)-x*x'/(2*sqrt(1+norm(x)^2/2)))/(sqrt(1+norm(x)^2/2))^2;

% Test Gradient
fdot_error = funApproxDer(@(t) f(x(t)),t) - gradf(x(t))'*funApproxDer(x,t);
if abs(fdot_error) > ERR_TOL
    error('Analytical gradient is incorrect');
end

% Test Hessian
hessf_error = funApproxDer(@(t) gradf(x(t)),t) - hessf(x(t))*funApproxDer(x,t);
if any(abs(hessf_error)>ERR_TOL)
    error('Analytical hessian is incorrect');
end

%% Assume xd = f(t), this is the trajectory following vector field
xd0 = 10*randn(3,1); % Initial desired pos
vd0 = 10*randn(3,1); % Constant desired vel
xd = @(t) xd0 + vd0*t;
% Define the new cost function
f2 = @(x,xd) 2*(sqrt(1+(norm(x-xd)^2/2))-1);
% Define what we think the new analytical grad is
gradf2 = @(x,xd) (x-xd)/sqrt(1+norm(x-xd)^2/2);
hessf2 = @(x,xd) (sqrt(1+norm(x-xd)^2/2)*eye(3)-(x-xd)*(x-xd)'/(2*sqrt(1+norm(x-xd)^2/2)))/(sqrt(1+norm(x-xd)^2/2))^2;

% Test Gradient
fdot_error2 = funApproxDer(@(t) f2(x(t),xd(t)),t) - gradf2(x(t),xd(t))'*(funApproxDer(x,t)-funApproxDer(xd,t));
if abs(fdot_error2) > ERR_TOL
    error('Analytical gradient2 is incorrect');
end

% Test Hessian
hessf_error2 = funApproxDer(@(t) gradf2(x(t),xd(t)),t) - hessf2(x(t),xd(t))*(funApproxDer(x,t)-funApproxDer(xd,t));
if any(abs(hessf_error2)>ERR_TOL)
    error('Analytical hessian2 is incorrect');
end