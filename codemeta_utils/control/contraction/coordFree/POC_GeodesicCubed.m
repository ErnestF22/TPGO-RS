% Test for eigenvalues of D^2(I,R) and D^3(I,R) on SO(3)

close all; clear all; clc;

w = cnormalize(randn(3,1)); % random vector
R = @(t) rot_expVec(eye(3),t*w); %random geodesic
t = linspace(1e-5,pi-1e-5); % define the distance
symm = @(A) (A+A')/2; % function to find symmetric part

theta = @(t) norm(rot3_log(R(t))); % theta is the norm of w
u = @(t) rot_vee(eye(3),rot3_log(R(t))/theta(t)); % normalize logR

HessD2 = @(t) rot3_logDiffMat(eye(3),R(t)); % Hessian of D^2
HessD3 = @(t) rot_vee(eye(3),rot3_log(R(t)))*u(t)' + theta(t)*HessD2(t); % Hessian of D^3...?

figure
subplot(1,2,1)
funPlot(@(t) sort(eig(symm(HessD2(t)))), t);
title('Distance Squared')
xlabel('distance')
ylabel('eVal')
subplot(1,2,2)
funPlot(@(t) sort(eig(symm(HessD3(t)))), t);
title('Distance Cubed')
xlabel('distance')
ylabel('eVal')