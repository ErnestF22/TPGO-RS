function [ Ft, dFt, appder ] = rot_metric_CompatibilityTest(R, dR, X, Y )
% Check if d/dt <X,Y> = <D_{dR}X, Y> + <X, D_{dR}Y>. This test will help
% determine whether or not the covariant derivative is properly
% defined/implemented on SO(n).
% INPUTS:
%   R(t) := A geodesic on SO(n)
%   dR(R) := The derivative of R as a function handle of paramter R
%   X(R), Y(R) := Vector fields on SO(n) given as function handles of R
% OUTPUTS:
%   Ft := numerical evaluation of <X,Y>
%   dFt := analytical evaluation of <D_dR{X},Y> + <X, D_dGamma{Y}>
%   appder := numerical approximation of d/dt <X,Y>

% Left Term
Left_Integral = @(R) rot_metric(R, X, Y);

% Right Terms
Right_1 = @(R) rot_metric(R, rot_covar(dR,X), Y);
Right_2 = @(R) rot_metric(R, X, rot_covar(dR,Y));
Right = @(R) Right_1(R) + Right_2(R);

% Perform check
[Ft, dFt, appder] = funCheckDer(@(t) Left_Integral(R(t)), @(t) Right(R(t)), ...
    linspace(0,1));

end

