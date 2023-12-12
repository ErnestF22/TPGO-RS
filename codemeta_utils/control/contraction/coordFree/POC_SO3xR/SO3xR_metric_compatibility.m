function [Ft, dFt, appder] = SO3xR_metric_compatibility(R,dR,p,dp,X,t,varargin)
% Check if d/dt<X,Y> = <D_{Z_dot}X,Y> + <X,D_{Z_dot}Y> evaluated at
% (R,p,t). Assume Y = Z_dot.
% INPUTS:
%   R(t) := Curve on SO(3)
%   dR(R) := Vector field along R(t)
%   p(t) := Curve on R
%   dp(p) := Vector field along p(t)
%   X(R,p) := Vector field defined on SO3xR
%   t := time of evaluation
% OUTPUTS:
%   Ft := numerical evaluation of <X,Y>
%   dFt := analytical evaluation of <D_{Z_dot}X, Y> + <X, D_{Z_dot}Y>
%   appder := numerical approximation of d/dt <X,Y>

%optional parameters
bCoupledDynamics = false;
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'coupled'
            bCoupledDynamics = true;
        otherwise
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

% Define useful parameters
Z = @(t) [R(t);p(t)*ones(1,3)]; % Curve on SO3xR
Z_dot = @(R,p) [dR(R);dp(p)*ones(1,3)]; % d/dt(Z)
Y = Z_dot;

% Left side d/dt <X,Y>
Left = @(t) SO3xR_metric(Z(t),X(R(t),p(t)),Y(R(t),p(t)));

% Right side SO3 comp
X_SO3 = @(R) extractComp(X(R,p(t)),1,3,1,3);
Y_SO3 = @(R) extractComp(Y(R,p(t)),1,3,1,3);
D_SO3_1 = @(R) rot_metric(R, rot_covar(dR,X_SO3), Y_SO3);
D_SO3_2 = @(R) rot_metric(R, X_SO3, rot_covar(dR,Y_SO3));
Right_SO3 = @(t) D_SO3_1(R(t)) + D_SO3_2(R(t));
% Right side R comp
X_p = @(p) extractComp(X(R(t),p),4,4,1,1);
X_p_dot = @(t) funApproxDer(@(t0) X_p(p(t0)),t);
Y_p = @(p) extractComp(Y(R(t),p),4,4,1,1);
Y_p_dot = @(t) funApproxDer(@(t0) Y_p(p(t0)),t);
Right_R = @(t) 1/2*( X_p_dot(t)*Y_p(p(t)) + X_p(p(t))*Y_p_dot(t) );
% Resulting Right Side
Right = @(t) Right_SO3(t) + Right_R(t);
% Perform check
[Ft, dFt, appder] = funCheckDer(Left, Right, t, 'nodisplay');

% Additional Term for coupled dynamics
if bCoupledDynamics
    D_SO3_Extra = @(R,p) rot_metric(R,dp(p)*rot_log(R,eye(3)),Y_SO3(R));
    dFt = dFt + D_SO3_Extra(R(t),p(t));
end
end

