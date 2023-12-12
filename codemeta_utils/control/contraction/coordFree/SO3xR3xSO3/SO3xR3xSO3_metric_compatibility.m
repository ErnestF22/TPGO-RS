function [Ft, dFt, appder] = SO3xR3xSO3_metric_compatibility(R1,dR1,p,dp,R2,dR2,X,t,varargin)
% Check if d/dt<X,Y> = <D_{Z_dot}X,Y> + <X,D_{Z_dot}Y> evaluated at
% (R1,p,R2,t). Assume Y = Z_dot.
% INPUTS:
%   R1(t) := Curve on SO(3)
%   dR1(R1) := Vector field along R1(t)
%   p(t) := Curve on R3 in skew symm form
%   dp(p) := Vector field along p(t)
%   R2(t) := Curve on SO(3)
%   dR2(R2) := Vector field along R2(t)
%   X(R1,p,R2) := Vector field defined on SO3xR3xSO3
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
Z = @(t) [R1(t);p(t);R2(t)]; % Curve on product manifold
Z_dot = @(R1,p,R2) [dR1(R1);dp(p);dR2(R2)];
Y = Z_dot;

% Left Side d/dt <X,Y>
Left = @(t) SO3xR3xSO3_metric(Z(t),X(R1(t),p(t),R2(t)),...
    Y(R1(t),p(t),R2(t)));

% Right side SO3 for R1
X_SO3_1 = @(R1) extractComp(X(R1,p(t),R2(t)),1,3,1,3);
Y_SO3_1 = @(R1) extractComp(Y(R1,p(t),R2(t)),1,3,1,3);
D_SO3_1_1 = @(R1) rot_metric(R1, rot_covar(dR1,X_SO3_1), Y_SO3_1);
D_SO3_1_2 = @(R1) rot_metric(R1, X_SO3_1, rot_covar(dR1,Y_SO3_1));
Right_SO3_1 = @(t) D_SO3_1_1(R1(t)) + D_SO3_1_2(R1(t));
% Right side R3
X_R3 = @(p) extractComp(X(R1(t),p,R2(t)),4,6,1,3);
X_R3_dot = @(t) funApproxDer(@(t0) X_R3(p(t0)),t);
Y_R3 = @(p) extractComp(Y(R1(t),p,R2(t)),4,6,1,3);
Y_R3_dot = @(t) funApproxDer(@(t0) Y_R3(p(t0)),t);
Right_R = @(t) 1/2*trace(X_R3_dot(t)'*Y_R3(p(t)) + X_R3(p(t))'*Y_R3_dot(t));
% Right side SO3 for R2
X_SO3_2 = @(R2) extractComp(X(R1(t),p(t),R2),7,9,1,3);
Y_SO3_2 = @(R2) extractComp(Y(R1(t),p(t),R2),7,9,1,3);
D_SO3_2_1 = @(R2) rot_metric(R2, rot_covar(dR2,X_SO3_2), Y_SO3_2);
D_SO3_2_2 = @(R2) rot_metric(R2, X_SO3_2, rot_covar(dR2,Y_SO3_2));
Right_SO3_2 = @(t) D_SO3_2_1(R2(t)) + D_SO3_2_2(R2(t));
% Resulting right side
Right = @(t) Right_SO3_1(t) + Right_R(t) + Right_SO3_2(t);

% Perform check
[Ft, dFt, appder] = funCheckDer(Left, Right, t, 'nodisplay');

if bCoupledDynamics
    D_Extra_SO3_1 = 1/2*trace(dR1(R1(t))'*R1(t)*dp(p(t)));
    D_Extra_R = 1/2*trace(dp(p(t))'*...
        ( funApproxDer( @(t0) R1(t0)'*rot_log(R1(t0),R2(t0)),t) ) );
    dFt = dFt + D_Extra_SO3_1 + D_Extra_R;
end
end

