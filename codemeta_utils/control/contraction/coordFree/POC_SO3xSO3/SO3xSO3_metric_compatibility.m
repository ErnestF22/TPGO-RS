function [Ft, dFt, appder] = SO3xSO3_metric_compatibility(R1,dR1,R2,dR2,X,t,varargin)
% Check if d/dt<X,Y> = <D_{Z_dot}X,Y> + <X,D_{Z_dot}Y> evaluated at
% (R1,R2,t). Assume Y = Z_dot.
% INPUTS:
%   R1(t) := Curve on SO(3)
%   dR1(R) := Vector field along R1(t)
%   R2(t) := Curve on SO(3)
%   dR2(R) := Vector field along R2(t)
%   X(R1,R2) := Vector field defined on SO3xSO3
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
Z = @(t) [R1(t);R2(t)]; % Curve on SO3xSO3
Z_dot = @(R1,R2) [dR1(R1);dR2(R2)]; % d/dt(Z)
Y = Z_dot;

% Left side d/dt <X,Y>
Left = @(t) SO3xSO3_metric(Z(t), X(R1(t),R2(t)), Y(R1(t),R2(t)) );


% Right side SO3 for R1
X_SO3_1 = @(R1) extractComp(X(R1,R2(t)),1,3,1,3);
Y_SO3_1 = @(R1) extractComp(Y(R1,R2(t)),1,3,1,3);
D_SO3_1_1 = @(R1) rot_metric(R1, rot_covar(dR1,X_SO3_1), Y_SO3_1);
D_SO3_1_2 = @(R1) rot_metric(R1, X_SO3_1, rot_covar(dR1,Y_SO3_1));
Right_SO3_1 = @(t) D_SO3_1_1(R1(t)) + D_SO3_1_2(R1(t));
% Right side SO3 for R2
X_SO3_2 = @(R2) extractComp(X(R1(t),R2),4,6,1,3);
Y_SO3_2 = @(R2) extractComp(Y(R1(t),R2),4,6,1,3);
D_SO3_2_1 = @(R2) rot_metric(R2, rot_covar(dR2,X_SO3_2), Y_SO3_2);
D_SO3_2_2 = @(R2) rot_metric(R2, X_SO3_2, rot_covar(dR2,Y_SO3_2));
Right_SO3_2 = @(t) D_SO3_2_1(R2(t)) + D_SO3_2_2(R2(t));
% Resulting right side
Right = @(t) Right_SO3_1(t) + Right_SO3_2(t);

% Perform check
[Ft, dFt, appder] = funCheckDer(Left, Right, t, 'nodisplay');

% Additional Term for coupled dynamics
if bCoupledDynamics
    D_SO3_1_Extra = 1/2*trace( dR1(R1(t))'*...
        ( R1(t)*hat3(rot3_logDiffMat(R1(t),R2(t))*rot_vee(R2(t),dR2(R2(t)))) ) );
    dFt = dFt + D_SO3_1_Extra;
end

