function [Ft, dFt, appder] = rotRef_metric_compatibility_nonNatural(R,U,RRef,X,Y,Z,M_nonnatural,t,varargin)
% Check if d/dt<X,Y> = <D_{Z_Dot}X,Y> + <X,D_{Z_dot}Y> evaluated at
% (R,U,RRef,t).
% ASSUMES THAT VECTOR FIELDS ARE DEFINED IN THE "NONNATURAL" COORDINATES. 
% TEST TO SEE IF THE NONNATURAL RESULTS CAN BE DERIVIED FROM THE NATURAL 
% ONES
% INPUTS:
%   R(t) := Curve on SO(3) for the TSO(3) componenet
%   U(R,t) := Tangent vector on SO(3) for the TSO(3) component
%   RRef(t) := Curve on SO(3) for SO(3) component
%   X(R,U,RRef), Y(R,U,RRef), Z(R,U,RRef) := Vector fields on TSO(3)xSO(3),
%       NOTE Z=d/dt Z_curve
%   M_nonnatural := A matrix of scaling factors (should be symmetric), [3x3] matrix
%   t := time to evaluate the metric
% OUTPUTS:
%   Ft := numerical evaluation of <X,Y>
%   dFt := analytical evaluation of <D_{Z_dot}X, Y> + <X, D_{Z_dot}Y>
%   appder := numerical approximation of d/dt <X,Y>

% Define Useful parameters
Z_curve = @(t) [RRef(t);R(t);U(R(t),t)]; % Point of evaluation
U_R = @(R) U(R,t); % Redefined as function of R for covar on TSO3

% Optional Parameters
flagRigidRotControl = false; % Assumes Y is the controled closed loop VF for rigid body rotations
kv = 0;
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'rigidrot'
            % Add adjustment terms since Y has some dependency on its
            % position along the fiber. IE Y is not completely
            % lift-decomposable.
            % Assumes VF is rigid body rotation with controller 
            % u = kd*log(R,eye(3)) - kv*Yh
            flagRigidRotControl=true;
            % Must also send in a value for kv gain!
            ivarargin=ivarargin+1;
            kv = varargin{ivarargin};
        otherwise    
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

% Define the metric as a function of t
Xt = @(t) X(R(t),U(R(t),t),RRef(t));
Yt = @(t) Y(R(t),U(R(t),t),RRef(t));
Left = @(t) rotRef_metric_nonNatural(Z_curve(t),Xt(t),Yt(t),M_nonnatural);

% Compute the right side <D_Y_X,Y> + <D_Y_Y,X>
% <D_Y_X,Y>
if flagRigidRotControl
    [D_Z_X_nonnatural,D_Z_X_natural] = rotRef_covar_nonNatural(R,U,RRef,t,Z,X,M_nonnatural,'rigidrot',kv,'metriccompatibility');
else
    [D_Z_X_nonnatural,D_Z_X_natural] = rotRef_covar_nonNatural(R,U,RRef,t,Z,X,M_nonnatural,'metriccompatibility');
end
g_DZX_Y = @(t) rotRef_metric_nonNatural(Z_curve(t),D_Z_X_nonnatural,Yt(t),M_nonnatural);
[D_Z_Y_nonnatural,D_Z_Y_natural] = rotRef_covar_nonNatural(R,U,RRef,t,Z,Y,M_nonnatural,'metriccompatibility');
g_DZY_X = @(t) rotRef_metric_nonNatural(Z_curve(t),D_Z_Y_nonnatural,Xt(t),M_nonnatural);
g_natural_total_der_fun = @(t) g_DZX_Y(t) + g_DZY_X(t);
%% Perform check
[Ft, dFt, appder] = funCheckDer(Left, g_natural_total_der_fun, t, 'nodisplay');
[D_Y_X_nonnatural,D_Y_X_natural] = rotRef_covar_nonNatural(R,U,RRef,t,Y,X,M_nonnatural,'rigidrot',kv);
end