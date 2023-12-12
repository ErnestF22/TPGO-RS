function [Ft, dFt, appder] = rotRef_metric_compatibility_nonNatural2(R,U,RRef,X,Y,Z,M_nonnatural,t,varargin)
% Check if d/dt<X,Y> = <D_{Z_Dot}X,Y> + <X,D_{Z_dot}Y> evaluated at
% (R,U,RRef,t). ASSUME Y = Z_Dot FOR OUR CONTRACTION FORM
% ASSUMES THAT VECTOR FIELDS ARE DEFINED IN THE "NONNATURAL" COORDINATES. 
% TEST TO SEE IF THE NONNATURAL RESULTS CAN BE DERIVIED FROM THE NATURAL 
% ONES
% NOTE: This compat. test is different from 
% rotRef_metric_compatibility_nonNatural.m in that the tangent vector on 
% the reference manifold is at the bottom instead of the top of
% the [9x9] vector.
% INPUTS:
%   R(t) := Curve on SO(3) for the TSO(3) componenet
%   U(R,t) := Tangent vector on SO(3) for the TSO(3) component
%   RRef(t) := Curve on SO(3) for SO(3) component
%   X(R,U,RRef), Y(R,U,RRef), Z(R,U,RRef) := Vector fields on TSO(3)xSO(3),
%       NOTE: Z=d/dt Z_curve; Also X is the closeloop vector field
%   M_nonnatural := A matrix of scaling factors (should be symmetric), [3x3] matrix
%   t := time to evaluate the metric
% OUTPUTS:
%   Ft := numerical evaluation of <X,Y>
%   dFt := analytical evaluation of <D_{Z_dot}X, Y> + <X, D_{Z_dot}Y>
%   appder := numerical approximation of d/dt <X,Y>

% Define Useful parameters
Z_curve = @(t) [R(t);U(R(t),t);RRef(t)]; % Point of evaluation
U_R = @(R) U(R,t); % Redefined as function of R for covar on TSO3

% Optional Parameters (No action, pass onto rotRef_covar_nonNatural2)
flagUseContractionMatrix = false;
ivarargin=1;
len_varargin = length(varargin);
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'matrix'
            % Require gains as vector [kd;kv;kp]
            % Compute <D_Y_X,Y> using rotRef_contractionMat2 assuming that
            % Y=Z=d/dt Z_curve
            flagUseContractionMatrix = true;
            ivarargin=ivarargin+1;
            gains = varargin{ivarargin};
        otherwise
            % ignore the entry and increase index until next entry is
            % a string
            while ivarargin+1 < len_varargin && ~ischar(varargin{ivarargin+1})
                ivarargin=ivarargin+1;
            end
    end
    ivarargin=ivarargin+1;
end

% Define the metric as a function of t
Xt = @(t) X(R(t),U(R(t),t),RRef(t));
Yt = @(t) Y(R(t),U(R(t),t),RRef(t));
Left = @(t) rotRef_metric_nonNatural2(Z_curve(t),Xt(t),Yt(t),M_nonnatural);

% Compute the right side <D_Y_X,Y> + <D_Y_Y,X>
% <D_Y_X,Y>
if flagUseContractionMatrix
    M_contraction = rotRef_contractionMat2(0,gains(1),gains(2),gains(3),M_nonnatural);
    Z_t = Z(R(t),U(R(t),t),RRef(t));
    zeta=rot_vee(R(t),extractComp(Z_t,1,3,1,3));
    eta=rot_vee(R(t),extractComp(Z_t,4,6,1,3));
    nu=rot_vee(RRef(t),extractComp(Z_t,7,9,1,3));
    g_DZX_Y = @(t) [zeta;eta;nu]'*M_contraction(R(t),rot_vee(R(t),U_R(R(t))),RRef(t))*[zeta;eta;nu];
else
    [D_Z_X_nonnatural,D_Z_X_natural] = rotRef_covar_nonNatural2(R,U,RRef,t,Z,X,M_nonnatural,varargin{:});
    g_DZX_Y = @(t) rotRef_metric_nonNatural2(Z_curve(t),D_Z_X_nonnatural,Yt(t),M_nonnatural);
end
% <D_Y_Y,X>
[D_Z_Y_nonnatural,D_Z_Y_natural] = rotRef_covar_nonNatural2(R,U,RRef,t,Z,Y,M_nonnatural);
g_DZY_X = @(t) rotRef_metric_nonNatural2(Z_curve(t),D_Z_Y_nonnatural,Xt(t),M_nonnatural);
g_natural_total_der_fun = @(t) g_DZX_Y(t) + g_DZY_X(t);
%% Perform check
[Ft, dFt, appder] = funCheckDer(Left, g_natural_total_der_fun, t, 'nodisplay');
end