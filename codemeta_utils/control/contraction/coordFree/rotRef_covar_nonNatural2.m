function [D_X_Y_nonnatural,D_X_Y_natural] = rotRef_covar_nonNatural2(R,U,RRef,t,X,Y,M_nonnatural,varargin)
% Compute the Levi-Civita connection compatible with the nonnatural metric
% M on TSO(3)xSO(3) at (R,U,RRef)|_{t=t}. To do so, we perform a coordinate
% change from "nonnatural" to "natural" using Schur's complement to find
% the transformation matrix.
% NOTE: This covar is different from rotRef_SchurComplement.m in that the
% tangent vector on the reference manifold is at the bottom instead of the
% top of the [9x9] vector.
% INPUTS:
%   R(t) := The evaluation point on the TSO(3) component
%   U(R,t) := The tangent vector at R which creates the curve on TSO(3) in
%       the form of [R; U]. NOTE: the tangent vector must remain in
%       T_{R}SO(3).       
%   RRef(t) := The evaluation point on the SO(3) component. dRRef/dt should
%       produce the tangent vector in the "natural" coordinates.
%   t := time of evaluation
%   X(R,U,RRef), Y(R,U,RRef) := Vector fields on TSO(3)xSO(3), represented
%       as [9x3] vector
%   M_nonnatural := A [3x3] pos. def. matrix representing the metric gains
% OUTPUTS:
%   D_X_Y_nonnatural := The covariant derivative of Y along X evaluated at
%       (R,U,RRef)|_{t=t}

% Compute the point (R,U,RRef)
R_t = R(t); U_t = U(R_t,t); RRef_t = RRef(t);
% Redefine U as function of R
U_R = @(R) U(R,t);
% Extract the 3 tangent vector of each vector field
X_R_nonnatural = @(R,U,RRef) extractComp(X(R,U,RRef),1,3,1,3); % \in T_{R}SO3
X_U_nonnatural = @(R,U,RRef) extractComp(X(R,U,RRef),4,6,1,3); % \in T_{R}SO3
X_RRef_nonnatural = @(R,U,RRef) extractComp(X(R,U,RRef),7,9,1,3); % \in T_{RRef}SO3
Y_R_nonnatural = @(R,U,RRef) extractComp(Y(R,U,RRef),1,3,1,3); % \in T_{R}SO3
Y_U_nonnatural = @(R,U,RRef) extractComp(Y(R,U,RRef),4,6,1,3); % \in T_{R}SO3
Y_RRef_nonnatural = @(R,U,RRef) extractComp(Y(R,U,RRef),7,9,1,3); % \in T_{RRef}SO3

% Optional Parameters
flagMetricNatural = false;
ivarargin=1;
len_varargin = length(varargin);
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'naturalmetric'
            % Use this option if M_nonnatural is actually a natural metric
            % with transformation matrix J (MUST PROVIDE)
            flagMetricNatural = true;
            ivarargin=ivarargin+1;
            J = varargin{ivarargin};
            M_natural = M_nonnatural;
        otherwise
            % ignore the entry and increase index until next entry is
            % a string
            while ivarargin+1 < len_varargin && ~ischar(varargin{ivarargin+1})
                ivarargin=ivarargin+1;
            end
    end
    ivarargin=ivarargin+1;
end

% Find the natural metric and its coordinate transformation matrix using
% Schur's complement
if ~flagMetricNatural
    % Find transformation matrix since M_nonnatural is in "nonnatural"
    % coordinates
    [J,M_natural] = rotRef_SchurComplement2(M_nonnatural);
end
% extact the metric on TSO(3) for rotBundle_metric_nonNatural.m
m_natural = [M_natural(1,1);M_natural(1,2);M_natural(2,2)];

% Compute nabla_X_Y on the SO3 component
% First define the transformed tangent vectors
X_RRef_natural = @(R,U,RRef) X_RRef_nonnatural(R,U,RRef)...
    +J(3,1)*RRef*R'*X_R_nonnatural(R,U,RRef)...
    +J(3,2)*RRef*R'*X_U_nonnatural(R,U,RRef);
X_SO3_natural_RRef = @(RRef) X_RRef_natural(R_t,U_t,RRef);
Y_RRef_natural = @(R,U,RRef) Y_RRef_nonnatural(R,U,RRef)...
    +J(3,1)*RRef*R'*Y_R_nonnatural(R,U,RRef)...
    +J(3,2)*RRef*R'*Y_U_nonnatural(R,U,RRef);
Y_RRef_natural_RRef = @(RRef) Y_RRef_natural(R_t,U_t,RRef);

% Compute Y_RRef_natural_RRef along RRef and the symmetric part as a function
% of X_RRef_natural (the transformed vector) using rot_covar(X,Y,dY)
d_XSO3_dRRef = funApproxDer(@(t) Y_RRef_natural(R_t,U_t,RRef(t)),t);
D_X_Y_SO3 = rot_covar(X_SO3_natural_RRef,Y_RRef_natural_RRef,@(RRef) d_XSO3_dRRef);
% Chain rule wrt R(t)|_{t=t}
dJY_dt_SO3 = funApproxDer(@(t) Y_RRef_natural(R(t),U(R(t),t),RRef_t),t);
D_X_Y_SO3 = @(RRef) D_X_Y_SO3(RRef) + dJY_dt_SO3;

% Compute nabla_X_Y on the TSO3 component (The tangent vectors are
% unchanged with this transformation)
X_TSO3_natural = @(R,U,RRef) [X_R_nonnatural(R,U,RRef); ...
    X_U_nonnatural(R,U,RRef)];
% Redefine as function of (R,U) to compute the covar on TSO(3)
X_TSO3_natural_RU = @(R,U) X_TSO3_natural(R,U,RRef_t);
Y_TSO3_natural = @(R,U,RRef) [Y_R_nonnatural(R,U,RRef);...
    Y_U_nonnatural(R,U,RRef)];
% Redefine as function of (R,U) to compute the covar on TSO(3)
Y_TSO3_natural_RU = @(R,U) Y_TSO3_natural(R,U,RRef_t);
% Compute the covar on TSO(3)
D_X_Y_TSO3 = rotBundle_covar_nonNatural(X_TSO3_natural_RU,...
    Y_TSO3_natural_RU, U_R, m_natural,varargin{:});

% Chain rule wrt RRef(t)|_{t=t}
dJY_dt = funApproxDer(@(t) Y_TSO3_natural(R_t,U_t,RRef(t)),t);
D_X_Y_TSO3 = @(R) D_X_Y_TSO3(R) + dJY_dt;

% Construct the "natural" covar
D_X_Y_natural = [D_X_Y_TSO3(R_t);D_X_Y_SO3(RRef_t)];
% Do an inverse coordinate transform to find the "nonnatural" covar
J_inv = inv(J);
D_X_Y_R_natural = extractComp(D_X_Y_natural,1,3,1,3);
D_X_Y_U_natural = extractComp(D_X_Y_natural,4,6,1,3);
D_X_Y_RRef_natural = extractComp(D_X_Y_natural,7,9,1,3);
D_X_Y_nonnatural = [D_X_Y_R_natural;...
    D_X_Y_U_natural;...
    D_X_Y_RRef_natural...
    + J_inv(3,1)*RRef_t*R_t'*D_X_Y_R_natural...
    + J_inv(3,2)*RRef_t*R_t'*D_X_Y_U_natural];
end

