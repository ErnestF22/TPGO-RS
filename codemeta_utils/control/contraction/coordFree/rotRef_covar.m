function [nabla_X_Y] = rotRef_covar(X,Y,R,U,RRef,M,t,varargin)
% Compute the Levi-Civita connection compatible with the natural metric on
% TSO(3)xSO(3) evaluated at (R,U,RRef)
% INPUTS:
%   X(R,U,RRef), Y(R,U,RRef) := Vector fields on TSO(3)xSO(3)
%   R(t) := The evaluation point on the TSO(3) component
%   U(R) := the tangent vector at R which creates the curve on TSO(3) in
%       the form of [R; U]
%   RRef(t) := The evaluation point on the SO(3) component as a function of
%       time
%   M := A matrix of scaling factors (should be symmetric), [3x3] matrix,
%       where M(1,2:3) = M(2:3,1) = 0 (a natural metric)
%   t := time to evaluate the metric
% OUTPUTS:
%   nabla_X_Y(R) := covariant derivative of Y in the direction of X at 
%       point (R,U,RRef) at time t

% Extract the metric gains for each manifold
m1 = M(2,2); m2 = M(2,3); m3 = M(3,3);
m_nonnatural_TSO3 = [m1;m2;m3]; % Defined for rotBundle_covar_nonNatural()
% m4 = M(1,1); %This is unused since the connection on SO3 is wrt the natural metric with coeffiction 1. m4 plays a role in the metric as a scale factor.

%optional parameters
flagRigidRotControl = false; % Assumes Y is the controled closed loop VF for rigid body rotations
flagEvalAtRRef = false; % Compute covar on TSO3 at (RRef,URef) instead of (R,U)
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
        case 'evalatrrefontso3'
            % Compute the covar on TSO3 at (RRef,URRef) instead of (R,U). In
            % this case, this function needs to know the vector field Y
            % defined as a function of time with left-translation to RRef.
            flagEvalAtRRef = true;
            % Can compute along any point on the fiber at RRef since Y(R) 
            % and X(R) are viewed as left-invariant and lift-decomposable
            % at RRef (due to motion strictly on the reference manifold),
            % so choose 0
            U = @(RRef) zeros(3);
            % Store the implicitly depedenent Y(R(t),t) to be used
            % by the chain rule
            ivarargin=ivarargin+1;
            Y_LT_RRef = varargin{ivarargin};
        otherwise
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end
% Find the connection on TSO(3) (point on SO(3) fixed at RRef(t))
X_TSO3 = @(R,U) extractComp(X(R,U,RRef(t)),4,9,1,3);
Y_TSO3 = @(R,U) extractComp(Y(R,U,RRef(t)),4,9,1,3);
if flagRigidRotControl
    nabla_X_Y_TSO3 = rotBundle_covar_nonNatural(X_TSO3, ...
        Y_TSO3, U, m_nonnatural_TSO3,'rigidrot',kv);
else
    nabla_X_Y_TSO3 = rotBundle_covar_nonNatural(X_TSO3, ...
        Y_TSO3, U, m_nonnatural_TSO3);
end

% Add term for changes wrt dynamic reference trajectories
if flagEvalAtRRef
    % Computing at (RRef,URRef), require chain rule for dependencies on
    % (R(t),U(t))
    dY_dt = funApproxDer(Y_LT_RRef,t);
else
    % Computing at (R,U), require chain rule for dependencies on RRef(t)
    R0 = R(t); U0 = U(R0);
    dY_dt = funApproxDer(@(t) extractComp(Y(R0,U0,RRef(t)),4,9,1,3),t);
end
nabla_X_Y_TSO3 = @(R) nabla_X_Y_TSO3(R) + dY_dt;

% Connection on SO(3) (does not depend on m4, point on TSO(3) fixed at (R(t),U(R(t),t)) ))
X_SO3 = @(RRef) extractComp(X(R(t),U(R(t)),RRef),1,3,1,3);
Y_SO3 = @(RRef) extractComp(Y(R(t),U(R(t)),RRef),1,3,1,3);
nabla_X_Y_SO3 = rot_covar(X_SO3,Y_SO3);

% Resulting connection
nabla_X_Y = @(R) [nabla_X_Y_SO3(RRef(t));...
    nabla_X_Y_TSO3(R)];
end