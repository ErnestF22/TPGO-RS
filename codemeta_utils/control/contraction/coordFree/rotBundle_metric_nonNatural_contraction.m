function [D, output] = rotBundle_metric_nonNatural_contraction(U, Vx, X, lambda, kd, kv, m, varargin)
% Compute the non natural metric as a factored quadratic form
% INPUTS:
%   U(R) := the tangent vector at R which creates the curve on TSO(n) in
%       the form of [R; U]
%   Vx(R,U) := Arbitrary tangent vector field on TSO(3) given as matrices 
%       (IE [6 x 3])
%   X(R,U) := System dynamics as a vector field on TSO(3) given as matrices 
%       (IE [6 x 3])
%   X_dot(R,U) := System dynamics as a vector field on TSO(3) given as matrices 
%       (IE [6 x 3])
%   m := an array of 3 values for the non natural metric on TSO(3)
% OUTPUTS:
%   D(R) := Result of g_bar(Del_{Vx}X + lambda*Vx, Vx) on TSO(3)

flagAddOutput=false;

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'extras'
            flagAddOutput=true;
        otherwise
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

%% Standard Variables
if (length(m) ~= 3)
    error('m must be an array of length 3')
end
m1 = m(1); m2 = m(2); m3 = m(3); % These should produce a pos. def. symm. matrix M


%% Split Vx, X into vertical and horizontal components
PU = @(R) [R; U(R)];
Vxh = @(R) rotBundle_extractHoriz(PU(R), Vx(R, U(R)));
zeta = @(R) rot_vee(R, Vxh(R));
Vxv = @(R) rotBundle_extractVert(PU(R), Vx(R, U(R)));
eta = @(R) rot_vee(R,Vxv(R));
Xh = @(R) rotBundle_extractHoriz(PU(R), X(R, U(R)));
omega = @(R) rot_vee(R, Xh(R));
% Xv = @(R) rotBundle_extractVert(PU(R), X(R, U(R)));
% theta = @(R) rot_vee(R, Xv(R));

% gammaVxh = @(R,t) rot_exp(R,t*Vxh(R));
% Xdoth = @(R) funApproxDer( @(t) Xh(gammaVxh(R,t)), 0);
% Xdotv = @(R) funApproxDer( @(t) Xv(gammaVxh(R,t)), 0);

%% Compute components of the metric
gbar_Del_Xh_Yh_Zh= @(R) +m2/4*zeta(R)'*hat3(omega(R))*hat3(omega(R))*zeta(R);
gbar_Del_Xh_Yh_Zv= @(R) - 1/2*m2*eta(R)'*hat3(omega(R))*zeta(R) +m3/8*eta(R)'*hat3(omega(R))*hat3(omega(R))*zeta(R);

% gbar_Del_Xh_Yv_Zh = @(R) m2/2*trace(hat3(zeta(R))'*(R'*Xdotv(R))); % Numeric solution works...
gbar_Del_Xh_Yv_Zh = @(R) m2*kd*zeta(R)'*rot3_logDiffMat(eye(3),R)*zeta(R);

% gbar_Del_Xh_Yv_Zv = @(R) m3/2*trace(hat3(eta(R))'*(R'*Xdotv(R))); % Numeric solution works...
gbar_Del_Xh_Yv_Zv = @(R) m3*kv/2*eta(R)'*hat3(omega(R))*zeta(R)...
    - m3*kd/2*eta(R)'*rot_log(eye(3),R)*zeta(R)...
    + m3*kd*eta(R)'*rot3_logDiffMat(eye(3),R)*zeta(R);
    
gbar_Del_Xv_Yh_Zh = @(R) +m3/8*eta(R)'*hat3(omega(R))*hat3(omega(R))*zeta(R);
gbar_LamVxh_Vxh = @(R) lambda*m1*norm(zeta(R))^2;
gbar_LamVxv_Vxh = @(R) 2*lambda*m2*zeta(R)'*eta(R);
gbar_LamVxv_Vxv = @(R) lambda*m3*norm(eta(R))^2;

%% Resulting metric is sum of all components above
D = @(R) gbar_Del_Xh_Yh_Zh(R) + gbar_Del_Xh_Yh_Zv(R)...
    + gbar_Del_Xh_Yv_Zh(R) + gbar_Del_Xh_Yv_Zv(R)...
    + gbar_Del_Xv_Yh_Zh(R) + gbar_LamVxh_Vxh(R)...
    + gbar_LamVxv_Vxh(R) + gbar_LamVxv_Vxv(R);

if flagAddOutput
    output.gbar_Del_Xh_Yh_Zh = gbar_Del_Xh_Yh_Zh;
    output.gbar_Del_Xh_Yh_Zv = gbar_Del_Xh_Yh_Zv;
    output.gbar_Del_Xh_Yv_Zh = gbar_Del_Xh_Yv_Zh;
    output.gbar_Del_Xh_Yv_Zv = gbar_Del_Xh_Yv_Zv;
    output.gbar_Del_Xv_Yh_Zh = gbar_Del_Xv_Yh_Zh;
    % This D is computed without the +lambda*Vx term
    output.D = @(R) gbar_Del_Xh_Yh_Zh(R) + gbar_Del_Xh_Yh_Zv(R)...
    + gbar_Del_Xh_Yv_Zh(R) + gbar_Del_Xh_Yv_Zv(R)...
    + gbar_Del_Xv_Yh_Zh(R);

%     output.Xdotv = Xdotv;
end
