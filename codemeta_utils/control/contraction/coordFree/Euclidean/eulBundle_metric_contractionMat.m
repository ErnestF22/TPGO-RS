function [M] = eulBundle_metric_contractionMat(kd, kv, m, n)
% Compute the metric on the tangent bundle of R^n ie: <D_{X}Y,X>
% NOTE: Assumes that X,Y are a vertical vector s.t. two tangent vectors on 
% R^n are stacked as [R^n; R^n].
% INPUTS:
%   X := vector field on TR^n evaluated at (p,u).
%   Y := vector field on TR^n evaluated at (p,u). Assume Closed Loop
%       VF
%   u := velocity vector function handle in R^n evaluated at p
%   kd, kv := scalar controller gains
%   m := an array of 3 values for the non natural metric on TR^n
%   n := scalar representing the dimension of the manifold
% OUTPUTS:
%   M := M matrix such that [zeta eta]*M*[zeta;eta] = <D_{X}Y,X>

%% Extract useful components
m1 = m(1); m2 = m(2); m3 = m(3);

%% Return the metric (for the cost function = norm(x))

if true
    %%% Expected Results on R^12
    M11 = -m2*kd*eye(n); M12 = ((-m3*kd+m1-m2*kv)*eye(n))/2;
    M21 = M12; M22 = (m2-m3*kv)*eye(n);
    %%% Expected Results on R^12
else
    %%%% Only induced components
    M11 = -m2*kd*eye(n); M12 = -m3*kd*eye(n);
    M21 = zeros(n); M22 = zeros(n);
    %%%% Only induced components
end

M = [M11 M12;...
    M21 M22];

end

