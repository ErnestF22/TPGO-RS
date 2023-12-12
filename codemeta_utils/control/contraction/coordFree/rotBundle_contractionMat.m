function [M, output] = rotBundle_contractionMat(U, lambda, kd, kv, m, varargin)
% Computes the matrix form of the M matrix such that 
% [zeta eta]*M*[zeta;eta] == rotBundle_metric_nonNatural_contraction, where
% zeta and eta are the R^3 representation of the horizontal and vertical
% tangent vectors defined in <Del_{Vx}X + lambda*Vx, Vx>
% ASSUMES: Vx = [R*hat3(zeta);R*hat3(eta)];
% ASSUMES: X = [U;kd*rot_log(R,eye(3)-kv*U];
% INPUTS:
%   U(R) := the tangent vector at R which creates the curve on TSO(n) in
%       the form of [R; U]
%   lambda := scalar for the exponential convergence rate
%   kd := scalar position gain
%   kv := scalar velocity gain
%   m := an array of 3 values for the non natural metric on TSO(3)
% OUTPUTS:
%   M(R) := M matrix such that 
%       [zeta eta]*M*[zeta;eta] == rotBundle_metric_nonNatural_contraction

flagAddOutput = false;
flagOutputSymmetric = false;

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'extras'
            flagAddOutput = true;
        case 'sym'
            flagOutputSymmetric = true;
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
omega = @(R) rot_vee(R, U(R));

%% Compute the components of the M matrix
if flagOutputSymmetric
    M11 = @(R) m2/4*hat3(omega(R))^2 ...
        - m2*kd/2*(rot3_logDiffMat(eye(3),R)+rot3_logDiffMat(eye(3),R)')...
        + lambda*m1*eye(3);
else
    M11 = @(R) m2/4*hat3(omega(R))^2 ...
        - m2*kd*rot3_logDiffMat(eye(3),R)...
        + lambda*m1*eye(3);
end
M21 = @(R) -m2/4*hat3(omega(R))...
    + m3/8*hat3(omega(R))^2 ...
    + m3*kv/4*hat3(omega(R))...
    + m3*kd/4*rot_log(eye(3),R)...
    - m3*kd/2*rot3_logDiffMat(eye(3),R)...
    + lambda*m2*eye(3);
M22 = @(R) lambda*m3*eye(3);

N21 = @(R) (m1-m2*kv)/2*eye(3);
N22 = @(R) (m2-m3*kv)*eye(3);
N = @(R) [zeros(3) N21(R)';...
    N21(R) N22(R)];

%% Define outputs
M = @(R) [M11(R) M21(R)';...
    M21(R) M22(R)] + N(R);

if flagAddOutput
    output.M11 = M11;
    output.M21 = M21;
    output.M22 = M22;
end

end

