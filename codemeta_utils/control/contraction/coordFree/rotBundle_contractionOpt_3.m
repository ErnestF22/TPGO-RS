function [ m, flag ] = rotBundle_contractionOpt_3( maxD, kv, kd, beta, magW )
% Gershgorin disc method to bound [2 x 2] matrix, bound the D\omega_D3
% eigenvalue by:
% sqrt(3)/8*magW*(m3*sqrt(4*kv^2+3*magW^2)+abs(m2)*sqrt(3*magW^2+4))-3/8*m2*magW^2
% INPUT:
%   maxD := largest distance error on SO(3)
%   kv := positive scalar velocity gain
%   kd := positive scalar position gain
%   beta := scalar min exponential convergence rate (ie exp(-betaIn*t))
%   magW := largest linear speed, i.e. abs(w1)
% OUTPUT:
%   m := [m1 m2; m2 m3] values
%   flag := true means the system is feasible

%% Default Parameters
ZERO_TOL = 1e-5;
flag = false; %set to false by default

%% Optimize using gershgorin bounds
cvx_begin
    variable m(2,2) semidefinite 
    
    % Define some useful constants
    b = (beta*m(1,2) + (m(1,1)-m(1,2)*kv)/2);
    theta = maxD; % max geodesic distance on SO(3)
    Mw_D5 = sqrt(3)/8*magW*(m(2,2)*sqrt(4*kv^2+3*magW^2)+abs(m(1,2))*sqrt(3*magW^2+4))-3/8*m(1,2)*magW^2;
    MR_D1 = beta*m(1,1)-m(1,2)*kd + abs(b-m(2,2)*kd/2);
    MR_D2 = beta*m(1,1)-m(1,2)*kd*theta/2*cot(theta/2) + ...
        abs(b-m(2,2)*kd/2*theta/2*cot(theta/2));
    MR_D4 = (beta*m(2,2)+m(1,2)-m(2,2)*kv) + abs(b-m(2,2)*kd/2);
    MR_D5 = (beta*m(2,2)+m(1,2)-m(2,2)*kv) + ...
        abs(b-m(2,2)*kd/2*theta/2*cot(theta/2));
    % Define the different combinations of the possible discs for Mw + MR
    M = [   Mw_D5 + MR_D1;...
            Mw_D5 + MR_D2;...
            Mw_D5 + MR_D4;...
            Mw_D5 + MR_D5;...
        ];
    minimize (max(M))
    subject to
        m(1,1)==1 % Fix to a particular solution since any scale of m's is 
                  % also a valid solution
cvx_end
% M

% Check if this is a feasible solution
if (max(M) < ZERO_TOL && det(m) > ZERO_TOL)
    flag = true;
end
end

