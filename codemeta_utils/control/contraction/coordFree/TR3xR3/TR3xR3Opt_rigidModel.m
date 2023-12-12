function [ dx ] = TR3xR3Opt_rigidModel( x, xd, m, F, gradf_ref, M_contract, M_contract_der_gains, M_nn )
% Last Edited: Oct. 27 2020 by Bee Vang
%Create a rigid point mass model in R^3 with a reference trajectory. The
%model also solves an optimzation problem to determine how the gains should
%change.
%Model is given by m*ddx=F
% INPUT:
%   x := current state [12x1], [curr. pos.; curr. vel.;ref. pos;kd;kv;kref]
%   xd := desired final pos [3x1]
%   m := mass
%   F := control input
%   gradf_ref(x,xd) := function handle to gradient of cost function for
%       reference trajectory
%   M_contract(x,v,xref,kd,kv,kref) : function handle to the contraction metric, the
%       \beta parameter should be set to 0
%   M_contract_der_gains(x,v,xref,kd_dot,kv_dot,kref_dot) := function
%       handle to the derivative of the contraction metric wrt the gains
%   M_nn := A [3x3] pos. def. matrix representing the metric gains.

% Define parameters limits
GAINS_CBF_RATE = 1e-3;
GAIN_UPPER_BOUND = 150;
GAIN_LOWER_BOUND = 0;

% Extract states
xpos = x(1:3); xvel = x(4:6); xref = x(7:9);
kd = x(10); kv = x(11); kref = x(12);
% Compute der of augmented system
dx = [zeros(3) eye(3) zeros(3);zeros(3,9);zeros(3,9)]*x(1:9) ... % state
    + [zeros(3);eye(3);zeros(3)]*F*1/m ... % control
    + [zeros(6,1);-kref*gradf_ref(x(7:9),xd)]; % reference traj.

% Compute der of gains by solving optimization problem
% Compute contraction matrix at current state
M_contract_curr = M_contract(xpos,xvel,xref,kd,kv,kref);
M_der_fun = @(kd_dot,kv_dot,kref_dot) M_contract_der_gains(xpos,xvel,xref,kd_dot,kv_dot,kref_dot);
cvx_begin sdp quiet
    variable kd_dot
    variable kv_dot
    variable kref_dot
    variable b nonnegative
    
    maximize ( b )
    subject to
        M_contract_curr + M_der_fun(kd_dot,kv_dot,kref_dot) + kron(b*M_nn,eye(3)) <= 0
        % CBF lower bound == 0
        -kd_dot <= (kd-GAIN_LOWER_BOUND)^3*GAINS_CBF_RATE
        -kv_dot <= (kv-GAIN_LOWER_BOUND)^3*GAINS_CBF_RATE
        -kref_dot <= (kref-GAIN_LOWER_BOUND)^3*GAINS_CBF_RATE
        % Bound rate of change
        abs(kd_dot) <= 100
        abs(kv_dot) <= 100
        abs(kref_dot) <= 100
        % CBF upper bound
        kd_dot <= (GAIN_UPPER_BOUND-kd)^3*GAINS_CBF_RATE
        kv_dot <= (GAIN_UPPER_BOUND-kv)^3*GAINS_CBF_RATE
        kref_dot <= (GAIN_UPPER_BOUND-kref)^3*GAINS_CBF_RATE       
cvx_end
b
% if no solution, set rates = 0 
if  ~contains(lower(cvx_status),'solved','IgnoreCase',true)
    % Otherwise, keep the latest gains
    kd_dot = 0;
    kv_dot = 0;
    kref_dot = 0;
end

% update the gain rates
dx = [dx;kd_dot;kv_dot;kref_dot];
end

