function [dparams, delta] = rotRefOptAll_optProblem(x,varargin)
% Last Edited: Nov 13, 2020 by Bee Vang
% Solve optimization problem to find changes for gains, metric, and
% convergence rate
% INPUTS:
%   x := [34x1] vector with states [R,w,RRef,kd,kv,kref,M_nn,beta]
% OUTPUTS:
%   dparams := two possible options
%       1) [13x1] vector representing the optimal change in the [gains,
%           metric, and convergence rate]
%       2) [1x1] the optimal beta assuming fixed gains and metric
%   delta := relaxation parameter for control CLF

%% Set Parameters
GAINS_CLF_RATE = 1e-3;
GAINS_CBF_SCALE = 1;
GAIN_UPPER_BOUND = 150;
GAIN_LOWER_BOUND = 0;
BETA_LOWER_BOUND = 0; % CBF to prevent beta <= 0 (system diverge)
BETA_CBF_SCALE = 1;
flagDynamicGains = true;
flagDynamicMetric = true;
flagOptConvergenceRateOnly = false; % Solve for beta assuming gains and metric are fixed
EIG_TOL = 1e-4; % Set the min eigenvalue for metric tensor
METRIC_SCALE_FACTOR = 1.0;
%% Load Optional Parameters
ivarargin=1;
len_varargin = length(varargin);
while ivarargin<=len_varargin
    switch lower(varargin{ivarargin})
        case 'contraction_matrix_handle'
            % Extracts the contraction matrix handle (instance of
            % rotRefOpt_contractionMat) to save generation time
            ivarargin=ivarargin+1;
            contractionMat_handle = varargin{ivarargin};
            % Also require the handle to the derivative of the contraction
            % matrix
            ivarargin=ivarargin+1;
            contractionMat_kdot_handle = varargin{ivarargin};
        case 'gains_cbf_scale'
            % Update CBF desired convergence rate
            ivarargin=ivarargin+1;
            GAINS_CBF_SCALE = varargin{ivarargin};
        case 'gains_upper_bound'
            % Update max gain bounds
            ivarargin=ivarargin+1;
            GAIN_UPPER_BOUND = varargin{ivarargin};
        case 'gains_lower_bound'
            % Update min gain bounds, must be positive
            ivarargin=ivarargin+1;
            GAIN_LOWER_BOUND = varargin{ivarargin};
        case 'gains_clf_rate'
            % Update desired exp conv rate of norm CLF
            ivarargin=ivarargin+1;
            GAINS_CLF_RATE = varargin{ivarargin};
        case 'static_metric'
            % Do not allow the metric parameters to change
            flagDynamicMetric = false;
        case 'static_gains'
            % Do not allow the gains to change
            flagDynamicGains = false;
        case 'opt_convergence_rate_only'
            % Flag to solve for optimial convergence rate \beta at the
            % particular state assuming the gains and metric are fixed
            flagOptConvergenceRateOnly = true;
        otherwise
            % ignore the entry and increase index until next entry is
            % a string
            while ivarargin+1 <= len_varargin && ~ischar(varargin{ivarargin+1})
                ivarargin=ivarargin+1;
            end
    end
    ivarargin=ivarargin+1;
end

if GAIN_LOWER_BOUND < 0
    GAIN_LOWER_BOUND = 0;
end

%% Extract current state and generate function handles if needed
[R,w,RRef,kd,kv,kref,M_nonnatural,b] = rotDynOptAll_stateUnpack(x);
gains = [kd;kv;kref];

% Generate the contraction matrix handle if one wasn't provided (not ideal)
if ~exist('contractionMat_handle','var')
    % Current contraction matrix (R,w,RRef,gains,M_nn,beta)
    contractionMat_handle = rotRefOptAll_contractionMat('sym');
    % Der of contraction matrix with respect to the gains (R,w,RRef,gains_dot,M_nn)
    contractionMat_kdot_handle = rotRefOpt_contractionMat_kdot('sym');
    % Der of contraction matrix with respect to the metric parameters (R,w,RRef,gains,M_nn,M_nn_dot,beta)
    contractionMat_mdot_handle = rotRefOpt_contractionMat_mdot('sym');
    % Note the der of contraction matrix wrt \beta is defined below (M_nn,b_dot)
end

% Define states as function of time
R_t = @(t) rot_exp(R,t*R*hat3(w));
w_t = @(t) w + t*rot_vee(R,kd*rot_log(R,RRef)-kv*R*hat3(w));
RRef_t = @(t) rot_exp(RRef,t*kref*rot_log(RRef,eye(3)));


%% Solve optimization problem
% Solve optimization problem: max b_dot, allow gains, metric, and beta to
% change
if flagOptConvergenceRateOnly
    % Solve for optimal beta with fixed gains and metric
    % Set convergence rate to zero and solve for it in optimization problem
    contractionMatrix_current = contractionMat_handle(R,w,RRef,gains,M_nonnatural,0);
    cvx_begin sdp quiet
        variable b nonnegative
        
        maximize (b)
        
        subject to
            contractionMatrix_current + kron(b*M_nonnatural,eye(3)) <= 0
    cvx_end
    % if no solution, set rates to zero
    if  ~contains(lower(cvx_status),'solved','IgnoreCase',true)
        fprintf('Infeasible\n')
        b = 0;
    end
    
    dparams = b;
else
    % Compute the current contraction matrix (including beta)
    contractionMatrix_current = contractionMat_handle(R,w,RRef,gains,M_nonnatural,b);
    cvx_begin sdp quiet
        if ~flagDynamicGains
            k_dot = zeros(3,1);
            delta = 0;
        else
            variable k_dot(3,1) % gain rates
            variable delta % CLF relaxation factor (positive means larger control effort, negative means less control effort)
        end
        
        if ~flagDynamicMetric
            m_dot = zeros(5,1);
        else
            variable m_dot(5,1) % [m1_dot;m2_dot;m3_dot;m5_dot;m6_dot], assume m4 is constant
        end
        variable b_dot % convergence change rate

        % Define the positive definite metric parameters rate of change matrix
        dMnn_dt = [m_dot(1), m_dot(2), m_dot(5);...
            m_dot(2), m_dot(3), m_dot(4);...
            m_dot(5), m_dot(4), 0];

%         optVars = [k_dot;m_dot;b_dot;delta];
%         Q = zeros(10); % Positive Semi-definite cost
%         P = [0;0;0;0;0;0;0;0;1;-1];
%         maximize ( optVars'*Q*optVars + P'*optVars ) 
        maximize ( b_dot - delta*1e-5) % looking to increase convergence rate
%         maximize ( b_dot - delta*1e-5 - 2*k_dot'*k_dot) % looking to increase convergence rate

        subject to
            % contraction constraint for M matrix:
            % M + dM/d(R,w,RRef) + dM/dk + dM/dm + dM/db <= 0
            contractionMatrix_current ... % M
                + funApproxDer(@(t) contractionMat_handle(R_t(t),w_t(t),RRef_t(t),gains,M_nonnatural,b),0) ... % dM/d(R,w,RRef)
                + contractionMat_kdot_handle(R,w,RRef,k_dot,M_nonnatural) ... % dM/dk
                + contractionMat_mdot_handle(R,w,RRef,gains,M_nonnatural,dMnn_dt,b) ... % dM/dm
                + kron(b_dot*M_nonnatural,eye(3)) <= 0 % dM/db
            
            % CBF to prevent beta from going below 0
            b_dot + BETA_CBF_SCALE*(b - BETA_LOWER_BOUND) >= 0
%             abs(b_dot) <= 1e4
            % If dynamic gains
            if flagDynamicGains
                % CLF for min norm controller
                2*[kd*rot_vee(R,rot_log(R,RRef))'*rot_vee(R,rot_log(R,RRef)) - kv*w'*rot_vee(R,rot_log(R,RRef)),...
                    -kd*rot_vee(R,rot_log(R,RRef))'*w + kv*(w'*w),...
                    0]*k_dot...
                    + funApproxDer(@(t) norm(kd*rot_log(R_t(t),RRef_t(t))-kv*R_t(t)*hat3(w_t(t)))^2,0) <= delta
%                     + GAINS_CLF_RATE*norm(kd*rot_log(R,RRef) - kv*R*hat3(w))^2 - delta <= 0

                % CBF lower bound == 0
                k_dot + GAINS_CBF_SCALE*(gains - GAIN_LOWER_BOUND) >= 0
                % Bound rate of change
%                 abs(k_dot) <= 13 % This is already bounded by CBF
                % CBF upper bound
                -k_dot + GAINS_CBF_SCALE*sqrt(GAIN_UPPER_BOUND - gains) >= 0
            end
            % If dynamic metric parameters
            if flagDynamicMetric
                % Positive definite constraint
                METRIC_SCALE_FACTOR*M_nonnatural + dMnn_dt >= EIG_TOL*eye(3) % Metric tensor must be PD
                % Bound max change rates
%                 abs(m_dot) <= 10
            end
            
%             if norm(rot_log(RRef,eye(3))) < 1e-3
%                 k_dot(3) == 0;
%             end
    cvx_end
    % if no solution, set rates to zero
    if  ~contains(lower(cvx_status),'solved','IgnoreCase',true)
        fprintf('Infeasible...\n')
        % Otherwise, keep the latest gains
        k_dot = zeros(3,1);
        delta = 0;
        m_dot = zeros(5,1);
        dMnn_dt = [m_dot(1), m_dot(2), m_dot(5);...
            m_dot(2), m_dot(3), m_dot(4);...
            m_dot(5), m_dot(4), 0];
        b_dot = 0;
    end
    if max(eig(contractionMatrix_current)) > 0
        fprintf('Positive contraction matrix: %f\n', max(eig(contractionMatrix_current)));
    %     error('contraction matrix not neg semi def.');
    end
    dparams = [k_dot;dMnn_dt(:);b_dot]
    delta
end

end