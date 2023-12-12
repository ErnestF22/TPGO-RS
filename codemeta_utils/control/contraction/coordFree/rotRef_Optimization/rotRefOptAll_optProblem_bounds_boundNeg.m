function [dparams, delta] = rotRefOptAll_optProblem_bounds_boundNeg(x,varargin)
% Last Edited: Jan. 7, 2021 by Bee Vang
% Solve optimization problem using the gershgorin bounds to find new gains,
% metric, and convergence rate
% NOTE: Only use feasible solutions that also have negative bounds to
% ensure stabilty proof
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
GAIN_UPPER_BOUND = 110;
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
[R,w,RRef,kd,kv,kref,M_nn,b] = rotDynOptAll_stateUnpack(x);
mag_R = norm(rot_log(R,RRef)); mag_RRef = norm(rot_log(RRef,eye(3)));
mag_W = norm(w);
gains = [kd;kv;kref];

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
    contractionMat_handle = rotRefOptAll_contractionMat('sym');
    contractionMatrix_current = contractionMat_handle(R,w,RRef,gains,M_nonnatural,0);
    cvx_begin sdp quiet
        variable b nonnegative
        
        maximize (b)
        
        subject to
            contractionMatrix_current + kron(b*M_nn,eye(3)) <= 0
    cvx_end
    % if no solution, set rates to zero
    if  ~contains(lower(cvx_status),'solved','IgnoreCase',true)
        fprintf('Infeasible\n')
        b = 0;
    end
    
    dparams = b;
else
    % Compute the current contraction matrix (including beta)
    cvx_begin sdp quiet
        if ~flagDynamicGains
            gains_dot = zeros(3,1);
            delta = 0;
        else
            variable gains_dot(3,1) % gain rates
            variable delta % CLF relaxation factor (positive means larger control effort, negative means less control effort)
        end
        
        if ~flagDynamicMetric
            m_dot = zeros(5,1);
        else
            variable m_dot(5,1) % [m1_dot;m2_dot;m3_dot;m5_dot;m6_dot], assume m4 is constant
        end
        variable b_dot % convergence change rate

        % Define the positive definite metric parameters rate of change matrix
        Mnn_dot = [m_dot(1), m_dot(2), m_dot(5);...
            m_dot(2), m_dot(3), m_dot(4);...
            m_dot(5), m_dot(4), 0];

%         optVars = [k_dot;m_dot;b_dot;delta];
%         Q = zeros(10); % Positive Semi-definite cost
%         P = [0;0;0;0;0;0;0;0;1;-1];
%         maximize ( optVars'*Q*optVars + P'*optVars ) 
%         maximize ( b_dot ) % looking to increase convergence rate
%         minimize ( 2*[kd*rot_vee(R,rot_log(R,RRef))'*rot_vee(R,rot_log(R,RRef)) - kv*w'*rot_vee(R,rot_log(R,RRef)),...
%                     -kd*rot_vee(R,rot_log(R,RRef))'*w + kv*(w'*w),...
%                     0]*gains_dot - b_dot )
        maximize ( 10*b_dot - delta)
%         maximize ( b_dot - delta*1e-5 - 2*k_dot'*k_dot) % looking to increase convergence rate

        subject to
            % Define the contraction constraints CBFs
            % d(bound)\d(R,w,RRef) + d(bound)\dk + d(bound)dm +
            % d(bound)\db <= 0 (a ZBF) where we use linear class K function
            allBounds_func = @(t) rotRefOpt_gen_allBounds( norm(rot_log(R_t(t),RRef_t(t))),...
                norm(rot_log(RRef_t(t),eye(3))),... % mag_RRef
                kd,kv,kref,.... % kref
                norm(w_t(t)),.... % mag_W
                M_nn(1,1),... % m1
                M_nn(1,2),... % m2
                M_nn(2,2),... % m3
                M_nn(3,3),... % m4
                M_nn(2,3),... % m5
                M_nn(1,3),... % m6
                b);
            boundVals = allBounds_func(0);
            funApproxDer(allBounds_func,0)...
                + rotRefOpt_gen_allBounds_der_A_mat(mag_R,...
                mag_RRef,... % mag_RRef
                kd,kv,kref,...
                mag_W,.... % mag_W
                M_nn(1,1),... % m1
                M_nn(1,2),... % m2
                M_nn(2,2),... % m3
                M_nn(3,3),... % m4
                M_nn(2,3),... % m5
                M_nn(1,3),... % m6
                b)*[gains_dot;m_dot(1);m_dot(2);m_dot(3);0;m_dot(4);m_dot(5);b_dot]...
                - rotRefOpt_gen_allBounds_der_B_mat(mag_R,...
                mag_RRef,... % mag_RRef
                kd,kv,kref,...
                mag_W,.... % mag_W
                M_nn(1,1),... % m1
                M_nn(1,2),... % m2
                M_nn(2,2),... % m3
                M_nn(3,3),... % m4
                M_nn(2,3),... % m5
                M_nn(1,3),... % m6
                b)...
                + boundVals <=0
            
%             % CBF to prevent beta from going below 0
            b_dot + BETA_CBF_SCALE*(b - BETA_LOWER_BOUND) >= 0
% %             abs(b_dot) <= 1e4
%             % If dynamic gains
            if flagDynamicGains
                % CLF for min norm controller
                2*[kd*rot_vee(R,rot_log(R,RRef))'*rot_vee(R,rot_log(R,RRef)) - kv*w'*rot_vee(R,rot_log(R,RRef)),...
                    -kd*rot_vee(R,rot_log(R,RRef))'*w + kv*(w'*w),...
                    0]*gains_dot...
                    + funApproxDer(@(t) norm(kd*rot_log(R_t(t),RRef_t(t))-kv*R_t(t)*hat3(w_t(t)))^2,0) <= delta

                % CBF lower bound == 0
                gains_dot + GAINS_CBF_SCALE*(gains - GAIN_LOWER_BOUND) >= 0
%                 % Bound rate of change
                abs(gains_dot) <= 30 % This is already bounded by CBF
%                 % CBF upper bound
                -gains_dot + GAINS_CBF_SCALE*(GAIN_UPPER_BOUND - gains) >= 0
            end
            % If dynamic metric parameters
            if flagDynamicMetric
                % Positive definite constraint
                METRIC_SCALE_FACTOR*M_nn + Mnn_dot >= EIG_TOL*eye(3) % Metric tensor must be PD
                % Bound max change rates
%                 abs(m_dot) <= 10
            end
    cvx_end
    % if no solution, set rates to zero
    maxBound = max(boundVals);
    if  ~contains(lower(cvx_status),'solved','IgnoreCase',true) || maxBound > 0
        fprintf('Infeasible...\n')
        % Otherwise, keep the latest gains
        gains_dot = zeros(3,1);
        delta = 0;
        m_dot = zeros(5,1);
        Mnn_dot = [m_dot(1), m_dot(2), m_dot(5);...
            m_dot(2), m_dot(3), m_dot(4);...
            m_dot(5), m_dot(4), 0];
        b_dot = 0;
    end
    contractionMat_handle = rotRefOptAll_contractionMat('sym');
    contractionMatrix_current = contractionMat_handle(R,w,RRef,gains,M_nn,b);
    if max(eig(contractionMatrix_current)) > 0
        fprintf('Positive contraction matrix: %f\n', max(eig(contractionMatrix_current)));
    %     error('contraction matrix not neg semi def.');
    end
    dparams = [gains_dot;Mnn_dot(:);b_dot]
    fprintf('delta: %f\n', delta)
    fprintf('max bound: %f\n', maxBound)
    fprintf('Optimal Val: %f\n',cvx_optval)
%     AA=2*[kd*rot_vee(R,rot_log(R,RRef))'*rot_vee(R,rot_log(R,RRef)) - kv*w'*rot_vee(R,rot_log(R,RRef)),...
%         -kd*rot_vee(R,rot_log(R,RRef))'*w + kv*(w'*w),...
%         0]*gains_dot...
%         + funApproxDer(@(t) norm(kd*rot_log(R_t(t),RRef_t(t))-kv*R_t(t)*hat3(w_t(t)))^2,0);
%     fprintf('Actual total rate of change: %f\n', AA);
end

end