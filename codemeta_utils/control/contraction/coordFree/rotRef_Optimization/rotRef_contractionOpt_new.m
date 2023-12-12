function [M_nn, flagFeasible,varargout] = rotRef_contractionOpt_new(mag_R,mag_RRef,kd,kv,kref,mag_W,beta,varargin)
% LAST EDITED: Jan. 28, 2021 by Bee Vang
% This is an improved version of rotRef_contractionOpt.m. We no longer
% require function handles to get the bounds, but instead rely on a
% generated function (see rotRef_allBounds_matrixForm.m) which is faster. In
% addition, we construct the LMI's in batches to reduce modeling overhead
% cost.
% INPUTS:
%   mag_R := positive scalar representing max distance from R to RRef
%   mag_RRef := positive scalar representing max distance from RRef to I
%   kd, kv, kref := scalar positive gains
%   magW := positive scalar representing maximum velocity magnitude
%   beta := positive scalar representing minimum convergence rate
% OUTPUTS:
%   m := [m1 m2 m6;m2 m3 m5;m6 m5 m4] metric tensor
%   flagFeasible :=  boolean indicating if the system is feasible (true)

%% Default Parameters
flagFeasible = false; % set to false by default
EIG_TOL = 1e-4; % Set the min eigenvalue for metric tensor
m4 = 1; % fix the scale of the metric on TSO3xSO3
BATCH_SIZE = 20; % Number of bounds to consider per batch (20 seems to be experimentally best)
%% Optional parameters
ivarargin=1;
len_varargin = length(varargin);
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'm4'
            ivarargin = ivarargin + 1;
            m4 = varargin{ivarargin};
        case 'eigtol'
            % Set the min eigenvalue for metric tensor
            ivarargin = ivarargin + 1;
            EIG_TOL = varargin{ivarargin};
        case 'batchsize'
            % Set the number of bounds to process per batch
            % NOTE: Larger size means more memory to store large matrix and
            % may result in less efficient overall computation times
            ivarargin = ivarargin + 1;
            BATCH_SIZE = varargin{ivarargin};
        otherwise
            % ignore the entry and increase index until next entry is
            % a string
            while ivarargin+1 <= len_varargin && ~ischar(varargin{ivarargin+1})
                ivarargin=ivarargin+1;
            end
    end
    ivarargin=ivarargin+1;
end
%% Solve the Optimization problem for the metric M_nn
% Compute the "constants" matrices, note A*x-B = bounds
A = rotRef_gen_allBounds_A(mag_R,mag_RRef,kd,kv,kref,mag_W,m4,beta);
ConstantTerms = rotRef_gen_allBounds_B(mag_R,mag_RRef,kd,kv,kref,mag_W,m4,beta); % Terms multlipled by m4 and any other constants (none)
NUM_BOUNDS = size(A,1);
% Break matrices down to proper components
LinearTerms = A(:,1:5); % Terms multipled by [m1;m2;m3;m5;m6]
QuadTerms = A(:,6:8); % Terms multipled by [m5*m6;m5^2;m6^2]

cvx_begin sdp quiet
    % Define x = [m1;m2;m3;m5;m6]
    variable X(5,1)
    variable theta_relax % use to at least return best result even if not all constraints are satisfied
    % Re-assign vars
    m1 = X(1); m2 = X(2); m3 = X(3); m5 = X(4); m6 = X(5);
    % Define the metric tensors
    M_nn = [m1 m2 m6;m2 m3 m5;m6 m5 m4];

    minimize theta_relax
    
    % Metric tensor must be PD
    M_nn >= EIG_TOL*eye(3)
    % NOTE: coverting LMI's to cvxr constraints is a bottle neck
    % ~18s
    %             tic
    %             for ii = 1:length(A_all_convex_2)
    %                 % Use identity and sqrt of Ai method
    %                 Ai = A_all_convex_2{ii};
    %                 % Add to cvx constraint but do not save
    %                 [-B_all_2(:,ii)'*X+theta_relax-C_all_2(ii), [m5,m6]*Ai';Ai*[m5;m6], eye(2)] >= 0
    %             end
    %             fprintf('Old LMI Method: ');toc

    % Construct large LMI matrix using kron
    for zz = 1:BATCH_SIZE:NUM_BOUNDS
        % Check for end of data
        if zz+BATCH_SIZE > NUM_BOUNDS + 1
            WORKING_SIZE = NUM_BOUNDS - zz;
        else
            WORKING_SIZE = BATCH_SIZE;
        end
        
        % Construct the LMI component matrices
        B_temp = sparse(WORKING_SIZE*3,6); % Linear terms multipled by [m1;m2;m3;m5;m6;1]
        iCounter = zz;
        for ii = 1:3:size(B_temp,1)
            B_temp(ii,:) = [-LinearTerms(iCounter,:), ConstantTerms(iCounter)];
            iCounter = iCounter + 1;
        end
        iCounter = zz;
        A_temp = cell(WORKING_SIZE,1);
        for ii = 1:WORKING_SIZE
            % Create positive semi-definite sqrtm(A) matrix
            Ai = [QuadTerms(iCounter,2) QuadTerms(iCounter,1)/2;QuadTerms(iCounter,1)/2 QuadTerms(iCounter,3)];
            [V,D] = eig(Ai);
            % Find new PSD Bi by removing negative eigenvalues. In the case
            % that both eigenvalues are <= 0, new eigenvalues should both be 0.
            % For the case of 1 neg. eigenvalue, the new eigenvalues should be
            % 0 and the original nonnegative one.
            Ai = V*( D + [max([0,-D(1,1)]) 0;0 max([0,-D(2,2)])] )*V';
            % Take the matrix square root
            if ~all(abs(Ai)<1e-6)
                tempSqrtM = sqrtm(Ai); % [2x2] matrix, need to make [5x5] later for use
            else
                tempSqrtM = zeros(2);
            end
             % Remove imaginary if its related to numerical errors
            if ~any(any(imag(tempSqrtM)>1e-5))
                Ai = real(tempSqrtM);
            else
                Ai = tempSqrtM;
            end
            
            % Create A matrix for LMI
            A_temp{ii} = [zeros(2,1),Ai;zeros(4,3)];
            iCounter = iCounter + 1;
        end
        % Convert to block diagonal matrix where each entry is the large
        % zero matrix with Ai
        A_temp = sparse(blkdiag(A_temp{:}));

        % Make idenity matrix to add for bottom right matrix
        E = sparse(eye(WORKING_SIZE));
        ID_LMI = kron(E,diag([0;1;1]));
        % And matrix to add theta_relax to the linear terms
        THETA_LMI = theta_relax*repmat([1;0;0],WORKING_SIZE,1);

        % Make cvx parameters then build LMI
        D=[ [m5 m6] zeros(1,4);zeros(2,6)];
        A_LMI = kron(E,D)*A_temp;
        A_LMI_Transpose = A_LMI';
        B_LMI = diag(B_temp*[X;1]+THETA_LMI);
        B_LMI + A_LMI + A_LMI_Transpose + ID_LMI >= 0
    end
cvx_end

% Additional checks
if  contains(lower(cvx_status),'solved','IgnoreCase',true)
    eM = eig(M_nn);
    if (all(eM > 1e-6) && max(rotRef_contractionMatrix_greshBound(mag_R,mag_RRef,kd,kv,kref,mag_W,beta,M_nn)) <= 1e-5)
        flagFeasible = true;
    end
end

%% Additional outputs
varargout{1} = cvx_optval; % This should be the maximum bound value (ideally < 0)
end