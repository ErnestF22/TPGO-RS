function [M_nn, flagFeasible,varargout] = rotRef_contractionOpt(mag_R,mag_RRef,kd,kv,kp,mag_W,beta,varargin)
% Bound eigenvalues of contraction matrix on TSO3xSO3 (implemented in
% rotRef_contractionMat2) by relaxing the disc radii using the fact that
% the main block diagonals are symmetric, thus diagonalizable by
% orthonomral matrices. Then the radii of the gersh discs are givng as
% norm(x,1), where x = P1'*(off diagonal block matrix)*Q and Q=[P1,P2,P3]
% such that Q*Q'=I and Q'*(main block diagonal matrix)*Q = (some diagonal
% matrix of eigenvalues). Then norm(x,1) <= sqrt(3)*norm(x,2) <= 
% sqrt(3)*singluarValue(off diagonal block matrix)
% NOTE: Method used is Semi-definite relaxation (same as S-procedure where
% we view X = x*x' and x = [1;m1;m2;m3;m4;m5;m6]. If the solution to X is
% optimal then rank(X) == 1, otherwise we'll need to approximate a X_approx
% with rank 1 close to X. One method is the SVD where x_approx =
% eig_max(X)*v_max where v_max is the eigenvector associated with
% eig_max(X).
% INPUTS:
%   mag_R := positive scalar representing max distance from R to RRef
%   mag_RRef := positive scalar representing max distance from RRef to I
%   kd, kv, kp := scalar positive gains
%   magW := positive scalar representing maximum velocity magnitude
%   beta := positive scalar representing minimum convergence rate
% OUTPUTS:
%   m := [m1 m2 m6;m2 m3 m5;m6 m5 m4] metric tensor
%   flagFeasible :=  boolean indicating if the system is feasible (true)

%% Default Parameters
flagFeasible = false; % set to false by default
OPT_FORM = 'sdp_relaxation'; % Solve using the nonconvex equations
EIG_TOL = 1e-4; % Set the min eigenvalue for metric tensor
m4 = 1; % fix the scale of the metric on TSO3xSO3
%optional parameters
ivarargin=1;
len_varargin = length(varargin);
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'm4'
            ivarargin = ivarargin + 1;
            m4 = varargin{ivarargin};
        case 'sdp_relaxation_convex'
            % Solve by first converting all function handles of Ai matrix
            % in trace(Ai*X) <= 0 such that the quadratic matrix Bi \in Ai
            % is positive semi-definite. Bi = Ai(2:end,2:end);
            OPT_FORM = varargin{ivarargin};
            ivarargin = ivarargin + 1;
            A_all_func = varargin{ivarargin};
        case 'schur_comp'
            % Solve by first converting all function handles of Ai matrix
            % in trace(Ai*X) <= 0 such that the quadratic matrix Bi \in Ai
            % is positive semi-definite. Bi = Ai(2:end,2:end);
            % Then the convex constraint is put into LMI form via Schur
            % complement lemma
            OPT_FORM = varargin{ivarargin};
            ivarargin = ivarargin + 1;
            A_all_func = varargin{ivarargin};
        case 'eigtol'
            % Set the min eigenvalue for metric tensor
            ivarargin = ivarargin + 1;
            EIG_TOL = varargin{ivarargin};
        otherwise
            % ignore the entry and increase index until next entry is
            % a string
            while ivarargin+1 <= len_varargin && ~ischar(varargin{ivarargin+1})
                ivarargin=ivarargin+1;
            end 
    end
    ivarargin=ivarargin+1;
end
%% Optimize 
switch OPT_FORM
    case 'sdp_relaxation_convex'
        % Solve by converting nonconvex constraints to be convex then setup
        % sdp relaxation where X=x'*x and x=[m1;m2;m3;m5;m6]
        A_all_convex = cell(size(A_all_func));
        for ii = 1:length(A_all_func)
            tempA_func = A_all_func{ii};
            tempA_eVal = tempA_func(mag_R,mag_RRef,kd,kv,kp,mag_W,beta,m4);
            % extract the Bi matrix in [m5 m6]*Bi*[m5;m6] matrix. These are the
            % only non-zero terms of the quadratic matrix
            Ai = tempA_eVal(5:6,5:6);
            [V,D] = eig(Ai);
            % Find new PSD Bi by removing negative eigenvalues. In the case
            % that both eigenvalues are <= 0, new eigenvalues should both be 0.
            % For the case of 1 neg. eigenvalue, the new eigenvalues should be
            % 0 and the original nonnegative one.
            Ai = V*( D + [max([0,-D(1,1)]) 0;0 max([0,-D(2,2)])] )*V';
            tempA_eVal(5:6,5:6) = Ai;
            A_all_convex{ii} = tempA_eVal;
        end

        cvx_begin sdp quiet
            % Define x = [1;m1;m2;m3;m5;m6] then X=x*x'
            variable X(6,6) semidefinite
            % Re-assign vars
            m1 = X(1,2); m2 = X(1,3); m3 = X(1,4); m5 = X(1,5); m6 = X(1,6);
            % Define the metric tensors
            M_nn = [m1 m2 m6;m2 m3 m5;m6 m5 m4];
            % Define gershgorin bounds
            set_A = cvx(zeros(size(A_all_convex))); % Must be cvx data
            for ii = 1:length(A_all_convex)
                set_A(ii) = trace(A_all_convex{ii}*X);
            end

            % Setup convex SDP
            minimize ( trace(X) + max(set_A))
    %         maximize ( trace(m) )
            X(1,1) == 1 % Set to 1 to give us non quadratic terms
    %         max(set_A) <= 0 % Bound eigenvalues to be negative semi-definite by greshgorin discs
            M_nn >= EIG_TOL*eye(3) % Metric tensor must be PD
        cvx_end
    case 'schur_comp'
        % Solve by converting nonconvex constraints to be convex then use
        % schur complement to convert to LMI constraints
        % Given quadratic constraint x'*A_all_convex*A_all_convex*x +
        % B_all'*x + C_all <= 0
        A_all_convex = cell(size(A_all_func));
        B_all = cell(size(A_all_func));
        C_all = cell(size(A_all_func));
        for ii = 1:length(A_all_func)
            tempA_func = A_all_func{ii};
            tempA_eVal = tempA_func(mag_R,mag_RRef,kd,kv,kp,mag_W,beta,m4);
            % extract the Bi matrix in [m5 m6]*Bi*[m5;m6] matrix. These are the
            % only non-zero terms of the quadratic matrix
            Ai = tempA_eVal(5:6,5:6);
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
                A_all_convex{ii} = real(tempSqrtM);
            else
                A_all_convex{ii} = tempSqrtM;
            end


            B_all{ii} = 2*tempA_eVal(1,2:end)'; % [5x1] vector
            C_all{ii} = tempA_eVal(1,1); % [1x1] vector for factors multipled by m4
        end
        cvx_begin sdp quiet
            % Define x = [m1;m2;m3;m5;m6]
            variable X(5,1)
            variable theta_relax % use to at least return best result even if not all constraints are satisfied
%             theta_relax = 0;
            % Re-assign vars
            m1 = X(1); m2 = X(2); m3 = X(3); m5 = X(4); m6 = X(5);
            % Define the metric tensors
            M_nn = [m1 m2 m6;m2 m3 m5;m6 m5 m4];
            % Define convex LMI bounds
%             set_A = cvx(6,6,zeros(size(A_all_convex))); % Must be cvx data
%             maximize ( trace(M_nn) )
            minimize theta_relax
            % NOTE: coverting LMI's to cvxr constraints is a bottle neck
            % ~18s
            for ii = 1:length(A_all_convex)
                % Use identity and sqrt of Ai method
                Ai = A_all_convex{ii};
                % Add to cvx constraint but do not save
                [-B_all{ii}'*X+theta_relax-C_all{ii}, [m5,m6]*Ai';Ai*[m5;m6], eye(2)] >= 0
            end
            M_nn >= EIG_TOL*eye(3) % Metric tensor must be PD
        cvx_end
        set_A = -inf; % Set to nothing for now   

    otherwise % 'sdp_relaxation' solve nonconvex sdp relaxation 
    % Solve nonconvex QCQP with approximate solution
        cvx_begin sdp quiet
            % Define x = [1;m1;m2;m3;m5;m6] then X=x*x' and let m4==1
            variable X(6,6) semidefinite
        %     variable x(7,1)
        %     X = x*x';
            % Re-assign vars
            m1 = X(1,2); m2 = X(1,3); m3 = X(1,4); m5 = X(1,5); m6 = X(1,6);
            m5_squared = X(5,5); m6_squared = X(6,6); m5m6 = X(5,6);
            % Define the metric tensors
            M_nn = [m1 m2 m6;m2 m3 m5;m6 m5 m4];
            % Define the gershgorin disc bounds
            M_21_Bound = sqrt(3)*max(abs( -kd/2*m3*[1;mag_R/2*cot(mag_R/2)] + kd*m5_squared/(4*m4)*mag_R*[0;1i] + 1/2*(m1-m2*kv) + beta*m2 ))...
                +sqrt(3)*abs(1/8*(-m3+m5_squared/m4)*mag_W^2 + 1/4*(m2-m5m6/m4-m3*kv+m5_squared/m4*kv)*mag_W*1i);
            M_31_Bound = sqrt(3)*max(abs( kd/2*(m2-m5)*[1;mag_R/2*(cot(mag_R/2)+1i)] + m5m6/(4*m4)*kd*mag_R*[0;-1i] ))...
                +sqrt(3)*abs(-m6*kp/2 + beta*m6)...
                +sqrt(3)*abs( (m6_squared-m5m6*kv)/(4*m4)*mag_W );
            M_32_Bound = sqrt(3)*max(abs( m3*kd/2*[1;mag_R/2*(cot(mag_R/2)+1i)] + m5_squared/(4*m4)*kd*mag_R*[0;-1i] + (m6-m5*kv)/2 ))...
                +abs(-m5*kp/2 + beta*m5 )...
                +sqrt(3)*abs( (-m5m6+m5_squared*kv)/(4*m4)*mag_W );

            % Row 1-3 discs
            D1 = -m2*kd*mag_R/2*cot(mag_R/2) + sqrt(3)*max(abs(-1/4*(m2-m5m6/m4)*mag_W^2*[0;1] + beta*m1 ))... % centroid M(1,1)
                +M_21_Bound + M_31_Bound;
            D2 = m2-m3*kv + beta*m3...  % centroid M(2,2)
                +M_21_Bound + M_32_Bound;
            D3 = -kp*m4*mag_RRef/2*cot(mag_RRef/2) + sqrt(3)*abs(m5)*kd + beta*m4... % centroid M(3,3)
                + M_31_Bound + M_32_Bound;

            % Define set to max over
            set_A = [D1;D2;D3];

            % Setup feasibility SDP (no objective == min/max 0)
            minimize ( trace(X) )
        %     minimize ( max(set_A) )
        %     maximize ( trace(m) )
            X(1,1) == 1 % Set to 1 to give us non quadratic terms 
            max(set_A) <= 0 % Bound eigenvalues to be negative semi-definite by greshgorin discs
            M_nn >= EIG_TOL*eye(3) % Metric tensor should be PD
        cvx_end
end
% Additional checks
if  contains(lower(cvx_status),'solved','IgnoreCase',true)
    eM = eig(M_nn);
%     [kd,kv,kp,beta,max(set_A)]
    if (all(eM > 1e-6) && max(rotRef_contractionMatrix_greshBound(mag_R,mag_RRef,kd,kv,kp,mag_W,beta,M_nn)) <= 1e-5)
        flagFeasible = true;
    end
end

varargout{1} = X;
varargout{2} = cvx_optval;
varargout{3} = set_A;
end

