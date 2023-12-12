% POC_rotRef_BoundGershDisc_2.m
% Bound eigenvalues of contraction matrix on TSO3xSO3 (implemented in
% rotRef_contractionMat2) by relaxing the disc radii using the fact that
% the main block diagonals are symmetric, thus diagonalizable by
% orthonomral matrices. Then the radii of the gersh discs are givng as
% norm(x,1), where x = P1'*(off diagonal block matrix)*Q and Q=[P1,P2,P3]
% such that Q*Q'=I and Q'*(main block diagonal matrix)*Q = (some diagonal
% matrix of eigenvalues). Then norm(x,1) <= sqrt(3)*norm(x,2) <= 
% sqrt(3)*singluarValue(off diagonal block matrix)
% NOTE: same as POC_rotRef_BoundGershDisc.m but combine hessisan and
% gradient

close all; clear all; clc;

for ii = 1:1e5
    % Define parameters
    w = randn(3,1); kd = abs(randn); kv = abs(randn); kp = abs(randn); t = rand;
    beta = abs(randn);
    M_nn = randn(3,3); M_nn = M_nn'*M_nn; % Nonnatural metric
    m1 = M_nn(1,1); m2 = M_nn(1,2); m3 = M_nn(2,2);
    m4 = M_nn(3,3); m5 = M_nn(2,3); m6 = M_nn(1,3);
    % For now, since we're choosing S as matrix of eigenvectors that
    % diagonalizes DLog_R_RRef and DLog_RRef_I3, we need to bound m2>=0
    m2=abs(m2); M_nn(1,2)=m2; M_nn(2,1)=m2;
    
    R = rot_randn; RRef = rot_randn;
    Dlog_RRef_R = rot3_logDiffMat(RRef,R);
    Dlog_RRef_I = rot3_logDiffMat(RRef,eye(3));
    Log_I_RRefTransposeR = rot_log(eye(3),RRef'*R);
    Log_I_RRef = rot_log(eye(3),RRef);
    theta_R = norm(Log_I_RRefTransposeR);
    theta_RRef = norm(Log_I_RRef);
    mag_W = norm(w);

    % Compute contraction matrix
    M_contraction = rotRef_contractionMat2(beta,kd,kv,kp,M_nn);
    M_eval = M_contraction(R,w,RRef);
    % Symmertize M_eval
    M_eval = (M_eval+M_eval')/2;
    % NOTE: the gershgorin discs are different depending on the form of
    % M_eval so we must transform to get an accurate result
    [V1,e_RRef_R] = eig(Log_I_RRefTransposeR);
    [V2,e_I_RRef] = eig(Log_I_RRef);
    % Compute eigenvalues as vectors
    e_RRef_R = diag(e_RRef_R);
    e_I_RRef = diag(e_I_RRef);
    L_R = V1\Dlog_RRef_R*V1; L_R = diag(L_R);
    L_RRef = V2\Dlog_RRef_I*V2; L_RRef = diag(L_RRef);
    % Compute transformed contraction matrix
    S = [V1 zeros(3,6);zeros(3) V2 zeros(3);zeros(3,6) V2];
    M_eval = S'*M_eval*S;

    % Compute bounds, for the radii assume worst case
    % Define bound on each off-diagonal block matrix
    M_21_Bound = @(L_R,L_RRef)  sqrt(3)*max(abs(-kd/2*m3*L_R + 1/2*(m1-m2*kv) + beta*m2))...
        +sqrt(3)*abs(-1/4*(m2-m5*m6/m4)*norm(w)*1i - 1/8*(m3-m5^2/m4)*norm(w)^2 + (m3-m5^2/m4)*kv/4*norm(w)*1i)...
        +sqrt(3)*abs(kd/4*(m3-m5^2/m4)*norm(Log_I_RRefTransposeR)*1i);
    M_21_Bound_Combined = @(L_R,L_RRef)  sqrt(3)*max(abs(-kd/2*m3*real(L_R) + 1/2*(m1-m2*kv) + beta*m2 -kd/4*m5^2/m4*e_RRef_R))...
        +sqrt(3)*abs(-1/4*(m2-m5*m6/m4)*norm(w)*1i - 1/8*(m3-m5^2/m4)*norm(w)^2 + (m3-m5^2/m4)*kv/4*norm(w)*1i);
    M_21_Bound_Analytical = sqrt(3)*max(abs(-kd/2*m3*[1;theta_R/2*cot(theta_R/2)]+1/2*(m1-m2*kv) + beta*m2 - kd*m5^2/(4*m4)*theta_R*[0;1i]))...
        + sqrt(3)*abs(-1/8*(m3-m5^2/m4)*mag_W^2 - 1/4*(m2-m5*m6/m4 - m3*kv + m5^2/m4*kv)*mag_W*1i);
    M_31_Bound = @(L_R,L_RRef) sqrt(3)*max(abs(kd/2*(m2-m5)*L_R + beta*m6))...
        +sqrt(3)*max(abs(-kp*m6/2*L_RRef))...
        +sqrt(3)*abs(m6^2/(4*m4)*norm(w)*1i-m5*m6/(4*m4)*kv*norm(w)*1i)...
        +sqrt(3)*abs(m5*m6/(4*m4)*kd*norm(Log_I_RRefTransposeR)*1i)...
        +sqrt(3)*abs(-m6*kp/4*norm(Log_I_RRef)*1i);
    M_31_Bound_Combined = @(L_R,L_RRef) sqrt(3)*max(abs(kd/2*(m2-m5)*L_R + beta*m6 +m5*m6/(4*m4)*kd*conj(e_RRef_R)))...
        +sqrt(3)*max(abs(-kp*m6/2*real(L_RRef)))...
        +sqrt(3)*abs(m6^2/(4*m4)*norm(w)*1i-m5*m6/(4*m4)*kv*norm(w)*1i);
    M_31_Bound_Analytical = sqrt(3)*max(abs(kd/2*(m2-m5)*[1;theta_R/2*(cot(theta_R/2)+1i)] + beta*m6 + m5*m6/(4*m4)*kd*theta_R*[0;-1i]))...
        + sqrt(3)*abs(-m6*kp/2)...
        + sqrt(3)*abs((m6^2-m5*m6*kv)/(4*m4)*mag_W);
    M_32_Bound = @(L_R,L_RRef) sqrt(3)*max(abs(m3*kd/2*L_R + beta*m5 + (m6-m5*kv)/2))...
        +max(abs(-m5*kp/2*L_RRef))...
        +sqrt(3)*abs((m5*m6/(4*m4)-m5^2/(4*m4)*kv)*norm(w)*1i)...
        +sqrt(3)*abs((m5^2/(4*m4)*kd*norm(Log_I_RRefTransposeR)*1i))...
        +sqrt(3)*abs(-m5*kp/4*norm(Log_I_RRef)*1i);
    M_32_Bound_Combined = @(L_R,L_RRef) sqrt(3)*max(abs(m3*kd/2*L_R + beta*m5 + (m6-m5*kv)/2 + m5^2/(4*m4)*kd*conj(e_RRef_R)))...
        +max(abs(-m5*kp/2*real(L_RRef)))...
        +sqrt(3)*abs((m5*m6/(4*m4)-m5^2/(4*m4)*kv)*norm(w)*1i);
    M_32_Bound_Analytical = sqrt(3)*max(abs(m3*kd/2*[1;theta_R/2*(cot(theta_R/2)+1i)] + beta*m5 + (m6-m5*kv)/2 + m5^2*kd/(4*m4)*theta_R*[0;-1i]))...
        + abs(-m5*kp/2)...
        + sqrt(3)*abs( (m5*m6-m5^2*kv)/(4*m4)*mag_W );

    % Row 1-3 discs
    D1 = @(L_R,L_RRef) -m2*kd*min(real(L_R)) + sqrt(3)*max(abs(-1/4*(m2-m5*m6/m4)*norm(w)^2*[0;1] + beta*m1)) - (1/4*(m2-m5*m6/m4)*norm(w)^2)... % centroid M(1,1)
        +M_21_Bound(L_R,L_RRef) + M_31_Bound(L_R,L_RRef);
    D2 = @(L_R,L_RRef) m2-m3*kv+beta*m3...  % centroid M(2,2)
        +M_21_Bound(L_R,L_RRef) + M_32_Bound(L_R,L_RRef);
    D3 = @(L_R,L_RRef) sqrt(3)*max(abs(m5*kd*real(L_R))) + beta*m4 - m4*kp*min(real(L_RRef))... % centroid M(3,3)
        + M_31_Bound(L_R,L_RRef) + M_32_Bound(L_R,L_RRef);

    D1_Combined = @(L_R,L_RRef) -m2*kd*min(real(L_R)) + sqrt(3)*max(abs(-1/4*(m2-m5*m6/m4)*norm(w)^2*[0;1] + beta*m1)) - (1/4*(m2-m5*m6/m4)*norm(w)^2)... % centroid M(1,1)
        +M_21_Bound_Combined(L_R,L_RRef) + M_31_Bound_Combined(L_R,L_RRef);
    D2_Combined = @(L_R,L_RRef) m2-m3*kv+beta*m3...  % centroid M(2,2)
        +M_21_Bound_Combined(L_R,L_RRef) + M_32_Bound_Combined(L_R,L_RRef);
    D3_Combined = @(L_R,L_RRef) sqrt(3)*max(abs(m5*kd*real(L_R))) + beta*m4 - m4*kp*min(real(L_RRef))... % centroid M(3,3)
        + M_31_Bound_Combined(L_R,L_RRef) + M_32_Bound_Combined(L_R,L_RRef);
    
    D1_Analytical = -m2*kd*theta_R/2*cot(theta_R/2) + sqrt(3)*max(abs(-1/4*(m2-m5*m6/m4)*mag_W^2*[0;1] + beta*m1)) - (1/4*(m2-m5*m6/m4)*norm(w)^2)... % centroid M(1,1)
        + M_21_Bound_Analytical + M_31_Bound_Analytical;
    D2_Analytical = m2 - m3*kv + beta*m3 ... % centroid M(2,2)
        + M_21_Bound_Analytical + M_32_Bound_Analytical;
    D3_Analytical = sqrt(3)*abs(m5)*kd + beta*m4 - m4*kp*theta_RRef/2*cot(theta_RRef/2)... % centroid M(3,3)
        + M_31_Bound_Analytical + M_32_Bound_Analytical;
    % Compare bounds vs actual discs
    [disc,~] = gershgorinDisc(M_eval);
    disc(:,3) = disc(:,1)+disc(:,2);
    disc(1:3,4) = D1_Analytical;
    disc(4:6,4) = D2_Analytical;
    disc(7:9,4) = D3_Analytical;
%     fprintf('discs = [centroid, radius, maxEval, bound, bound-maxEval>0]\n');
    disc(:,5) = disc(:,4)-disc(:,3);
    
    % Error checking for disc bound and eigenvalue bound
    % Need to convert to real since M_eval = S'*M_eval*S makes the matrix 
    % complex
    if max(real(disc(:,4))) < max(real(eig(M_eval))) %|| any(real(disc(:,5)) < 0) % The 2nd condition is that our bounded disc should be >= the actual numerical discs (This may be not always be true)
        error('Bounds does not work');
    end
    
    % Print difference in bound if non-combined is smaller
    BoundError = [D1(L_R,L_RRef)-D1_Combined(L_R,L_RRef);...
        D2(L_R,L_RRef)-D2_Combined(L_R,L_RRef);...
        D3(L_R,L_RRef)-D3_Combined(L_R,L_RRef)];
    if any(BoundError < 0)
        fprintf('The combined bounds is greater than non-combined terms, should not happen, see negative bound error\n');
        disc
        BoundError
    end
    
    BoundError_Analytical = [D1_Combined(L_R,L_RRef)-D1_Analytical;...
        D2_Combined(L_R,L_RRef)-D2_Analytical;...
        D3_Combined(L_R,L_RRef)-D3_Analytical];
    if any(abs(BoundError_Analytical) > 1e-6)
        fprintf('The combined bounds is greater than the analytical bounds, should not happen...\n');
        disc
        BoundError_Analytical
    end
end