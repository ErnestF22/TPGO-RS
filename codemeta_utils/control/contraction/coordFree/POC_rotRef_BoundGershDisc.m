% POC_rotRef_BoundGershDisc.m
% Bound eigenvalues of contraction matrix on TSO3xSO3 (implemented in
% rotRef_contractionMat2) by relaxing the disc radii using the fact that
% the main block diagonals are symmetric, thus diagonalizable by
% orthonomral matrices. Then the radii of the gersh discs are givng as
% norm(x,1), where x = P1'*(off diagonal block matrix)*Q and Q=[P1,P2,P3]
% such that Q*Q'=I and Q'*(main block diagonal matrix)*Q = (some diagonal
% matrix of eigenvalues). Then norm(x,1) <= sqrt(3)*norm(x,2) <= 
% sqrt(3)*singluarValue(off diagonal block matrix)

close all; clear all; clc;

for ii = 1:1e5
    % Define parameters
    w = randn(3,1); kd = abs(randn); kv = abs(randn); kp = abs(randn); t = rand;
    beta = 0;
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
    L_R = eig(Dlog_RRef_R);
    L_RRef = eig(Dlog_RRef_I);

    % Compute contraction matrix
    M_contraction = rotRef_contractionMat2(beta,kd,kv,kp,M_nn);
    M_eval = M_contraction(R,w,RRef);
    % Symmertize M_eval
    M_eval = (M_eval+M_eval')/2;
    % NOTE: the gershgorin discs are different depending on the form of
    % M_eval so we must transform to get an accurate result
    [V1,~] = eig(Log_I_RRefTransposeR);
    [V2,~] = eig(Log_I_RRef);
    S = [V1 zeros(3,6);zeros(3) V2 zeros(3);zeros(3,6) V2];
    M_eval = S'*M_eval*S;

    % Compute bounds, for the radii assume worst case
    % Define bound on each off-diagonal block matrix
    M_21_Bound = @(L_R,L_RRef)  sqrt(3)*max(abs(-kd/2*m3*L_R + 1/2*(m1-m2*kv) + beta*m2))...
        +sqrt(3)*abs(-1/4*(m2-m5*m6/m4)*norm(w)*1i - 1/8*(m3-m5^2/m4)*norm(w)^2 + (m3-m5^2/m4)*kv/4*norm(w)*1i)...
        +sqrt(3)*abs(kd/4*(m3-m5^2/m4)*norm(Log_I_RRefTransposeR)*1i);
    M_31_Bound = @(L_R,L_RRef) sqrt(3)*max(abs(kd/2*(m2-m5)*L_R + beta*m6))...
        +sqrt(3)*max(abs(-kp*m6/2*L_RRef))...
        +sqrt(3)*abs(m6^2/(4*m4)*norm(w)*1i-m5*m6/(4*m4)*kv*norm(w)*1i)...
        +sqrt(3)*abs(m5*m6/(4*m4)*kd*norm(Log_I_RRefTransposeR)*1i)...
        +sqrt(3)*abs(-m6*kp/4*norm(Log_I_RRef)*1i);
    M_32_Bound = @(L_R,L_RRef) sqrt(3)*max(abs(m3*kd/2*L_R + beta*m5 + (m6-m5*kv)/2))...
        +max(abs(-m5*kp/2*L_RRef))...
        +sqrt(3)*abs((m5*m6/(4*m4)-m5^2/(4*m4)*kv)*norm(w)*1i)...
        +sqrt(3)*abs((m5^2/(4*m4)*kd*norm(Log_I_RRefTransposeR)*1i))...
        +sqrt(3)*abs(-m5*kp/4*norm(Log_I_RRef)*1i);

    % Row 1-3 discs
    D1 = @(L_R,L_RRef) -m2*kd*min(real(L_R)) + sqrt(3)*abs(-1/4*(m2-m5*m6/m4)*norm(w)^2 + beta*m1)... % centroid M(1,1)
        +M_21_Bound(L_R,L_RRef) + M_31_Bound(L_R,L_RRef);
    D2 = @(L_R,L_RRef) m2-m3*kv+beta*m3...  % centroid M(2,2)
        +M_21_Bound(L_R,L_RRef) + M_32_Bound(L_R,L_RRef);
    D3 = @(L_R,L_RRef) sqrt(3)*max(abs(m5*kd*L_R + beta*m4)) - m4*kp*min(real(L_RRef))... % centroid M(3,3)
        + M_31_Bound(L_R,L_RRef) + M_32_Bound(L_R,L_RRef);

    % Compare bounds vs actual discs
    [disc,~] = gershgorinDisc(M_eval);
    disc(:,3) = disc(:,1)+disc(:,2);
    disc(1:3,4) = D1(L_R,L_RRef);
    disc(4:6,4) = D2(L_R,L_RRef);
    disc(7:9,4) = D3(L_R,L_RRef);
    fprintf('discs = [centroid, radius, maxEval, bound, bound-maxEval>0]\n');
    disc(:,5) = disc(:,4)-disc(:,3)
    
    % Error checking for disc bound and eigenvalue bound
    % Need to convert to real since M_eval = S'*M_eval*S makes the matrix 
    % complex
    if any(real(disc(:,5)) < 0) || max(real(disc(:,4))) < max(real(eig(M_eval)))
        error('Bounds does not work');
    end
end