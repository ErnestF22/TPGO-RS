function [max_eVals,varargout] = rotRef_contractionMatrix_greshBound(mag_R,mag_RRef,kd,kv,kp,mag_W,beta,M_nn)
% Compute the greshgorin disc bounds given by relaxing the disc radii 
% using the fact that the main block diagonals are symmetric, thus diagonalizable by
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
%   M_nn := [3x3] non-natural metric tensor on TSO3xSO3
% OUTPUTS:
%   max_eVals := max eigenvalue bounded by the 3 greshgorin discs
%   varargout := {1} centriods
%                {2} radii

% Extract the m's
m1 = M_nn(1,1); m2 = M_nn(1,2); m3 = M_nn(2,2); m4 = M_nn(3,3);
m5 = M_nn(2,3); m6 = M_nn(1,3);

% Compute the bounds
M_21_Bound = sqrt(3)*max(abs( -kd/2*m3*[1;mag_R/2*cot(mag_R/2)] - kd*m5^2/(4*m4)*mag_R*[0;1i] + 1/2*(m1-m2*kv) + beta*m2 ))...
    +sqrt(3)*abs(-1/8*(m3-m5^2/m4)*mag_W^2 - 1/4*(m2-m5*m6/m4-m3*kv+m5^2/m4*kv)*mag_W*1i);
M_31_Bound = sqrt(3)*max(abs( kd/2*(m2-m5)*[1;mag_R/2*(cot(mag_R/2)+1i)] + m5*m6/(4*m4)*kd*mag_R*[0;-1i] ))...
    +sqrt(3)*abs(-m6*kp/2 + beta*m6)...
    +sqrt(3)*abs( (m6^2-m5*m6*kv)/(4*m4)*mag_W );
M_32_Bound = sqrt(3)*max(abs( m3*kd/2*[1;mag_R/2*(cot(mag_R/2)+1i)] + m5^2/(4*m4)*kd*mag_R*[0;-1i] + (m6-m5*kv)/2 ))...
    +max(abs([-m5*kp/2+beta*m5, -m5*kp/2*mag_RRef/2*cot(mag_RRef/2)+beta*m5]) )...
    +sqrt(3)*abs( (m5*m6-m5^2*kv)/(4*m4)*mag_W );

% Row 1-3 discs
% Bound 1/4*m2_prime*hat3(w)^2 by 1/4*m2_prime*norm(w)^2*eye(3)
centroid_1_1 = -m2*kd*mag_R/2*cot(mag_R/2) + max(0*abs(1/4*(m2-m5*m6/m4)*mag_W^2*[0;1] + beta*m1 ));
D1 = centroid_1_1 + M_21_Bound + M_31_Bound;
centroid_2_2 = m2-m3*kv + beta*m3;
D2 = centroid_2_2 + M_21_Bound + M_32_Bound;
centroid_3_3 = -kp*m4*mag_RRef/2*cot(mag_RRef/2) + sqrt(3)*abs(m5)*kd + beta*m4;
D3 = centroid_3_3 + M_31_Bound + M_32_Bound;

% Results
max_eVals = [D1;D2;D3];
varargout{1} = [centroid_1_1;centroid_2_2;centroid_3_3];
varargout{2} = [M_21_Bound + M_31_Bound;...
    M_21_Bound + M_32_Bound;...
    M_31_Bound + M_32_Bound];
end

