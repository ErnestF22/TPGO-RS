% Test jordan forms of contraction metric

close all; clear all; clc;

% Generate random pos def metric gains
A = rand(2,2);
B = (A+A')/2+2*eye(2);
m = [B(1,1);B(1,2);B(2,2)];

% Specify system parameters
R = rot_randn;
w=cnormalize(randn(3,1));
% Gains and contraction rate
lambda = 1;
kd = randn;
kv = 1;
t = 1;
U = @(R) R*hat3(t*w);

% Matrices of constants/gains
M1 = kron([lambda*m(1) lambda*m(2)+(m(1)-m(2)*kv)/2;...
    lambda*m(2)+(m(1)-m(2)*kv)/2 lambda*m(3)+m(2)-m(3)*kv], eye(3));

% Matrices depending on omega
M2 = [zeros(3) (-m(2)+m(3)*kv)/4*hat3(t*w)';(-m(2)+m(3)*kv)/4*hat3(t*w) zeros(3)];
M4 = [m(2)/4*hat3(t*w)^2 m(3)/8*hat3(t*w)^2;m(3)/8*hat3(t*w)^2 zeros(3)];
% Get eigenvectors of hat3(t*w)
% Note eigenvectors of hat3(t*w)^2 are same as hat3(t*w), however the
% transformation matrix using the squared matrix does not work for the 1st
% order matrix
[Vw,Dw] = eig(hat3(t*w));
% Convert M2+M4 to diagonal form
Sw = [Vw zeros(3); zeros(3) Vw];
Diag_w = Sw\(M2+M4)*Sw
% Find the real Jordan form
[Vw_real,Dw_real] = cdf2rdf(Vw,Dw);
Sw_real = [Vw_real zeros(3);zeros(3) Vw_real];
Diag_w_real = Sw_real\(M2+M4)*Sw_real
% As long as eigen values align within the block matrices, the whole matrix
% should give the same eigenvalues
Diag_w_analytic = kron(norm(t*w)^2*[m(2)/4 m(3)/8; m(3)/8 0],...
    [0 0 0;0 -1 0;0 0 -1]) + kron(norm(t*w)*[0 -(-m(2)+m(3)*kv)/4;(-m(2)+m(3)*kv)/4 0],...
    [0 0 0;0 +i 0;0 0 -i]);
    

% Matrices depending on R
M3 = [zeros(3) m(3)*kd/4*rot3_log(R)';m(3)*kd/4*rot3_log(R) zeros(3)];
M5 = [-m(2)*kd*rot3_logDiffMat(eye(3),R) -m(3)*kd/2*rot3_logDiffMat(eye(3),R)';...
    -m(3)*kd/2*rot3_logDiffMat(eye(3),R) zeros(3)];
% rot3_log and rot3_logDiffMat should have the same eigenvectors??
[VR,DR] = eig(rot3_log(R));
% Convert M3+M5 to diagonal form
SR = [VR zeros(3);zeros(3) VR];
Diag_R = SR\(M3+M5)*SR
% For some reason, the top right/bottom left block matrices are
% symmerical!! (rot3_logDiffMat(eye(3),R)-1/2*rot3_log(R) is symmertical!)
% Find the real Jordan form
[VR_real,DR_real] = cdf2rdf(VR,DR);
SR_real = [VR_real zeros(3);zeros(3) VR_real];
Diag_R_real = SR_real\(M3+M5)*SR_real
eDlog = eig(rot3_logDiffMat(eye(3),R));
Diag_R_analytic = [-m(2)*kd*diag(eDlog) -m(3)*kd/2*diag(real(eDlog));...
    -m(3)*kd/2*diag(real(eDlog)) zeros(3)];

% Check that M1+...+M5 = rotBundle_contractionMat
M = rotBundle_contractionMat(U, lambda, kd, kv, m);
% LogicalStr = {'true','false'};
fprintf('M1+M2+M3+M4+M5 == rotBundle_contractionMat: %s \n', ...
    mat2str(~any(any(M(R)-(M1+M2+M3+M4+M5)>1e-9))) );

% Show different eigenvalues
[eig(M1), eig(M2+M4), eig(M3+M5)]
eig(M(R))