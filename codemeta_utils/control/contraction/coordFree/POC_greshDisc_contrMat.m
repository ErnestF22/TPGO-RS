% Test the gershgorin discs of the contraction matrix
% close all; clear all; clc;

% load('data_ACC2019.mat');
m1 = m_contract(1,1);m2 = m_contract(1,2); m3 = m_contract(2,2);
% m1 = randn; m2 = randn; m3 = randn;
kd = abs(kd); kv = abs(kv); % make sure we're using positive gains for this implementation
w = cnormalize(rand(3,1));
% w_t = @(t) w;%/norm(w)*t*sqrt(kd/2);
% w_t = @(t) sqrt(kd/2*(maxD^2-t^2))*w/norm(w);
w_t = @(t) t*w;
t = linspace(1e-3,4);
R_t = @(t) rot_expVec(eye(3),t*w);
symm = @(A) (A+A')/2; % return symmetric matric
figure;
for i = 1:length(t)
    % Find P_D bounds
%     theta = norm(rot3_log(R_t(t(i))));
    theta = t(i);
    Lambda = diag([1, mod(theta,pi)/2*cot(mod(theta,pi)/2), mod(theta,pi)/2*cot(mod(theta,pi)/2)]);
    P_D11 = -m2*kd*Lambda + m1*beta*eye(3);
    P_D21 = -m3*kd/2*Lambda + (m2*beta + (m1-m2*kv)/2)*eye(3);
    P_D12 = P_D21;
    P_D22 = 1/2*(m3*beta + m2 - m3*kv)*eye(3);
    P_D = [P_D11 P_D12; P_D21 P_D22];
    [A,maxPDEval(i)] = gershgorinDisc(P_D);
    
    subplot(2,2,1);
    plot(theta,sum(A,2),'rx');
    title('PD Matrix')
    hold on
    
    % Find Omega_D bounds
    alpha = norm(w_t(theta))*(m3*kv-m2)/4;
    O_D11 = m2/4*norm(w_t(theta))^2*diag([0, -1, -1]);
    O_D21 = m3/8*norm(w_t(theta))^2*diag([0, -1, -1]) + alpha*diag([0,1j,-1j]);
    O_D12 = conj(O_D21);
    O_D22 = 1/2*(m3*beta + m2 -m3*kv)*eye(3);
    O_D = [O_D11 O_D12; O_D21 O_D22];
    [B,maxODEval(i)] = gershgorinDisc(O_D);
    
    subplot(2,2,2);
    plot(theta,sum(B,2),'rx');
    title('OD Matrix');
    hold on
    plot(theta,-m2/4*theta^2+m3/8*theta^2+abs(m2-m3*kv)/4*theta,'go')
    
    C = gershgorinDisc(O_D11);
    subplot(2,2,3);
    plot(theta,sum(C,2),'rx');
    title('OD11');
    hold on
    
    D = gershgorinDisc(abs(O_D12));
    subplot(2,2,4);
    plot(theta,sum(D,2),'rx');
    title('OD12');
    hold on
    
    % Check contraction matrix
    U = @(R) R*hat3(w_t(theta));
    M = rotBundle_contractionMat(U,beta,kd,kv,[m1;m2;m3]);
    maxE(i) = max(eig(symm(M(R_t(t(i))))));
%     
%     if (maxPDEval(i) + maxODEval(i) > 1e-5)
%         fprintf('large jump here\n');
%     end
end

figure
plot(t,maxPDEval,t,maxODEval,t,maxPDEval+maxODEval,t,maxE,'LineWidth',3);
legend('max PD Eval','max OD Eval','max Eval','actual Eval');