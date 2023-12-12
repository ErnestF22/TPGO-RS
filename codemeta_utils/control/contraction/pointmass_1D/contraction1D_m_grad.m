function [ m ] = contraction1D_m_grad( kv, kd, beta )
%% TEST WITH THIS
% m =contraction1D_m_grad(1,1,.0001);
% contraction1D_CheckResults(1,1,.0001,[1;m(:,end)])
%%

%For given kv, kd, beta find m such that m'Qm is minimized by gradient
%descent
%f = 1/2*m'Qm
%grad(f) = Qm
%m_bar(k) = -eps*Q*m(k-1)
%min norm(m_bar - m)^2
%such that
%   m'C1 < 0
%   m'C2 < 0
%   m'C3 < 0
%   abs(m) < alpha (bound in a box)
%Output: flag - a boolean indicating if det(G) > 0
%Output: m - the matrix [m12; m22], m11 is fixed = 1
close all; clc;

%define limits/parameters
MAX_ITER = 100;
alpha = 10;
eps = 0.001;
ERR_TOL = 1e-5;

%comparing m'Qm+m'C4 - 1 to det(G), we get the following
a = -beta^2+2*beta*kv-4*kd-kv^2;
b = kd*kv;
c = -kd^2;
d = 2*kv;
e = beta^2-2*beta*kv+2*kd;

%redefine Q, C using the new parameters
Q = [a b;b c]; % a symmetric matrix
C4 = [d;e];

%define C1
C1 = [-2*kd; 0];

%define C2
C2 = [2;beta-2*kv];

%define C3
C3 = [2-2*kd;beta-2*kv];

%initial guess for m
%m = rand(3,1);
[S,D]=eig(Q);
d=diag(D);
idx=find(d>0,1);
if ~isempty(idx)
    m=S(:,idx);
else
    m=rand(2,1);
end

f(1) = 1/2*(m(:,1)'*Q*m(:,1) + m(:,1)'*C4 - 1);

%grad ascent
for i = 1:MAX_ITER
    gradf(:,i) = Q*m(:,i)+C4/2;
    m_bar = m(:,i) + eps*gradf(:,i);   
    
    %conditions
    if (m_bar'*C1 + beta >= 0 || m_bar'*C2 >= 0 || m_bar'*C3 + beta >= 0)% || any(abs(m_bar>alpha) )
        %find mnew s.t. it is within our bounds since m_bar is not
        cvx_begin
            variable mnew(2)
            minimize ( norm(mnew-m_bar,2) )
            subject to
                %constaints for G
				mnew'*C1 + beta<=0      %G(1,1)
				mnew'*C2<=0             %G(2,2)
				mnew'*C3 + beta<=0      %trace(G)
                %M must be pos def, so if -M is neg def, we have the
                %following
                [0 1]*mnew(2) >= 0 %trace(M)
%                 mnew(2) >= mnew(1)^2 %det(M)>0
                
        cvx_end
        %adjust m(i+1) using mnew
        m(:,i+1)=mnew;
    else
        %adjust m(i+1) using m_bar
        m(:,i+1)=m_bar;
    end
    
    %update cost
    f(:,i+1) = 1/2*(m(:,i+1)'*Q*m(:,i+1) + m(:,i+1)'*C4 - 1);
    
    if (all(abs(m(:,i) - m(:,i+1)) < ERR_TOL))
        break; %converged
    end
end

%print/plot results
fprintf('Iterations: %d\n',i);

figure
plot(f);
title('Cost f = 1/2*(m''*Q*m + m''*C4 - 1)');

figure
plot(m(1,:),m(2,:));
title('m Values');

disp('initial cost')
disp(f(1))

if all(d<0)
    disp('all evals of Q are negative')
end

%return the last m estimate
m_last = m(:,end);
if (m_last'*C1+beta < 0)
    G11_Condition = 'Passed';
else
    G11_Condition = 'Failed';
end
fprintf('m''C1+beta (G11): %f, %s\n', m_last'*C1+beta, G11_Condition);
if (m_last'*C2 < 0)
    G22_Condition = 'Passed';
else
    G22_Condition = 'Failed';
end
fprintf('m''C2 (G22): %f, %s\n', m_last'*C2, G22_Condition);
if (m_last'*C3+beta < 0)
    TraG_Condition = 'Passed';
else
    TraG_Condition = 'Failed';
end
fprintf('m''C3+beta (tra(G)): %f, %s\n', m_last'*C3+beta, TraG_Condition);
if (1/2*(m_last'*Q*m_last+m_last'*C4 - 1) > 0)
    DetG_Condition = 'Passed';
else
    DetG_Condition = 'Failed';
end
fprintf('1/2*(m''*Q*m+m''*C4 - 1) (det(G)): %f, %s\n', 1/2*(m_last'*Q*m_last+m_last'*C4 - 1), DetG_Condition);
M = [1 m_last(1); m_last(1) m_last(2)]
if (trace(M) > 0)
    TraM_Condition = 'Passed';
else
    TraM_Condition = 'Failed';
end
fprintf('tra(M): %f, %s\n', trace(M), TraM_Condition);
if (det(M) > 0)
    DetM_Condition = 'Passed';
else
    DetM_Condition = 'Failed';
end
fprintf('det(M): %f, %s\n', det(M), DetM_Condition);

end

