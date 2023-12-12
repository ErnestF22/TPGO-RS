function POC_BoundMats
% Test bounding of contraction metric by spliting into 3 matrices
% 1) depends on constants (including gains and metric scalars)
% 2) depends on omega
% 3) depends on R (position on SO(3))

close all;

% Generate random pos def metric gains
A = rand(2,2);
B = (A+A')/2+2*eye(2);
m = [B(1,1);B(1,2);B(2,2)];

% Specify system parameters
R = eye(3);
w=cnormalize(randn(3,1));
wR=pi*w;

% Set t parameters
maxT = 10;
minT = -10;
steps = 100;

% Plot eigen value of contraction metric
figure
subplot(2,2,1)
funPlot(@(t) f(t,R,w,m), linspace(minT,maxT,steps))
title('eVals of total Contraction Metric(w)')

subplot(2,2,2)
funPlot(@(t) f2(t,w,wR,m), linspace(minT,maxT,steps))
title('eVals of total Contraction Metric(R)')

% figure
subplot(2,2,3)
funPlot(@(t) fw(t,w,m), linspace(minT,maxT,steps))
title('eVals of omega parts')

% figure
subplot(2,2,4)
funPlot(@(t) fR(t,wR,m), linspace(minT,maxT,steps));
title('eVals of R parts');

function e=f(t,R,w,m)
% Compute eigenvalues of Contraction metric as omega changes
lambda = 1;
kd = 1;
kv = 1;
U = @(R) R*hat3(t*w);
M = rotBundle_contractionMat(U, lambda, kd, kv, m);
e=sort(eig(symm(M(R))));

function e=f2(t,w,wR,m)
% Compute eigenvalues of Contraction metric as omega changes
lambda = 1;
kd = 1;
kv = 1;
U = @(R) R*hat3(w);
M = rotBundle_contractionMat(U, lambda, kd, kv, m);
R=rot_expVec(eye(3),t*wR);
e=sort(eig(symm(M(R))));

function A=symm(A)
A=(A+A')/2;

function e=fw(t,w,m)
% Compute eigenvalues of 2) as omega changes
kv = 1;
M2 = [zeros(3) (-m(2)+m(3)*kv)/4*hat3(t*w)';(-m(2)+m(3)*kv)/4*hat3(t*w) zeros(3)];
M4 = [m(2)/4*hat3(t*w)^2 m(3)/8*hat3(t*w)^2;m(3)/8*hat3(t*w)^2 zeros(3)];
e = sort(eig(symm(M2+M4)));

function e=fR(t,wR,m)
% Compute eigenvalues of 3) as R changes
kd = 1;
R=rot_expVec(eye(3),t*wR);
M3 = [zeros(3) -m(3)*kd/4*rot3_log(R)';-m(3)*kd/4*rot3_log(R) zeros(3)];
M5 = [m(2)*kd*rot3_logDiffMat(eye(3),R) m(3)*kd/2*rot3_logDiffMat(eye(3),R)';...
    m(3)*kd/2*rot3_logDiffMat(eye(3),R) zeros(3)];
e = sort(eig(symm(M3+M5)))