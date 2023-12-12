% Check if the combined omega and R matrics [12x12] has greatest
% eigenvalues when log(R) and hat3(w) are colinear
close all; 

% Define constants
ANG_TOL = 1e-2;
A = rand(2,2);
B = (A+A')/2+2*eye(2);
m = [B(1,1);B(1,2);B(2,2)]
% m = [1;0;1];
% Gains and contraction rate
lambda = 1;
kd = 10*rand;
kv = 10*rand;
% Choose an arb dir for omega
w = cnormalize(randn(3,1));
t = linspace(0+ANG_TOL,pi-ANG_TOL);
% Define the omega and R matrices
R = @(t) rot_expVec(eye(3),t*w); % Find an R along the w dir
M2 = @(t) [zeros(3) (-m(2)+m(3)*kv)/4*hat3(t*w)';(-m(2)+m(3)*kv)/4*hat3(t*w) zeros(3)];
M4 = @(t) [m(2)/4*hat3(t*w)^2 m(3)/8*hat3(t*w)^2;m(3)/8*hat3(t*w)^2 zeros(3)];
M3 = @(t) [zeros(3) -m(3)*kd/4*rot3_log(R(t))';-m(3)*kd/4*rot3_log(R(t)) zeros(3)];
M5 = @(t) [m(2)*kd*rot3_logDiffMat(eye(3),R(t)) m(3)*kd/2*rot3_logDiffMat(eye(3),R(t))';...
    m(3)*kd/2*rot3_logDiffMat(eye(3),R(t)) zeros(3)];


h=figure('units','normalized','outerposition',[0 0 1 1])
eVal_max = @(t) max(eig(symm( M2(t)+M3(t)+M4(t)+M5(t) )));
funPlot(eVal_max, t);
hold on
Rrev = @(t) rot_expVec(eye(3),-t*w); % Find an R along the -w dir
M3rev = @(t) [zeros(3) -m(3)*kd/4*rot3_log(Rrev(t))';-m(3)*kd/4*rot3_log(Rrev(t)) zeros(3)];
M5rev = @(t) [m(2)*kd*rot3_logDiffMat(eye(3),Rrev(t)) m(3)*kd/2*rot3_logDiffMat(eye(3),Rrev(t))';...
    m(3)*kd/2*rot3_logDiffMat(eye(3),Rrev(t)) zeros(3)];
eVal_maxrev = @(t) max(eig(symm( M2(t)+M3rev(t)+M4(t)+M5rev(t) )));
funPlot(eVal_maxrev, t);
legend('colinear same dir', 'colinear opp dir');
set(findall(h, 'Type', 'Line'),'LineStyle', '--');
for i = 1:10
    angle = 0;
    % find a new direction
    while abs(angle) < ANG_TOL
        wTest = cnormalize(randn(3,1));
        angle = atan2(norm(cross(w,wTest)), dot(w,wTest));
    end
    
    % Define R matrices for some arbitary R
    RTest = @(t) rot_expVec(eye(3),t*wTest); % Find an R along another dir
    M3Test = @(t) [zeros(3) -m(3)*kd/4*rot3_log(RTest(t))';-m(3)*kd/4*rot3_log(RTest(t)) zeros(3)];
    M5Test = @(t) [m(2)*kd*rot3_logDiffMat(eye(3),RTest(t)) m(3)*kd/2*rot3_logDiffMat(eye(3),RTest(t))';...
        m(3)*kd/2*rot3_logDiffMat(eye(3),RTest(t)) zeros(3)];
    eValTest_max = @(t) max(eig(symm( M2(t)+M3Test(t)+M4(t)+M5Test(t) )));
        
    funPlot(eValTest_max,t);
end
% Make lines thicker
set(findall(h, 'Type', 'Line'),'LineWidth', 1);

function A = symm(A)
    A=(A+A')/2;
end