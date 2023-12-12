% Test bounding of the differential of the log function

close all; clear all; clc;

%% Test for the max/min value of the differential
if false
    for i = 1:1e5
        R = rot_randn; 
        eVal = eig((rot3_logDiffMat(eye(3),R)+rot3_logDiffMat(eye(3),R)')/2);
        maxMin(i,:) = [max(eVal), min(eVal)];
    end

    max(maxMin,[],1)
    min(maxMin,[],1)
end

%% Test if bounding matrix is correct
% Gen pos def m matrix
A = rand(2,2);
B = (A+A')/2+2*eye(2);
m = [B(1,1);B(1,2);B(2,2)];

% Specify system parameters
lambda = 1;
kd = 1;
kv = 1;

for i = 1:1e5
    w = randn(3,1);
    R = rot_randn;
    
    Dlog = rot3_logDiffMat(eye(3),R);
    M = [m(2)*kd*Dlog m(3)*kd/2*Dlog'; m(3)*kd*Dlog zeros(3)];
    MB = @(eVal) [m(2)*kd*eVal*eye(3) m(3)*kd/2*eVal*eye(3); m(3)*kd/2*eVal*eye(3) zeros(3)];
    
    M_sym = (M+M')/2;
    Test_eVal = eig(M_sym - MB(pi/2));
    
    if any(Test_eVal > 1e-5)
        Test_eVal
        error('Bound Mat <= M_symm');
    end
end

fprintf('Done!\n');