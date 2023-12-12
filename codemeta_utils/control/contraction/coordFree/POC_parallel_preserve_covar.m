% Test if parallel transport preserve the covariant derivative on SO3
clear all;

% Define tangent vectors at R1 and R2
x1 = randn(3,1);
x2 = randn(3,1);
R1 = rot_randn;
R2 = rot_randn;
X1 = @(R) R*hat3(x1);
X2 = @(R) R*hat3(x2);
parallelToR2 = @(R1,R2,X) rot_parallel(R1,R2,X(R1),'torotation');
parallelToR1 = @(R1,R2,X) rot_parallel(R2,R1,X(R2),'torotation');

% Compute covar at R1
covarR1 = rot_covar(X1,@(R1) parallelToR1(R1,R2,X1));
% and at R2
covarR2 = rot_covar(@(R2) parallelToR2(R1,R2,X1),X2);
covarNominal = rot_covar(X1,X2);
% Compare covar at R1
[rot_vee(R1,covarR1(R1)) rot_vee(R1,parallelToR1(R1,R2,covarR2))] % Not equal