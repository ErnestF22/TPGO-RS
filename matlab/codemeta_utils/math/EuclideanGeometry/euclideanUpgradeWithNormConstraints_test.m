function euclideanUpgradeWithNormConstraints_test
resetRands()

desiredNorm=9.8;
S=randn(3,4);
A=desiredNorm*squeeze(sphere_randn([1;0;0],[],30));
A=A+0.1*randn(size(A));
SA=S*homogeneous(A,4);

disp('Norms of vectors in A before rectification')
disp(cnorm(SA))

[invS,AEst]=euclideanUpgradeWithNormConstraints(A,desiredNorm);

disp('Norms of vectors in A after rectification')
disp(cnorm(AEst))

disp('Variance in norm after rectification due to noise')
disp(sqrt(var(cnorm(AEst))))

%try again with homogeneous coordinates
[invS,AEst]=euclideanUpgradeWithNormConstraints(homogeneous(A,4),desiredNorm);
disp(size(AEst))

disp('Norms of vectors in A after rectification')
disp(cnorm(homogeneous(AEst,3)))


