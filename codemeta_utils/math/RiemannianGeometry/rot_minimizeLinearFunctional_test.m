function rot_minimizeLinearFunctional_test
A=randn(20,9);
RTruth=rot_randn();
b=A*RTruth(:);
REst=rot_minimizeLinearFunctional(A,b,eye(3));
disp([RTruth REst])