%Evaluate energy associated with the 1-D regularized median problem
%function medianRegularizedEnergy(uMedian,u,lambda,uReg,lambdaReg)
function f=medianRegularizedEnergy(uMedian,u,lambda,uReg,lambdaReg)
f=lambdaReg*(uReg-uMedian)^2+lambda*norm(u-uMedian,1);