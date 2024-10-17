%Solves the 1-D regularized median problem
%function uMedian=medianRegularized(uhat,u)
%Finds the minimizer of lambda*norm(uHat-uMedian,1)+lambdaReg*(uReg-uMedian)^2.
%The length of uHat is assumed to be even. The correctness is proved in
%
% Y. Li, S. Osher. 
% A new median formula with applications to PDE based denoising.
% Commun. Math. Sci. 7 (2009), no. 3, 741--753
%
function uMedian=medianRegularized(u,lambda,uReg,lambdaReg)
NUHat=length(u);
uData=uReg+lambda/lambdaReg*(-NUHat/2:NUHat/2);
uMedian=median([uData u]);
