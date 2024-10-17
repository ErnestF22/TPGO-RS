function ikose_test
NInliers=50;
NOutliers=200;
r=abs([0.3*randn(1,NInliers) 10*rand(1,NOutliers)]);
%r=abs(randn(1,NInliers));
optsIkose={[],'factorSigmaInliers',2.5};
s=ikose(r,optsIkose{:},'showStats');
disp(s)

subplot(2,1,1)
hist(r,20)
subplot(2,1,2)
ikose_plot(r,optsIkose{:})
hold off





