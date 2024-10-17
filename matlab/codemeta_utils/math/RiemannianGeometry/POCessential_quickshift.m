function POCessential_quickshift
%resetRands();
N=100;
noiseSigma=1;
kernelSigma=noiseSigma;
Q01=essential_eye();
Q02=essential_randn();
Q=cat(3,essential_randn(Q01,1,N),essential_randn(Q02,200,N/2));

figure(1)
disp('Computing kernel with full distance matrix')
subplot(1,2,1)
[d,f,idxMax]=test(Q,Q01,Q02,kernelSigma,'full');
subplot(1,2,2)
disp('Computing kernel with reduced distance matrix')
test(Q,Q01,Q02,kernelSigma,'sparse');

figure(2)
t=essential_kernelQuickShiftTree(d,f);
essential_kernelQuickShiftTreeSpy(t)

figure(3)
cumDistPerc(d(idxMax,:))

figure(4)
[~,idxSort]=sort(d(idxMax,:));
imagesc(d(idxSort,idxSort))

function [d,f,idxMax]=test(Q,Q01,Q02,kernelSigma,method)
tic
[QMax,d,f,idxMax]=essential_kernelmode(Q,kernelSigma,'methodDistance',method);
toc
disp('Distance [QMax-Q01 QMax-Q02]')
disp([essential_dist(QMax,Q01) essential_dist(QMax,Q02)]);
imagesc(d)
