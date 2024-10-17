function proposal_plot
%load dynamicSfm_testNoise_data
load('proposal_trajectoryVisual_interpolation_testNoise_data')
%figDir='~/Documents/svnDocuments/proposals/DynamicSfM/figures';
%figDir='~/Documents/BU/proposals/NSFCRIIDynamicSfM/figures/';
figDir='/Users/tron/Documents/BU/papers/vision/2016-cvpr-dynamicSfM/figures/';

figDim=[3.7 2.7];
figDimSmall=[3 2.7];
figDimLarge=[7 2.5];
figure(1)
plotPoints(X,'ro')
hold on
plotPoints(XEst,'bx')
plotLines(X,XEst,'k-')
plotPoints(Tsb,'r:')
plotPoints(TsbEst,'b-')
plotPoints(TsbIntegrated,'k-')
plotPoints(TsbEstAffine,'m--')
draw3dcameraFromRT(invR(RbsEst(:,:,1)),TsbEst(:,1),'scale',0.2,'flagAlpha',false)
hold off
view(75,40)
axis equal
savefigure(fullfile(figDir,'dynamicSfM_reconstruction'),'epsc',figDim)
view(0,0)
savefigure(fullfile(figDir,'dynamicSfM_reconstruction_front'),'epsc',figDimSmall)
view(90,90)
savefigure(fullfile(figDir,'dynamicSfM_reconstruction_top'),'epsc',figDim)

figure(2)
subplot(1,2,1)
plot(t,rot_logVec(Rbs(:,:,1),Rbs)',':',t,rot_logVec(Rbs(:,:,1),RbsEst)',t,rot_logVec(Rbs(:,:,1),RbsEstAffine)','--')
xlabel('Time [s]')
ylabel('log(R_f) [rad]')
subplot(1,2,2)
plot(t,Tsb',':',t,TsbEst',t,TsbEstAffine','--')
xlabel('Time [s]')
ylabel('T_f [m]')
savefigure(fullfile(figDir,'dynamicSfM_trajectory'),'epsc',figDimLarge)

figure(3)
subplot(1,2,1)
plot(t,dTsb',':',t,dTsbEst')
title('dTsb')
subplot(1,2,2)
plot(t,ddTsb',':',t(idxValid),ddTsbEst(:,idxValid)')
title('ddTsb')
