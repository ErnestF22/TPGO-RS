function POCDynamicSfM_numericalDerivatives
resetRands()
load('sampleDynamicSfMDataset')
filterOrder=2;
filterWindow=21;
tFrames=min(t):1/1000:max(t);

xFrames=real_interp(t,xb,tFrames,'spline');
dxFrames=real_interp(t,dxb,tFrames,'spline');
ddxFrames=real_interp(t,ddxb,tFrames,'spline');
xFramesNoise=xFrames+0.005*randn(size(xFrames));

figure(1)
disp('Without noise')
showResults(tFrames,xFrames,filterOrder,filterWindow,xFrames,dxFrames,ddxFrames,t,xb,dxb,ddxb);
set(gcf,'name','Without noise')

figure(2)
disp('With noise')
showResults(tFrames,xFramesNoise,filterOrder,filterWindow,xFrames,dxFrames,ddxFrames,t,xb,dxb,ddxb);
set(gcf,'name','With noise')


function showResults(tFrames,xFramesNoise,filterOrder,filterWindow,xFrames,dxFrames,ddxFrames,t,xb,dxb,ddxb)
coordShow=2;
idxShow=3;
[xFilter,dxFilter,ddxFilter]=sgolayFilterDerivatives(mean(diff(tFrames)),xFramesNoise,filterOrder,filterWindow);
idxValid=round(filterWindow:length(tFrames)-filterWindow);
subplot(3,1,1)
plot(tFrames(idxValid),squeeze(xFilter(coordShow,idxShow,idxValid)),t,squeeze(xb(coordShow,idxShow,:)),...
    tFrames(idxValid),squeeze(xFilter(coordShow,idxShow,idxValid)));
legend('Filter','Original','Noise')
subplot(3,1,2)
plot(tFrames(idxValid),squeeze(dxFilter(coordShow,idxShow,idxValid)),t,squeeze(dxb(coordShow,idxShow,:)));
legend('Filter','Original')
subplot(3,1,3)
plot(tFrames(idxValid),squeeze(ddxFilter(coordShow,idxShow,idxValid)),t,squeeze(ddxb(coordShow,idxShow,:)));
legend('Filter','Original')
disp('RMSE: xFrames dxFrames ddxFrames')
disp([rmse(xFrames(:,:,idxValid),xFilter(:,:,idxValid)) rmse(dxFrames(:,:,idxValid),dxFilter(:,:,idxValid)) rmse(ddxFrames(:,:,idxValid),ddxFilter(:,:,idxValid))])
