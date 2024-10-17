function POCDynamicSFM_noisy
resetRands()
load('sampleDynamicSfMDataset')
flagNumericalDerivatives=false;
methodDecimation='sampling';

filterOrder=4;
filterWindow=21;
sigmaNoise=0.003;%0.005;

switch methodDecimation
    case 'sampling'
        idxSampling=[1:length(t)-9800];
        tFrames=t(idxSampling);
        xbFrames=xb(:,:,idxSampling);
        dxbFrames=dxb(:,:,idxSampling);
        ddxbFrames=ddxb(:,:,idxSampling);
        GammabFrames=Gammab(:,idxSampling);
        wIMUFrames=wIMU(:,idxSampling);
        alphaIMUFrames=alphaIMU(:,idxSampling);
    case 'interpolation'
        tFrames=min(t):1/500:max(t);
        xbFrames=real_interp(t,xb,tFrames,'spline');
        dxbFrames=real_interp(t,dxb,tFrames,'spline');
        ddxbFrames=real_interp(t,ddxb,tFrames,'spline');
        GammabFrames=real_interp(t,Gammab,tFrames,'spline');
        wIMUFrames=real_interp(t,wIMU,tFrames,'spline');
        alphaIMUFrames=real_interp(t,alphaIMU,tFrames,'spline');
end

if flagNumericalDerivatives
    %xbFramesNoise=xbFrames+sigmaNoise*randn(size(xbFrames));
    xbFramesNoise=xbFrames;
    [xbFilter,dxbFilter,ddxbFilter]=sgolayFilterDerivatives(mean(diff(tFrames)),xbFramesNoise,filterOrder,filterWindow);
    idxValid=filterWindow:length(tFrames)-filterWindow;
else
    xbFilter=xbFrames+sigmaNoise*randn(size(xbFrames));
    dxbFilter=dxbFrames+2*sigmaNoise*randn(size(xbFrames));
    ddxbFilter=ddxbFrames+4*sigmaNoise*randn(size(xbFrames));
    idxValid=1:length(tFrames);
end

disp('RMSE [xbFilter dxbFilter ddxbFilter]')
disp([rmse(xbFilter(:,:,idxValid), xbFrames(:,:,idxValid)) rmse(dxbFilter(:,:,idxValid), dxbFrames(:,:,idxValid)) rmse(ddxbFilter(:,:,idxValid), ddxbFrames(:,:,idxValid))])

[XEst,RbsEst,taubEst,nubEst,output]=dynamicSfm(xbFilter(:,:,idxValid),dxbFilter(:,:,idxValid),ddxbFilter(:,:,idxValid),...
    Gammab(:,idxValid),wIMU(:,idxValid),alphaIMU(:,idxValid),...
    'inertiamatrix',J,'bodyCameraRotation',Rbc,'showDiagnosticInfo','t',tFrames(idxValid));

figure(1)
plotPoints(X,'x')
hold on
plotPoints(XEst,'o')
hold off
axis equal

disp('RMSE X')
disp(rmse(X,XEst))
