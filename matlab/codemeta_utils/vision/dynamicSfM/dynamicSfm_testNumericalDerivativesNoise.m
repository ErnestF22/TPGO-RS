function dynamicSfm_testNumericalDerivativesNoise
%load('sampleDynamicSfMDataset_sampling')
load('sampleDynamicSfMDataset_interpolation')
%load('sampleDynamicSfMDataset_original')

fprintf('Frequency of measurements: %.1f Hz\n\n',1/mean(diff(t)))

sigmaNoise=0.001;
filterOrder=2;
filterWindow=41;
TSampling=mean(diff(t));
xbNoise=xb+sigmaNoise*randn(size(xb));
[xbFilter,dxbFilter,ddxbFilter]=sgolayFilterDerivatives(TSampling,xbNoise,filterOrder,filterWindow);

idxValid=filterWindow:length(t)-filterWindow;
subsample idxValid xb xbFilter dxb dxbFilter ddxb ddxbFilter ...
    Gammab wIMU alphaIMU Rbs Tsb taub nub MTruth t

disp('RMSE xbFilter dxbFilter ddxbFilter')
disp([rmse(xbFilter, xb) rmse(dxbFilter, dxb) rmse(ddxbFilter, ddxb)])
[XEst,RbsEst,taubEst,nubEst,output]=dynamicSfm(xbFilter,dxbFilter,ddxbFilter,Gammab,wIMU,alphaIMU,...
    'inertiamatrix',J,'RFirstFrame',Rbs(:,:,1),'showDiagnosticInfo',...
    'translationDerivativeConstraint','sgolay',1,21,TSampling,...
    'translationDerivativeConstraintWeight',10,'flagTranslationDerivativeConstraint',true,...
    'translationSecondDerivativeConstraint','sgolay',1,3,TSampling,...
    'translationSecondDerivativeConstraintWeight',10,'flagTranslationSecondDerivativeConstraint',true...
    );

disp('RMSE X')
disp(rmse(XEst,X))
disp('RMS Distance Rbs [deg]')
disp(rmse(rot_dist(RbsEst,Rbs,'vector')*180/pi))
disp('RMSE translation')
disp(rmse(taubEst,taub))
disp('RMSE linear velocity')
disp(rmse(nubEst,nub))

TsbEst=multiprodMatVec(invR(RbsEst),taubEst);
dTsbEst=multiprodMatVec(invR(RbsEst),nubEst);
figure(1)
plotPoints(X,'ro')
hold on
plotPoints(XEst,'bx')
plotPoints(Tsb,'r:')
plotPoints(TsbEst,'b-')
draw3dcameraFromRT(invR(RbsEst(:,:,1)),TsbEst(:,1),'scale',0.2)
hold off
view(0,90)
axis equal

figure(2)
subplot(3,1,1)
plot(t,rot_logVec(Rbs(:,:,1),Rbs)',':',t,rot_logVec(Rbs(:,:,1),RbsEst)')
subplot(3,1,2)
plot(t,Tsb',':',t,TsbEst')
subplot(3,1,3)
plot(t,dTsbEst')

function subsample(nameIdxVar,varargin)
for ivarargin=1:length(varargin)
    var=varargin{ivarargin};
    A=evalin('caller',var);
    sz=size(A);
    if sz(2)==1
        sz(2)=[];
    end
    switch length(sz)
        case 1
            evalin('caller',[var '=' var '(' nameIdxVar ');']);
        case 2
            evalin('caller',[var '=' var '(:,' nameIdxVar ');']);
        case 3
            evalin('caller',[var '=' var '(:,:,' nameIdxVar ');']);
    end
end
