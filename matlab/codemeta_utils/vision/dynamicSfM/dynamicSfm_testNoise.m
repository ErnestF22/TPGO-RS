function dynamicSfm_testNoise
resetRands()
%load('sampleDynamicSfMDataset_sampling')
%load('sampleDynamicSfMDataset_interpolation')
%load('sampleDynamicSfMDataset_original')
fileName='proposal_trajectoryVisual_interpolation';
load(fileName)

fprintf('Frequency of measurements: %.1f Hz\n\n',1/mean(diff(t)))

sigmaNoiseImage=0.01;
sigmaNoiseW=deg2rad(0.3);
sigmaNoiseAlpha=0.2;
%sigmaNoiseW=0.001;
%sigmaNoiseAlpha=0.01;
filterOrder=2;
filterWindow=21;
TSampling=mean(diff(t));
xbNoise=xb+sigmaNoiseImage*randn(size(xb));
[xbFilter,dxbFilter,ddxbFilter]=sgolayFilterDerivatives(TSampling,xbNoise,filterOrder,filterWindow);

wIMUNoise=wIMU+sigmaNoiseW*randn(size(wIMU));
alphaIMUNoise=alphaIMU+sigmaNoiseAlpha*randn(size(alphaIMU));

idxValid=filterWindow:length(t)-filterWindow;
subsample idxValid xb xbNoise  xbFilter dxb dxbFilter ddxb ddxbFilter ...
    Gammab wIMUNoise alphaIMUNoise Rbs Tsb dTsb ddTsb taub nub MTruth t

disp('RMSE xbFilter dxbFilter ddxbFilter')
disp([rmse(xbFilter, xb) rmse(dxbFilter, dxb) rmse(ddxbFilter, ddxb)])
[XEst,RbsEst,taubEst,nubEst,output]=dynamicSfm(xbFilter,dxbFilter,ddxbFilter,Gammab,wIMUNoise,alphaIMUNoise,...
    'inertiamatrix',J,'RFirstFrame',Rbs(:,:,1),'showDiagnosticInfo',...
    'projectmfactorrotations',...
    'rotationExtractionWeights',[3 1 1],...
    'rotationDerivativeConstraint',TSampling,...
    'rotationDerivativeConstraintWeight',10,...
    ...%'debugSubstituteRotations',Rbs,...
    'translationDerivativeConstraint','sgolay',1,3,TSampling,...
    'translationDerivativeConstraintWeight',1,'flagTranslationDerivativeConstraint',true,...
    'translationSecondDerivativeConstraint','sgolay',1,3,TSampling,...
    'translationSecondDerivativeConstraintWeight',1,'flagTranslationSecondDerivativeConstraint',true,...
    'flagEstimateAlignedGravity',true,...
    'retriangulateStructure'...
    );

disp('# Dynamic Affine SfM')
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
[ddTsbEst,idxValid]=convVectorized(dTsbEst,output.secondFilterCoefficients,'same');

[~,~,TrProcrustes]=procrustes([X Tsb]',[XEst TsbEst]','reflection',false,'scaling',false);
RProcrustes=TrProcrustes.b*TrProcrustes.T';
TProcrustes=TrProcrustes.c(1,:)';
TsbEst=rigidTransform(RProcrustes,TProcrustes,TsbEst);
dTsbEst=rigidTransform(RProcrustes,zeros(3,1),dTsbEst);
ddTsbEst=rigidTransform(RProcrustes,zeros(3,1),ddTsbEst);
XEst=rigidTransform(RProcrustes,TProcrustes,XEst);

[~,TsbIntegrated]=dynamicSfM_integrateIMU(t,wIMUNoise,alphaIMUNoise,...
    Rbs(:,:,1),Tsb(:,1),dTsbEst(:,1));

[XEstAffine,RbsEstAffine,taubEstAffine]=affineSfM(xbFilter,'flagFlip',true,'RFirstFrame',Rbs(:,:,1));
TsbEstAffine=multiprodMatVec(invR(RbsEstAffine),taubEstAffine);
[~,~,TrProcrustesAffine]=procrustes(X',XEstAffine','reflection',false,'scaling',false);
RProcrustesAffine=TrProcrustesAffine.b*TrProcrustesAffine.T';
TProcrustesAffine=TrProcrustesAffine.c(1,:)';
TsbEstAffine=rigidTransform(RProcrustesAffine,TProcrustesAffine,TsbEstAffine);
XEstAffine=rigidTransform(RProcrustesAffine,TProcrustesAffine,XEstAffine);

disp('# Affine SfM (images only)')
disp('RMSE X')
disp(rmse(XEstAffine,X))
disp('RMS Distance Rbs [deg]')
disp(rmse(rot_dist(RbsEstAffine,Rbs,'vector')*180/pi))
disp('RMSE translation')
disp(rmse(taubEstAffine,taub))


figure(1)
plotPoints(X,'ro')
hold on
plotPoints(XEst,'bx')
plotPoints(XEstAffine,'bd')
plotLines(X,XEst,'k-')
plotPoints(Tsb,'r:')
plotPoints(TsbEst,'b-')
plotPoints(TsbIntegrated,'k-')
plotPoints(TsbEstAffine,'m--')

draw3dcameraFromRT(invR(RbsEst(:,:,1)),TsbEst(:,1),'scale',0.2)
hold off
view(0,90)
axis equal

figure(2)
subplot(2,2,1)
plot(t,rot_logVec(Rbs(:,:,1),Rbs)',':',t,rot_logVec(Rbs(:,:,1),RbsEst)',...
    t,rot_logVec(Rbs(:,:,1),RbsEstAffine)','--')
title('Rbs')
subplot(2,2,2)
plot(t,Tsb',':',t,TsbEst',t,TsbEstAffine','--')
title('Tsb')
subplot(2,2,3)
plot(t,dTsb',':',t,dTsbEst')
title('dTsb')
subplot(2,2,4)
plot(t,ddTsb',':',t(idxValid),ddTsbEst(:,idxValid)')
title('ddTsb')

fileNameSave=[fileName '_testNoise_data'];
fprintf('Saving results to %s\n',fileNameSave)
save(fileNameSave)

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
