function dynamicSfm_test
%load('sampleDynamicSfMDataset_original')
%load('sampleDynamicSfMDataset_sampling')
load('sampleDynamicSfMDataset_interpolation')
[XEst,RbsEst,taubEst,nubEst]=dynamicSfm(xb,dxb,ddxb,Gammab,wIMU,alphaIMU,...
    'inertiamatrix',J,'RFirstFrame',Rbs(:,:,1),'showDiagnosticInfo',...
    'rotationExtractionWeights',[1 1 1],...%[1 0.1 0.01],...
    'rotationDerivativeConstraint',min(diff(t)),...
    'methodTranslationFactor','reduced',...
    'methodTranslationExtraction','singlesystem',...
    'flagEstimateAlignedGravity',false,...
    'retriangulateStructure');

disp('RMSE X')
disp(rmse(XEst,X))
disp('RMS Distance Rbs')
disp(rmse(rot_dist(RbsEst,Rbs,'vector')))
disp('RMSE translation')
disp(rmse(taubEst,taub))
disp('RMSE linear velocity')
disp(rmse(nubEst,nub))

[XEstAffine,RbsEstAffine,taubEstAffine]=affineSfM(xb,'flagFlip',true,'RFirstFrame',Rbs(:,:,1));
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

keyboard