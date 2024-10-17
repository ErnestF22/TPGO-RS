function dynamicSfm_testSmoothness
%load('sampleDynamicSfMDataset_sampling')
load('sampleDynamicSfMDataset_interpolation')
%methodFilter='given';
methodFilter='auto';

filterOrder=4;
filterWindow=21;
TSampling=mean(diff(t));

switch methodFilter
    case 'given'
        [~,g] = sgolay(filterOrder,filterWindow);   
        g=-g(:,2)/TSampling;
        optsFilter={'translationDerivativeConstraint',g};
    case 'auto'
        optsFilter={'translationDerivativeConstraint','sgolay',filterOrder,filterWindow,TSampling};
end
           
[XEst,RbsEst,taubEst,nubEst]=dynamicSfm(xb,dxb,ddxb,Gammab,wIMU,alphaIMU,...
    'inertiamatrix',J,'RFirstFrame',Rbs(:,:,1),'showDiagnosticInfo',...
    'methodTranslationFactor','reduced',...
    'methodTranslationExtraction','singlesystem',...
    optsFilter{:},'flagTranslationDerivativeConstraint',true);

disp('RMSE X')
disp(rmse(XEst,X))
disp('RMS Distance Rbs')
disp(rmse(rot_dist(RbsEst,Rbs,'vector')))
disp('RMSE translation')
disp(rmse(taubEst,taub))
disp('RMSE linear velocity')
disp(rmse(nubEst,nub))
