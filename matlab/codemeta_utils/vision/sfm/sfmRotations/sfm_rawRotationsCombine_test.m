function sfm_rawRotationsCombine_test
resetRands()
load([mfilename 'Dataset'],'E','Rij','RiRef','RijNoise')

allOpts={
    {'optsGlobal',{'methodLowRank','none'}}
    {'optsGlobal',{'sdp'}}
    {'optsGlobal',{'alm'}}
    {'optsGlobal',{'alm','optsLowRank',{'criterion','l2'}}}
    {'optsGlobal',{'alm','optsLowRank',{'criterion','l2','innerProj',false}}}
    {'optsGlobal',{'alm','optsLowRank',{'criterion','l1'}}}
    {'optsGlobal',{'alm','optsLowRank',{'criterion','l1','innerProj',false}}}
    {'optsGlobal',{'alm','optsLowRank',{'criterion','l12'}}}
    {'optsGlobal',{'alm','optsLowRank',{'criterion','l12','innerProj',false}}}
    {'optsGlobal',{'normalize'}}
    {'optsGlobal',{'laplacian'}}
    {'optsGlobal',{'laplacian','normalize'}}
    {'inferenceOutliers','optsGlobal',{'methodLowRank','none'}}
    {'localRefine','optsGlobal',{'methodLowRank','none'}}
    {'localRefine','optsLocal',{'method','Weiszfeld'},'optsGlobal',{'methodLowRank','none'}}
    {'inferenceOutliers','localRefine','optsGlobal',{'methodLowRank','none'}}
    {'optsGlobal',{'sdp','projectBeforeFactorization'}}
    {'optsGlobal',{'alm','optsLowRank',{'criterion','l1','optsOptimizer',{'lowrankprior','SDP'}}}}
    {'optsGlobal',{'alm','optsLowRank',{'criterion','l1','innerProj',false,'optsOptimizer',{'lowrankprior','SDP'}}}}
    {'optsGlobal',{'alm','optsLowRank',{'criterion','l12','optsOptimizer',{'lowrankprior','SDP'}}}}
    {'optsGlobal',{'alm','optsLowRank',{'criterion','l12','innerProj',false,'optsOptimizer',{'lowrankprior','SDP'}}}}
    };

disp('Errors with clean measurements (absolute poses) [deg] (mean/median) (should be around zero)')
testMethods(Rij)
disp('Errors with noisy measurements (absolute poses) [deg] (mean/median) (should be around 6.0)')
testMethods(RijNoise)

    function testMethods(Rij)
    NOpts=length(allOpts);
    for iOpts=1:NOpts
        opts=allOpts{iOpts};
        Ri=sfm_rawRotationsCombine(Rij,E,opts{:});
        disp('  ---')
        disp(cellExpand(opts))
        displayErrors(Ri)
    end
    end
        
    function displayErrors(Ri)
    RiTransformed=rotationProcrustesAlign(RiRef,Ri,'left');
    e=rot_dist(RiRef,RiTransformed,'vector')'*180/pi;
    fprintf('%.4e / %.4e\n', mean(e), median(e))
    end
end
