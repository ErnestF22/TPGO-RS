function sfm_rawRotationsCombineFactorization_test
[E,Rij,RiRef,RijNoise]=sfm_utilityRotationDataset();

allOpts={
    {'methodLowRank','none'}
%     {'sdp'}
%     {'sdp','optsLowRank',{'spectrahedral'}}
%     {'sdp','projectBeforeFactorization'}
%     {'sdp','projectBeforeFactorization','optsLowRank',{'spectrahedral'}}
    {'alm','optsLowRank',{'criterion','l1','optsOptimizer',{'lowrankprior','SDP'}}}
    {'alm','optsLowRank',{'criterion','l1','innerProj',false,'optsOptimizer',{'lowrankprior','SDP'}}}
    {'alm','optsLowRank',{'criterion','l12','optsOptimizer',{'lowrankprior','SDP'}}}
    {'alm','optsLowRank',{'criterion','l12','innerProj',false,'optsOptimizer',{'lowrankprior','SDP'}}}
%     {'alm'}
%     {'alm','optsLowRank',{'criterion','l2'}}
%     {'alm','optsLowRank',{'criterion','l2','innerProj',false}}
%     {'alm','optsLowRank',{'criterion','l1'}}
%     {'alm','optsLowRank',{'criterion','l1','innerProj',false}}
%     {'alm','optsLowRank',{'criterion','l12'}}
%     {'alm','optsLowRank',{'criterion','l12','innerProj',false}}
%     {'normalize'}
%     {'laplacian'}
%     {'laplacian','normalize'}
    };

disp('Errors with clean measurements (absolute poses) [deg] (mean/median) (should be around zero)')
testMethods(Rij)
disp('Errors with noisy measurements (absolute poses) [deg] (mean/median) (should be less than 6.0')
testMethods(RijNoise)

    function testMethods(Rij)
    NOpts=length(allOpts);
    for iOpts=1:NOpts
        opts=allOpts{iOpts};
        Ri=sfm_rawRotationsCombineFactorization(Rij,E,opts{:});
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
