function POCNewMultiCalibration
c=multiCalbration_datasetToCorrespondencesFile('calibrationDataset');
c=multiCalibration_correspondencesProcess(c);
flagComputeRelPoses=false;
flagUseCovariances=false;
NTrials=1000;
optsRansac={'Ntrials',NTrials,'collectResiduals','inlierEstimate','waitBar'};

if flagComputeRelPoses
    %%% Relative pose estimation 
    disp('# Relative pose estimation')

    NCorrespondeces=length(c);
    poses=[];
    output=cell(NCorrespondeces,1);

    for iCorrespondence=1:NCorrespondeces
        cCurrent=c{iCorrespondence};
        pCurrent=[];
        switch cCurrent.type
            case 'plane-plane'
                ransacThreshold=0.1;
                [pCurrent.G,outputCurrent]=poseEstimationFromNormalsDistancesGRansac(...
                    cCurrent.n1,cCurrent.d1,...
                    cCurrent.n2,cCurrent.d2,...
                    ransacThreshold,'optsRansac',optsRansac);
                fInliers=outputCurrent.flagInliers;
                pCurrent.Sigma=poseEstimationFromNormalsDistancesCovarianceG(pCurrent.G,...
                    cCurrent.n1(:,fInliers),cCurrent.d1(fInliers),...
                    cCurrent.n2(:,fInliers),cCurrent.d2(fInliers),...
                    0.02,0.001);
            case 'plane-point'
                ransacThreshold=0.05;
                [pCurrent.G,outputCurrent]=poseEstimationFromNormalsPointsGRansac(...
                    cCurrent.n1,cCurrent.d1,cCurrent.x22D,ransacThreshold,...
                    'optsRansac',optsRansac,'refineEnd');

                fInliers=outputCurrent.flagInliers;
                pCurrent.Sigma=poseEstimationFromNormalsPointsCovarianceG(pCurrent.G,...
                    cCurrent.n1(:,fInliers),cCurrent.d1(fInliers),cCurrent.x22D(:,fInliers),...
                    0.02,0.02,0.02);
        end
        poses(iCorrespondence).G=pCurrent.G;
        poses(iCorrespondence).Sigma=pCurrent.Sigma;
        output{iCorrespondence}=outputCurrent;
    end
else
    load([mfilename '_data'])
end

disp('# Absolute pose estimation through averaging')
[t_node,nodeNames]=multiCalibration_networkPrepare(poses,c,flagUseCovariances);
optsLocalization={'showMessages'};
t_node=localization_MLE_rigid(t_node,optsLocalization{:});

testNetworkDisplayErrors(t_node,'rtn')

save([mfilename '_data'])
