function [posesAbsolute,posesRelativeMeasured,posesAbsoluteTree]=multiCalibrationNew(correspondences,varargin)
flagUseCovariances=false;
NTrials=1000;
optsRansac={'Ntrials',NTrials,'collectResiduals','inlierEstimate','waitBar'};

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'flagusecovariances'
            ivarargin=ivarargin+1;
            flagUseCovariances=varargin{ivarargin};
        otherwise
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

%%% Add relative pose estimation to correspondences structure
disp('# Relative pose estimation')

NCorrespondeces=length(correspondences);
posesRelativeMeasured=repmat(struct('G',eye(4),'Sigma',eye(6)),NCorrespondeces,1);
output=cell(NCorrespondeces,1);

for iCorrespondence=1:NCorrespondeces
    cCurrent=correspondences{iCorrespondence};
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
    posesRelativeMeasured(iCorrespondence).G=pCurrent.G;
    posesRelativeMeasured(iCorrespondence).Sigma=pCurrent.Sigma;
    output{iCorrespondence}=outputCurrent;
end

disp('# Absolute pose estimation through averaging')
[t_node,nodeNames]=multiCalibration_networkPrepare(posesRelativeMeasured,correspondences,flagUseCovariances);
ETree=testNetworkGenerateRandomTree(t_node);
t_node=testNetworkLocalizeTree(t_node,ETree,'References');
t_node.giInit=t_node.gi;

optsLocalization={'noinit'};%'showMessages'};
t_node=localization_MLE_rigid(t_node,optsLocalization{:});

% copy results in output variables
NNodes=testNetworkGetNumberOfNodes(t_node);
posesAbsolute=repmat(struct('G',eye(4),'nodeName',[]),NNodes,1);
posesAbsoluteTree=posesAbsolute;
for iNode=1:NNodes
    posesAbsolute(iNode).G=t_node.gi(:,:,iNode);
    posesAbsolute(iNode).nodeName=nodeNames{iNode};
    posesAbsoluteTree(iNode).G=t_node.giInit(:,:,iNode);
    posesAbsoluteTree(iNode).nodeName=nodeNames{iNode};
end
