%Prepare network for averaging poses
%function [t_node,idxCamera,idxVelodyne,idxHokuyo]=multiCalibrationPrepareNetwork(GCameraVelodyneEst,GHokuyoCameraEst,GHokuyoVelodyneEst)
function [t_node,idxCamera,idxVelodyne,idxHokuyo]=multiCalibrationPrepareNetwork(poses,flagUseCovariances)
idxCamera=1;
idxVelodyne=2;
idxHokuyo=3;
NNodes=3;

E=[idxCamera idxVelodyne; idxCamera idxHokuyo; idxVelodyne idxHokuyo];
E=[E;fliplr(E)];

t_node=testNetworkCreateStruct(E,'Edges',NNodes,'Type','Single');
E=t_node.E;

NEdges=size(E,1);
gij=zeros(4,4,NEdges);
Gammaij=zeros(6,6,NEdges);
for iEdge=1:NEdges
    switch toStr(E(iEdge,:))
        case toStr([idxCamera idxVelodyne])
            g=invg(poses.GCameraVelodyne);
            Sigma=poses.SigmaCameraVelodyne;
            
        case toStr([idxVelodyne idxCamera])
            g=poses.GCameraVelodyne;
            Sigma=poses.SigmaCameraVelodyne;
            
        case toStr([idxCamera idxHokuyo])
            g=poses.GHokuyoCamera;
            Sigma=poses.SigmaHokuyoCamera;
            
        case toStr([idxHokuyo idxCamera])
            g=invg(poses.GHokuyoCamera);
            Sigma=poses.SigmaHokuyoCamera;
            
        case toStr([idxVelodyne idxHokuyo])
            g=poses.GHokuyoVelodyne;
            Sigma=poses.SigmaHokuyoVelodyne;
            
        case toStr([idxHokuyo idxVelodyne])
            g=invg(poses.GHokuyoVelodyne);
            Sigma=poses.SigmaHokuyoVelodyne;
            
        otherwise
            error('Edge node combination not recognized');
    end
    gij(:,:,iEdge)=g;
    Gammaij(:,:,iEdge)=pinv(Sigma);
end

t_node=testNetworkAddMeasurements(t_node,'Method','Given',gij);
if ~flagUseCovariances
    t_node=testNetworkAddDispersionMatricesRT(t_node,'identity');
else
    t_node=testNetworkAddDispersionMatricesRT(t_node,'given',Gammaij);
end

giInit=noiserigid(cat(3,eye(4),poses.GCameraVelodyne,invg(poses.GHokuyoCamera)),0.2);
t_node=testNetworkInitializeStates(t_node,'G',invg(giInit));

function s=toStr(a)
s=num2str(a,'%02d');
