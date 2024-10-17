function POCCameraDepthSensorCalibration2
resetRands()
flagDisplay=false;
flagDisplayAvg=true;
flagCheckDepthPlaneEstimates=false;
flagNoisy=true;

if ~flagNoisy
    sigmaNoisePoses=0;
    sigmaNoiseNormals=0;
    sigmaNoiseDistances=0;
else
    sigmaNoisePoses=0.05;
    sigmaNoiseNormals=0.05;
    sigmaNoiseDistances=0.05;
end

%use three coplanar known targets and pose estimation
methodAbsolutePoses='references';

L=3;
offset=[L/2;-L/2;0];
GTag1=RT2G(eye(3),zeros(3,1)+offset);
GTag2=RT2G(eye(3),[-L;0;0]+offset);
GTag3=RT2G(eye(3),[0;L;0]+offset);

GTag=cat(3,GTag1,GTag2,GTag3);

NTags=size(GTag,3);
NPlanes=4;
NVec=[0;0;1;0];

%pose camera
GCamera=eye(4);
%pose depth sensor
GDepth=RT2G(rot_randn(eye(3),0.1),[0;3;0]);

%pose planes and tags on planes
GPlanes=zeros(4,4,NPlanes);
NPlaneTags=NPlanes*NTags;
GPlaneTags=zeros(4,4,NPlaneTags);
idxPlaneTags=reshape(1:NPlaneTags,NTags,NPlanes);
for iPlane=1:NPlanes
    R=rot_randn(eye(3),1);
    T=R*[0;0;5+2*rand];
    GPlanes(:,:,iPlane)=RT2G(R,T);
    for iTag=1:NTags
        GPlaneTags(:,:,idxPlaneTags(iTag,iPlane))=GPlanes(:,:,iPlane)*GTag(:,:,iTag);
    end
end
GPlaneTags=noiserigid(GPlaneTags,sigmaNoisePoses);

%plane equations in world's reference frame
NVecPlanes=squeeze(rigidTransformG(invg(GPlanes),NVec,'planes',methodAbsolutePoses));

%data in depth sensor's reference frame
NVecPlanesDepth=rigidTransformG(GDepth,NVecPlanes,'planes',methodAbsolutePoses,'wc');
NVecPlanesDepth=planeRandn(NVecPlanesDepth,sigmaNoiseNormals,sigmaNoiseDistances);

%relative depth sensor to plane poses
GijPlanesDepthEst=planeToG(NVecPlanesDepth,methodAbsolutePoses);
GPlanesDepth=computeRelativePoseFromG(GDepth,GPlanes,methodAbsolutePoses);
GammaPlanesDepth=zeros(6,6,NPlanes);
GammaPlanesDepthInv=zeros(6,6,NPlanes);
for iPlane=1:NPlanes
   GammaPlanesDepth(:,:,iPlane)=computeGammaPlane(GijPlanesDepthEst(:,:,iPlane));
   GammaPlanesDepthInv(:,:,iPlane)=computeGammaPlaneInverse(GijPlanesDepthEst(:,:,iPlane));
end

if flagCheckDepthPlaneEstimates
    for iPlane=1:NPlanes
       disp(computeError(GPlanesDepth(:,:,iPlane),GijPlanesDepthEst(:,:,iPlane),GammaPlanesDepth(:,:,iPlane))) 
       disp(computeError(invg(GPlanesDepth(:,:,iPlane)),invg(GijPlanesDepthEst(:,:,iPlane)),GammaPlanesDepthInv(:,:,iPlane))) 
    end
end

%relative poses for tags in the same plane
e3=[0;0;1];
GammaPlaneTags=repmat(blkdiag(eye(3),e3*e3'),[1,1,NPlaneTags]);
% GammaPlaneTags=zeros(6,6,NPlaneTags);
% for iPlaneTag=1:NPlaneTags
%     GammaPlaneTags(:,:,iPlaneTag)=computeGammaPlane(GPlaneTags(:,:,iPlaneTag));
% end
%disp(computeError(GPlaneTags(:,:,3),GPlaneTags(:,:,2),computeGammaPlane(GPlaneTags(:,:,1))));

%% Prepare data to build network
%Note: the nodes are organized in the following order
%   Plane 1, Tag 1
%   Plane 1, Tag 2
%   ...
%   Plane 2, Tag 1
%   ...
%   Camera
%   Depth sensor

%Intra-plane constraints
EOnePlaneTags=[reshape((1:NTags)'*ones(1,NTags),[],1) reshape(ones(NTags,1)*(1:NTags),[],1)];
EOnePlaneTags=EOnePlaneTags(EOnePlaneTags(:,1)~=EOnePlaneTags(:,2),:);

NEOnePlaneTags=size(EOnePlaneTags,1);
NEPlaneTags=NPlanes*NEOnePlaneTags;
EPlaneTags=zeros(NEPlaneTags,2);
idxPlanePlaneTags=1:NEOnePlaneTags:NEPlaneTags;
for iPlane=1:NPlanes
    EPlaneTags(idxPlanePlaneTags(iPlane):(idxPlanePlaneTags(iPlane)+NEOnePlaneTags-1),:)=...
        EOnePlaneTags;
    EOnePlaneTags=EOnePlaneTags+NTags;
end

GijPlaneTags=repmat(eye(4),[1,1,2*NPlaneTags]);%GPlaneTags(:,:,EPlaneTags(:,1));
GammaijPlaneTags=GammaPlaneTags(:,:,EPlaneTags(:,1));

%Camera-plane constraints
EPlanesCamera=[NPlaneTags+ones(NPlaneTags,1) (1:NPlaneTags)'];
EPlanesCamera=[EPlanesCamera;fliplr(EPlanesCamera)]; %symmetrize
GijPlanesCamera=cat(3,GPlaneTags,invg(GPlaneTags));
GammaijPlanesCamera=repmat(eye(6),[1,1,2*NPlaneTags]);

%Depth-plane constraints
EPlanesDepth=[NPlaneTags+1+ones(NPlaneTags,1) (1:NPlaneTags)'];
EPlanesDepth=[EPlanesDepth;fliplr(EPlanesDepth)]; %symmetrize
idxRepeat=reshape(repmat(1:NPlanes,NTags,1),[],1);
GijPlanesDepth=cat(3,GPlanesDepth(:,:,idxRepeat),invg(GPlanesDepth(:,:,idxRepeat)));
GammaijPlanesDepth=cat(3,GammaPlanesDepth(:,:,idxRepeat),GammaPlanesDepthInv(:,:,idxRepeat));

NVecPlaneTagsDepth=NVecPlanesDepth(:,reshape(repmat(1:NPlanes,NTags,1),[],1));
[nPlaneTagDepth,dPlaneTagDepth]=planeNVecToNd(NVecPlaneTagsDepth);
GDepthEst=poseEstimationFromNormalsPosesG(nPlaneTagDepth,dPlaneTagDepth,GPlaneTags);



%Combine ground truth, edges, relative poses and dispersion matrices
Gi=cat(3,GPlaneTags,GCamera,GDepth);
GiInit=cat(3,GPlaneTags,eye(4),GDepthEst);
E=[EPlaneTags; EPlanesCamera; EPlanesDepth];
Gij=cat(3,GijPlaneTags,GijPlanesCamera,GijPlanesDepth);
Gammaij=cat(3,GammaijPlaneTags,GammaijPlanesCamera,GammaijPlanesDepth);

%Build network structure
t_node=testNetworkCreateStruct(E,'edges',NPlaneTags+2);
t_node=testNetworkAddGroundTruth(t_node,Gi);
t_node=testNetworkAddMeasurements(t_node,'method','given',Gij);
t_node=testNetworkAddDispersionMatricesRT(t_node,'given',Gammaij);
t_node=testNetworkInitializeStates(t_node,'G',noiserigid(GiInit,0.01));

[Ri,Ti]=G2RT(GiInit);
[Rij,Tij]=G2RT(Gij);

disp(logLikelihoodNetwork(Ri,Ti,Rij,Tij,Gammaij,E))
optsLocalization={'noinit','displayIt','showMessages',...
    'optsLieMinimize',{'gradientOnly'}};
t_node=localization_MLE_rigid(t_node,optsLocalization{:});

save([mfilename '_data'])

GiAvg=t_node.gi;
GPlaneTagsAvg=GiAvg(:,:,1:NPlaneTags);
GCameraAvg=GiAvg(:,:,NPlaneTags+1);
GDepthAvg=GiAvg(:,:,NPlaneTags+2);

GComp=invg(GCameraAvg);
for iPlaneTags=1:NPlaneTags
    GPlaneTagsAvg(:,:,iPlaneTags)=GComp*GPlaneTagsAvg(:,:,iPlaneTags);
end
GCameraAvg=GComp*GCameraAvg;
GDepthAvg=GComp*GDepthAvg;

if flagDisplayAvg
    draw3dcameraFromG(GPlaneTagsAvg,methodAbsolutePoses,'shape','april')
    hold on
    draw3dcameraFromG(GCameraAvg,methodAbsolutePoses,'color','g')
    draw3dcameraFromG(GDepthAvg,methodAbsolutePoses,'color','b')
    draw3dPlane(NVecPlanes,'side',3)
    axis equal
end

if flagDisplay
    draw3dcameraFromG(GPlaneTags,methodAbsolutePoses,'shape','april')
    hold on
    draw3dcameraFromG(GCamera,methodAbsolutePoses,'color','g')
    draw3dcameraFromG(GDepth,methodAbsolutePoses,'color','b')
    draw3dPlane(NVecPlanes,'side',3)
    axis equal
end

function [e,v]=computeError(GTruth,G,Gamma)
v=rot3r3_vee(GTruth,rot3r3_log(GTruth,G));
e=(v'*Gamma*v)/2;

function Gamma=computeGammaEpipole(G)
T=G2T(G);
Gamma=blkdiag(eye(3),orthComplementProjector(T));

function Gamma=computeGammaPlane(G)
R=G2R(G);
e3=[0;0;1];
Gamma=blkdiag(orthComplementProjector(e3),R(:,3)*R(:,3)');

function Gamma=computeGammaPlaneInverse(G)
R=G2R(G);
e3=[0;0;1];
Gamma=blkdiag(orthComplementProjector(R(:,3)),e3*e3');
