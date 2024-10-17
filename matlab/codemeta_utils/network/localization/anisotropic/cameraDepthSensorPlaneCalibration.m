function [GCameraAvg,GDepthAvg,GPlaneTagsAvg,GDepthLinear]=cameraDepthSensorPlaneCalibration(GPlaneTags,NVecPlanesDepth,varargin)
optsLocalization={'noinit',...
    'optsLieMinimize',{'gradientOnly'}};
flagCompensateReference=true;
coeffPlaneConstraints=1;

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'optslocalization'
            ivarargin=ivarargin+1;
            optsLocalization=[optsLocalization varargin{ivarargin}];
        case 'flagcompensatereference'
            ivarargin=ivarargin+1;
            flagCompensateReference=varargin{ivarargin};
        case 'coeffplaneconstraints'
            ivarargin=ivarargin+1;
            coeffPlaneConstraints=varargin{ivarargin};
        otherwise
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

methodAbsolutePoses='references';
NPlanes=size(NVecPlanesDepth,2);
NPlaneTags=size(GPlaneTags,3);
NTags=NPlaneTags/NPlanes;
if NTags~=round(NTags)
    error('size(GPlaneTags,3) must be divisible by size(NVecPlanesDepth,2)')
end

%% Transform input in pose data
%relative depth sensor to plane poses
GPlanesDepth=planeToG(NVecPlanesDepth,methodAbsolutePoses);

%prepare singular dispersion matrices
GammaPlanesDepth=zeros(6,6,NPlanes);
GammaPlanesDepthInv=zeros(6,6,NPlanes);
for iPlane=1:NPlanes
   GammaPlanesDepth(:,:,iPlane)=computeGammaPlane(GPlanesDepth(:,:,iPlane));
   GammaPlanesDepthInv(:,:,iPlane)=computeGammaPlaneInverse(GPlanesDepth(:,:,iPlane));
end

e3=[0;0;1];
GammaPlaneTags=repmat(blkdiag(eye(3),e3*e3'),[1,1,NPlaneTags]);

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
GammaijPlaneTags=coeffPlaneConstraints*GammaPlaneTags(:,:,EPlaneTags(:,1));

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
GDepthLinear=poseEstimationFromNormalsPosesG(nPlaneTagDepth,dPlaneTagDepth,GPlaneTags);

%Combine ground truth, edges, relative poses and dispersion matrices
GiInit=cat(3,GPlaneTags,eye(4),GDepthLinear);
E=[EPlaneTags; EPlanesCamera; EPlanesDepth];
Gij=cat(3,GijPlaneTags,GijPlanesCamera,GijPlanesDepth);
Gammaij=cat(3,GammaijPlaneTags,GammaijPlanesCamera,GammaijPlanesDepth);

%% Build network
t_node=testNetworkCreateStruct(E,'edges',NPlaneTags+2);
t_node=testNetworkAddMeasurements(t_node,'method','given',Gij);
t_node=testNetworkAddDispersionMatricesRT(t_node,'given',Gammaij);
t_node=testNetworkInitializeStates(t_node,'G',GiInit);

%% Optimize
t_node=localization_MLE_rigid(t_node,optsLocalization{:});

GiAvg=t_node.gi;
GPlaneTagsAvg=GiAvg(:,:,1:NPlaneTags);
GCameraAvg=GiAvg(:,:,NPlaneTags+1);
GDepthAvg=GiAvg(:,:,NPlaneTags+2);

%% Put world's reference at the camera if requested
if flagCompensateReference
    GComp=invg(GCameraAvg);
    for iPlaneTags=1:NPlaneTags
        GPlaneTagsAvg(:,:,iPlaneTags)=GComp*GPlaneTagsAvg(:,:,iPlaneTags);
    end
    GCameraAvg=GComp*GCameraAvg;
    GDepthAvg=GComp*GDepthAvg;
end

function Gamma=computeGammaPlane(G)
R=G2R(G);
e3=[0;0;1];
Gamma=blkdiag(orthComplementProjector(e3),R(:,3)*R(:,3)');

function Gamma=computeGammaPlaneInverse(G)
R=G2R(G);
e3=[0;0;1];
Gamma=blkdiag(orthComplementProjector(R(:,3)),e3*e3');

