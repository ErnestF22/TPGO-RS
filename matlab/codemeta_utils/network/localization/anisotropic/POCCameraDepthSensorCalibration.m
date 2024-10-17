function POCCameraDepthSensorCalibration
resetRands(2)
flagCheckPlanes=false;
flagCheckHomographyEstimates=true;
flagCheckDepthPlaneEstimates=false;
flagDisplay=false;
flagDisplayDepth=false;                             
flagDisplayAveragedSolution=true;
flagNoisy=false;

if ~flagNoisy
    sigmaNoisePixel=0;
    sigmaNoiseNormals=0;
    sigmaNoiseDistances=0;
else
    sigmaNoisePixel=0.01;
    sigmaNoiseNormals=0.05;
    sigmaNoiseDistances=0.05;
end

methodAbsolutePoses='references';

NPoints=20;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
NPlanes=5;
NVec=[0;0;1;0];
X=planeGeneratePoints(NVec,NPoints);

GCamera=eye(4);
GDepth=RT2G(rot_randn(eye(3),0.1),[0;3;0]);
GPlanes=zeros(4,4,NPlanes);
for iPlane=1:NPlanes
    R=rot_randn(eye(3),1);
    T=R*[0;0;5+2*rand];
    GPlanes(:,:,iPlane)=RT2G(R,T);
end
%prepare dataset (in world's reference frame)
XPlanes=rigidTransformG(invg(GPlanes),X,methodAbsolutePoses);
NVecPlanes=squeeze(rigidTransformG(invg(GPlanes),NVec,'planes',methodAbsolutePoses));

%data in depth sensor's reference frame
NVecPlanesDepth=rigidTransformG(GDepth,NVecPlanes,'planes',methodAbsolutePoses,'wc');
NVecPlanesDepth=planeRandn(NVecPlanesDepth,sigmaNoiseNormals,sigmaNoiseDistances);

if flagDisplayDepth || flagDisplayAveragedSolution
    XPlanesDepth=rigidTransformG(GDepth,XPlanes,methodAbsolutePoses,'wc');
end

%check plane residuals
if flagCheckPlanes
    for iPlane=1:NPlanes
        disp(max(abs(planeResiduals(NVecPlanes(:,iPlane),XPlanes(:,:,iPlane)))))
    end
end


%relative plane poses from homography
xPlanes=projectFromG(GCamera,XPlanes,methodAbsolutePoses);
xPlanes=xPlanes+sigmaNoisePixel*randn(size(xPlanes));

EHomography=[reshape((1:NPlanes)'*ones(1,NPlanes),[],1) reshape(ones(NPlanes,1)*(1:NPlanes),[],1)];
EHomography=EHomography(EHomography(:,1)~=EHomography(:,2),:);

NEdgesHomography=size(EHomography,1);
GEdgesHomography=zeros(4,4,NEdgesHomography);
GammaEdgesHomography=zeros(6,6,NEdgesHomography);
for iEdge=1:NEdgesHomography
    iNode=EHomography(iEdge,1);
    jNode=EHomography(iEdge,2);
    xINode=xPlanes(:,:,iNode);
    xJNode=xPlanes(:,:,jNode);

    [GEdgesHomography(:,:,iEdge),NiEst]=poseEstimationHomographyG(xJNode,xINode);
    GEdgesHomography(1:3,4,iEdge)=GEdgesHomography(1:3,4,iEdge)/sign(NiEst(3));
    GammaEdgesHomography(:,:,iEdge)=computeGammaEpipole(GEdgesHomography(:,:,iEdge));
    if flagCheckHomographyEstimates
        Gij=GPlanes(:,:,iNode)*invg(GPlanes(:,:,jNode));
        disp(computeError(Gij,GEdgesHomography(:,:,iEdge),GammaEdgesHomography(:,:,iEdge)));
    end
end

%relative depth sensor to plane poses
GPlanesDepthEst=planeToG(NVecPlanesDepth,methodAbsolutePoses);
GPlanesDepth=computeRelativePoseFromG(GDepth,GPlanes,methodAbsolutePoses);
GammaPlanesDepth=zeros(6,6,NPlanes);
GammaPlanesDepthInv=zeros(6,6,NPlanes);
for iPlane=1:NPlanes
   GammaPlanesDepth(:,:,iPlane)=computeGammaPlane(GPlanesDepthEst(:,:,iPlane));
   GammaPlanesDepthInv(:,:,iPlane)=computeGammaPlaneInverse(GPlanesDepthEst(:,:,iPlane));
end

if flagCheckDepthPlaneEstimates
    for iPlane=1:NPlanes
       disp(computeError(GPlanesDepth(:,:,iPlane),GPlanesDepthEst(:,:,iPlane),GammaPlanesDepth(:,:,iPlane))) 
       disp(computeError(invg(GPlanesDepth(:,:,iPlane)),invg(GPlanesDepthEst(:,:,iPlane)),GammaPlanesDepthInv(:,:,iPlane))) 
    end
end

EPlanesDepth=[NPlanes+1*ones(NPlanes,1) (1:NPlanes)'];
%symmetrize
EPlanesDepth=[EPlanesDepth;fliplr(EPlanesDepth)];
GPlanesDepthEst=cat(3,GPlanesDepthEst,invg(GPlanesDepthEst));
GPlanesDepth=cat(3,GPlanesDepth,invg(GPlanesDepth));
GammaPlanesDepth=cat(3,GammaPlanesDepth,GammaPlanesDepthInv);

E=[EHomography; EPlanesDepth];
EType=[2*ones(size(EHomography,1),1); ones(size(EPlanesDepth,1),1)];
Gamma=cat(3,GammaEdgesHomography,GammaPlanesDepth);
Gij=cat(3,GEdgesHomography,GPlanesDepthEst);
Gi=cat(3,GPlanes,GDepth);
GiInit=cat(3,GPlanesDepth,eye(4));
t_node=testNetworkCreateStruct(E,'edges',NPlanes+1);
t_node=testNetworkSetEdgeType(t_node,EType);
t_node=testNetworkAddGroundTruth(t_node,Gi);
t_node=testNetworkAddMeasurements(t_node,'method','given',Gij);
t_node=testNetworkAddDispersionMatricesRT(t_node,'given',Gamma);
t_node=testNetworkInitializeStates(t_node,'G',GiInit);

optsLocalization={'noinit','displayIt','showMessages',...
    'optsLieMinimize',{'gradientOnly'}};
t_node=localization_MLE_rigid(t_node,optsLocalization{:});

if flagDisplay
    %draw scene
    figure(1)
    draw3dcameraFromG(GCamera)
    hold on
    draw3dcameraFromG(GDepth,'color','g',methodAbsolutePoses)
    draw3dcameraFromG(GPlanes,'color','b',methodAbsolutePoses)
    plotPoints(XPlanes)
    if flagDisplayDepth
        plotPoints(XPlanesDepth,'*')
        draw3dcameraFromG(GPlanesDepthEst(:,:,1:NPlanes),'color','g',methodAbsolutePoses)
    end
    hold off
    axis equal
end

if flagDisplayAveragedSolution
    %display solution after averaging
    %first align poses so that the dept sensor is at the origin
    GDepthAvg=t_node.gi(:,:,NPlanes+1);
    GPlanesDepthAvg=zeros(size(t_node.gi));
    for iNode=1:NPlanes+1;
        GPlanesDepthAvg(:,:,iNode)=invg(GDepthAvg)*t_node.gi(:,:,iNode);
    end
    figure(2)
    draw3dcameraFromG(GPlanesDepthAvg)
    hold on
    plotPoints(XPlanesDepth,'*')
    hold off
    axis equal
end

save([mfilename '_data'])

function e=computeErrorEpipole(GTruth,G)
Gamma=computeGammaEpipole(G);
e=computeError(GTruth,G,Gamma);

function e=computeErrorPlane(GTruth,G)
Gamma=computeGammaPlane(G);
e=computeError(GTruth,G,Gamma);


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
