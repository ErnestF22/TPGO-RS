function multiCalibration_test
resetRands();

%% Initialization of the dataset
methodAbsolutePoses='reference';
RCamera=rot(pi/4*[0;0;1])*rot(pi/2*[1;0;0]);
TCamera=[1;0;0];
GCamera=RT2G(RCamera,TCamera);

RVelodyne=rot(3*pi/7*[0;0;1])*rot(pi/2*[1;0;0])*rot(pi/8*[0;0;1]);
TVelodyne=[0;-1;-1];
GVelodyne=RT2G(RVelodyne,TVelodyne);

RHokuyo=rot(2*pi/7*[0;0;1])*rot(pi/2*[1;0;0])*rot(7*pi/16*[0;0;1]);
THokuyo=[0;1;-1];
GHokuyo=RT2G(RHokuyo,THokuyo);

hokuyoVectorPlane=rigidTransformG(GHokuyo,[1;0;0;0],'planes','methodAbsolutePoses',methodAbsolutePoses,'cw');
VectorPlanes=[   1  , 1   , 1  , 0.8, 0.8, 0.5 ;
                -0.5, 0   ,-0.5,-0.4, 0.1, -0.2;
                 0.5, 0   ,-0.5, 0.1,-0.4, 0.2;
                -4  ,-4.5 ,-3.5,-5  ,-5  ,-5];
VectorPlanes=planeNormalize(VectorPlanes);

VectorLines=planeIntersect(VectorPlanes,hokuyoVectorPlane);

NPlanes=size(VectorPlanes,2);

%One point in each plane and two points in each line
XPlanes=zeros(3,NPlanes);
XLines=zeros(3,2,NPlanes);
for iPlane=1:NPlanes
    XPlanes(:,iPlane)=homogeneousProjectPoints(zeros(3,1),VectorPlanes(:,iPlane));
    XLines(:,:,iPlane)=homogeneousProjectPoints(randn(3,2),VectorLines(:,:,iPlane));
end

%% Visualization
draw3dcameraFromRT(RCamera,TCamera,'methodAbsolutePoses',methodAbsolutePoses);
hold on
draw3dcameraFromRT(RVelodyne,TVelodyne,'methodAbsolutePoses',methodAbsolutePoses,'color1','c');
draw3dcameraFromRT(RHokuyo,THokuyo,'methodAbsolutePoses',methodAbsolutePoses,'color1','m','scale',[0 1]);
for iPlane=1:NPlanes
    draw3dPlane(VectorPlanes(:,iPlane),'center',TCamera,'side',2)
    draw3dLine(VectorLines(:,:,iPlane),'style','r','side',8)
end
draw3dPlane(hokuyoVectorPlane,'center',THokuyo,'side',3,'style','r')
plotPoints(XPlanes,{'color','g','Marker','*'})
plotPoints(XLines(1:3,:),{'color','g','Marker','+'})
plot3(0,0,0,'bo')
xlabel('x')
ylabel('y')
zlabel('z')
axis equal
hold off

%% Data transformation in local coordinates and checks
disp('Check points belong to planes')
disp(diag(VectorPlanes'*homogeneous(XPlanes,4))')
disp(computePointPlaneErrors(XLines,VectorPlanes))

disp('Check points belong to lines')
disp(computePointLineErrors(XLines,VectorLines));

VectorPlanesCamera=rigidTransformG(GCamera,VectorPlanes,'planes','methodAbsolutePoses',methodAbsolutePoses,'wc');
XPlanesCamera=rigidTransformG(GCamera,XPlanes,'points','methodAbsolutePoses',methodAbsolutePoses,'wc');
disp('Check points belong to planes in camera coordinates')
disp(diag(VectorPlanesCamera'*homogeneous(XPlanesCamera,4))')

VectorPlanesVelodyne=rigidTransformG(GVelodyne,VectorPlanes,'planes','methodAbsolutePoses',methodAbsolutePoses,'wc');
XPlanesVelodyne=rigidTransformG(GVelodyne,XPlanes,'points','methodAbsolutePoses',methodAbsolutePoses,'wc');
disp('Check points belong to planes in velodyne coordinates')
disp(diag(VectorPlanesVelodyne'*homogeneous(XPlanesVelodyne,4))')

VectorLinesHokuyo=rigidTransformG(GHokuyo,VectorLines,'planes','methodAbsolutePoses',methodAbsolutePoses,'wc');
XLinesHokuyo=rigidTransformG(GHokuyo,XLines,'points','methodAbsolutePoses',methodAbsolutePoses,'wc');
disp('Check points belong to lines in Hokuyo coordinates')
disp(computePointLineErrors(XLinesHokuyo,VectorLinesHokuyo));
disp('Check zero x coordinate for points on lines in Hokuyo''s reference frame')
disp(XLinesHokuyo(1,:))

normalPlanesCamera=VectorPlanesCamera(1:3,:);
normalPlanesVelodyne=VectorPlanesVelodyne(1:3,:);
distancePlanesCamera=-VectorPlanesCamera(4,:);
distancePlanesVelodyne=-VectorPlanesVelodyne(4,:);

GCameraVelodyne=computeRelativePoseFromG(GVelodyne,GCamera,'methodAbsolutePoses',methodAbsolutePoses);
[RCameraVelodyne,TCameraVelodyne]=G2RT(GCameraVelodyne);
disp('Check transformation GCameraVelodyne maps points in the two frames correctly')
disp(norm(XPlanesVelodyne-[eye(3) zeros(3,1)]*GCameraVelodyne*homogeneous(XPlanesCamera,4),'fro'))
disp(norm(XPlanesVelodyne-rigidTransformG(GCameraVelodyne,XPlanesCamera),'fro'))

disp('Check tranformation RCameraVelodyne maps normals between the two frames correctly')
disp(norm(RCameraVelodyne*normalPlanesCamera-normalPlanesVelodyne,'fro'))

GHokuyoCamera=computeRelativePoseFromG(GCamera,GHokuyo,'methodAbsolutePoses',methodAbsolutePoses);
[RHokuyoCamera,THokuyoCamera]=G2RT(GHokuyoCamera);
XLinesCamera=rigidTransformG(GCamera,XLines,'points','methodAbsolutePoses',methodAbsolutePoses,'wc');
disp('Check transformation GCameraHokuyo maps points in the two frames correctly')
disp(norm(XLinesCamera(1:3,:)-rigidTransformG(GHokuyoCamera,XLinesHokuyo(1:3,:)),'fro'))

GHokuyoVelodyne=computeRelativePoseFromG(GVelodyne,GHokuyo,'methodAbsolutePoses',methodAbsolutePoses);
%[RHokuyoVelodyne,THokuyoVelodyne]=G2RT(GHokuyoVelodyne);
XLinesVelodyne=rigidTransformG(GVelodyne,XLines,'points','methodAbsolutePoses',methodAbsolutePoses,'wc');
disp('Check transformation GVelodyneHokuyo maps points in the two frames correctly')
disp(norm(XLinesVelodyne(1:3,:)-rigidTransformG(GHokuyoVelodyne,XLinesHokuyo(1:3,:)),'fro'))

%Compute line directions from points
lHokuyo=zeros(3,NPlanes);
for iPlane=1:NPlanes
    lHokuyo(:,iPlane)=XLinesHokuyo(:,:,iPlane)*[1;-1];
end

disp('Check relation between normals in camera coordinates and line directions in Hokuyo coordinates')
e=zeros(1,NPlanes);
for iPlane=1:NPlanes
    e(iPlane)=normalPlanesCamera(:,iPlane)'*RHokuyoCamera*lHokuyo(:,iPlane);
end
disp(e)

disp('Check relation between normals in camera coordinates and points in Hokuyo coordinates')
normalPlanesCameraFlat=reshape(repmat(normalPlanesCamera,[2 1]),3,[]);
distancePlanesCameraFlat=reshape(repmat(distancePlanesCamera,[2 1]),1,[]);
XLinesHokuyoFlat2D=XLinesHokuyo(2:3,:);
e=poseEstimationFromNormalsPointsResidualsG(GHokuyoCamera,normalPlanesCameraFlat,distancePlanesCameraFlat,XLinesHokuyoFlat2D);
disp(e')

disp('Check relation between normals in Velodyne coordinates and points in Hokuyo coordinates')
normalPlanesVelodyneFlat=reshape(repmat(normalPlanesVelodyne,[2 1]),3,[]);
distancePlanesVelodyneFlat=reshape(repmat(distancePlanesVelodyne,[2 1]),1,[]);
e=poseEstimationFromNormalsPointsResidualsG(GHokuyoVelodyne,normalPlanesVelodyneFlat,distancePlanesVelodyneFlat,XLinesHokuyoFlat2D);
disp(e')

%% Test Camera to Velodyne transformation estimation
GCameraVelodyneEst=poseEstimationFromNormalsDistancesG(normalPlanesVelodyne,distancePlanesVelodyne,normalPlanesCamera,distancePlanesCamera);
disp('[GCameraVelodyneEst GCameraVelodyne GCameraVelodyneEst-GCameraVelodyne]')
disp([GCameraVelodyneEst GCameraVelodyne GCameraVelodyneEst-GCameraVelodyne])

%compute and show residuals
[RCameraVelodyneEst,TCameraVelodyneEst]=G2RT(GCameraVelodyneEst);
[eRot,eTransl]=poseEstimationFromNormalsDistancesResiduals(RCameraVelodyneEst,TCameraVelodyneEst,normalPlanesVelodyne,distancePlanesVelodyne,normalPlanesCamera,distancePlanesCamera);
disp('Rotation and translation residuals from camera to Velodyne estimation')
disp(eRot')
disp(eTransl')

%test with outliers
normalPlanesVelodyneOutliers=normalPlanesVelodyne;
normalPlanesVelodyneOutliers(:,2)=sphere_randn(eye(3,1));
%distancePlanesCamera(3)=10*rand;
[GCameraVelodyneEstRansac,output]=poseEstimationFromNormalsDistancesGRansac(normalPlanesVelodyneOutliers,distancePlanesVelodyne,normalPlanesCamera,distancePlanesCamera,3,'optsRansac',{'collectResiduals'});
disp('[GCameraVelodyneEstRansac GCameraVelodyne GCameraVelodyneEstRansac-GCameraVelodyne]')
disp([GCameraVelodyneEstRansac GCameraVelodyne GCameraVelodyneEstRansac-GCameraVelodyne])

%% Test Hokuyo to Camera transformation estimation
GHokuyoCameraEst=poseEstimationFromNormalsPointsG(normalPlanesCameraFlat,distancePlanesCameraFlat,XLinesHokuyoFlat2D);
disp('[GHokuyoCameraEst GHokuyoCamera GHokuyoCameraEst-GHokuyoCamera]')
disp([GHokuyoCameraEst GHokuyoCamera GHokuyoCameraEst-GHokuyoCamera])
save('poseEstimationFromNormalsPointsRefineG_test_inputData','normalPlanesCameraFlat','distancePlanesCameraFlat','XLinesHokuyoFlat2D','GHokuyoCamera')

%compute and show errors
e=poseEstimationFromNormalsPointsResidualsG(GHokuyoCameraEst,normalPlanesCameraFlat,distancePlanesCameraFlat,XLinesHokuyoFlat2D);
disp('Residuals from camera to Hokuyo estimation')
disp(e)

% %GHokuyoCameraEst=poseEstimationFromNormalsDistancesPointsG(normalPlanesCamera,distancePlanesCamera,XLinesHokuyo);
% %compute and show residuals
% [RHokuyoCameraEst,THokuyoCameraEst]=G2RT(GHokuyoCameraEst);
% [eRot,eTransl]=poseEstimationFromNormalsDistancesPointsResiduals(RHokuyoCameraEst,THokuyoCameraEst,...
%     normalPlanesCamera,distancePlanesCamera,XLinesHokuyo);
% disp('Rotation and translation residuals from camera to Hokuyo estimation')
% disp(eRot')
% disp(eTransl')

%% Test Hokuyo to Velodyne transformation estimation
GHokuyoVelodyneEst=poseEstimationFromNormalsPointsG(normalPlanesVelodyneFlat,distancePlanesVelodyneFlat,XLinesHokuyoFlat2D);
disp('[GHokuyoVelodyneEst GHokuyoVelodyne GHokuyoVelodyneEst-GHokuyoVelodyne]')
disp([GHokuyoVelodyneEst GHokuyoVelodyne GHokuyoVelodyneEst-GHokuyoVelodyne])

%compute and show errors
e=poseEstimationFromNormalsPointsResidualsG(GHokuyoVelodyneEst,normalPlanesVelodyneFlat,distancePlanesVelodyneFlat,XLinesHokuyoFlat2D);
disp('Residuals from camera to Hokuyo estimation')
disp(e)
% GHokuyoVelodyneEst=poseEstimationFromNormalsDistancesPointsG(normalPlanesVelodyne,distancePlanesVelodyne,XLinesHokuyo);
% disp('[GHokuyoVelodyneEst GHokuyoVelodyne GHokuyoVelodyneEst-GHokuyoVelodyne]')
% disp([GHokuyoVelodyneEst GHokuyoVelodyne GHokuyoVelodyneEst-GHokuyoVelodyne])

%% Prepare dataset for averaging
% Prepare network
poses.GCameraVelodyne=GCameraVelodyneEst;
poses.SigmaCameraVelodyne=eye(6);
poses.GHokuyoCamera=GHokuyoCamera;
poses.SigmaHokuyoCamera=eye(6);
poses.GHokuyoVelodyne=GHokuyoVelodyne;
poses.SigmaHokuyoVelodyne=eye(6);

[t_node,idxCamera,idxVelodyne,idxHokuyo]=multiCalibrationPrepareNetwork(poses,false);

% Add ground truth
NNodes=3;
giTruth=zeros(4,4,NNodes);
for iNode=1:NNodes
    switch iNode
        case idxCamera
            g=GCamera;
        case idxVelodyne
            g=GVelodyne;
        case idxHokuyo
            g=GHokuyo;
        otherwise
            error('NNodes not valid');
    end
    giTruth(:,:,iNode)=g;
end
t_node=testNetworkAddGroundTruth(t_node,giTruth,'methodAbsolutePoses',methodAbsolutePoses);
% Average
optsLocalization={'displayIt','showMessages'};
t_node=localization_MLE_rigid(t_node,optsLocalization{:},'noinit');
% Compute poses
testNetworkDisplayErrors(t_node,'rt','methodAbsolutePoses',methodAbsolutePoses)

%%

function e=computePointLineErrors(XLines,VectorLines)
NPoints=size(XLines,2);
NLines=size(VectorLines,3);
e=zeros(2,NPoints,NLines);
for iPlane=1:NLines
    e(:,:,iPlane)=VectorLines(:,:,iPlane)'*homogeneous(XLines(:,:,iPlane),4);
end
e=reshape(e,2,[]);


function e=computePointPlaneErrors(XLines,VectorPlanes)
NPoints=size(XLines,2);
NLines=size(VectorPlanes,2);
e=zeros(1,NPoints,NLines);
for iPlane=1:NLines
    e(:,:,iPlane)=VectorPlanes(:,iPlane)'*homogeneous(XLines(:,:,iPlane),4);
end
e=reshape(e,2,[]);


