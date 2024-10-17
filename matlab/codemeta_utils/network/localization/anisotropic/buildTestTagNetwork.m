function [t_node,A,G]=buildTestTagNetwork(varargin)
seed=3;
L=15;   %size of square where tags are positioned
h=10;   %height of cameras
NCameras=2;
NTags=8;
flagConnectCameras=false;

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(varargin{ivarargin})
        case 'seed'
            ivarargin=ivarargin+1;
            seed=varargin{ivarargin};
        case 'ncameras'
            ivarargin=ivarargin+1;
            NCameras=varargin{ivarargin};
        case 'ntags'
            ivarargin=ivarargin+1;
            NTags=varargin{ivarargin};
        case 'flagconnectcameras'
            ivarargin=ivarargin+1;
            flagConnectCameras=varargin{ivarargin};
        otherwise
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

if ~isempty(seed)
    resetRands(seed)
end
N=NCameras+NTags;
idxCameras=1:NCameras;
idxTags=NCameras+(1:NTags);

%adjacency matrix
A=zeros(N);
if flagConnectCameras
    for iCamera=1:NCameras-1
        A(idxCameras(iCamera),idxCameras(iCamera+1))=1;
    end
end
for iCamera=1:NCameras
    for iTag=1:NTags
        A(idxCameras(iCamera),idxTags(iTag))=1;
    end
end
A=A+A';

%poses
RCamera=[
    1  0  0;
    0 -1  0;
    0  0 -1
    ];
RTags=eye(3);
TCameras=[L*rand(2,NCameras);h*ones(1,NCameras)];
TTags=[L*rand(2,NTags); zeros(1,NTags)];

G=zeros(4,4,N);
G(4,4,:)=1;
for iCamera=1:NCameras
    G(1:3,1:3,idxCameras(iCamera))=RCamera*randXYRotation();
    G(1:3,4,idxCameras(iCamera))=TCameras(:,iCamera);
end
for iTag=1:NTags
    G(1:3,1:3,idxTags(iTag))=RTags*randXYRotation();
    G(1:3,4,idxTags(iTag))=TTags(:,iTag);
end

%structure
t_node=testNetworkCreateStruct(A,'type','single');
t_node=testNetworkAddGroundTruth(t_node,G);

%Relative poses
%t_node=testNetworkAddDispersionMatricesRT(t_node,'varT',0.2*[0.1 1],'varR',0.2*[0.1 1]);
%t_node=testNetworkAddDispersionMatricesRT(t_node,'methodT','zero','methodR','identity','methodCoupling','zero');
t_node=testNetworkAddDispersionMatricesRT(t_node);
%t_node=testNetworkAddMeasurements(t_node,'method','truth','unnormalized');
t_node=testNetworkAddMeasurements(t_node,'method','noisy','joint','unnormalized');

%Essential matrices
% t_node=testNetworkAddEssentialMatricesGroundTruth(t_node);
% t_node=testNetworkAddDispersionMatricesE(t_node);
% %t_node=testNetworkAddMeasurementsEssential(t_node,'method','noisy','sigmaNoise',[]);
% t_node=testNetworkAddMeasurementsEssential(t_node,'method','truth');


%State initialization
t_node=testNetworkInitializeStates(t_node,'MethodR','Truth','MethodT','Truth','MethodScale','Ignore');
%t_node=testNetworkInitializeStates(t_node,'MethodR','Rand','MethodT','Rand','MethodScale','Ignore');


function R=randXYRotation()
theta=rand*2*pi;
R=[cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];