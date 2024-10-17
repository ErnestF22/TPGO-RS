function t_node=bearingNetworkBuildTestNetwork(NNodes,D,varargin)
flagEdgesProvided=false;
flagEdgesRangesProvided=false;
flagFoVs=false;

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'edges'
            ivarargin=ivarargin+1;
            E=varargin{ivarargin};
            flagEdgesProvided=true;
        case 'edgesranges'
            ivarargin=ivarargin+1;
            Er=varargin{ivarargin};
            flagEdgesRangesProvided=true;
        case 'flagfovs'
            ivarargin=ivarargin+1;
            flagFoVs=varargin{ivarargin};
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

if ~exist('NNodes','var') || isempty(NNodes)
    NNodes=7;
end

if ~exist('D','var') || isempty(D)
    D=2;
end

if ~flagEdgesProvided
    switch(D)
        case 2
            A=adjgallery(NNodes,'banded',2);
        case 3
            A=zeros(11);
            A(1,6)=1;
            A(6,11)=1;
            A(3,4)=1;
            A(4,5)=1;
            A(1,2:5)=1;
            A(6,[2 4])=1;
            for iNode=2:5
                A(iNode,iNode+5)=1;
            end
            A(7,10)=1;
            A(9,10)=1;
            A(6,[8 10])=1;
            A(11,7:10)=1;
            A=A+A';
            A=A(1:NNodes,1:NNodes);
    end
    E=adj2edges(A);
end

t_node=testNetworkCreateStruct(E,'type','single','Edges',NNodes);
if ~flagEdgesRangesProvided
    if NNodes>4
        t_node.Er=[3 4; 4 3];
    else
        t_node.Er=[1 2; 2 1];
    end
else
    t_node.Er=Er;
end

switch D
    case 2
        theta=0:2*pi/NNodes:2*pi*(NNodes-1)/NNodes;
        Titruth=6*[cos(theta); sin(theta)];
    case 3
        if NNodes>11
            error('3-D configuration for more than 11 nodes not implemented')
        end
        Titruth=6*sphereGrid('NRadii',1,'NLatitudes',4,'NLongitudes',4,'sorted');
        Titruth=rot([0;0;1],pi/6)*Titruth;
        Titruth=Titruth(:,1:NNodes);
end
t_node=bearingNetworkAddGroundTruth(t_node,Titruth);
t_node.Titruth=t_node.Titruth-mean(t_node.Titruth,2)*ones(1,NNodes);
t_node=bearingNetworkInitializeStates(t_node);
t_node=bearingNetworkAddMeasurements(t_node);

%field of views
if D==2 && flagFoVs
    t_node.fov=[160*pi/180 140*pi/180];
    t_node.y0=[1;0];
    t_node.Ritruth=fitFov(t_node,t_node.Titruth);
    t_node.Ri=fitFov(t_node,t_node.Ti);
end

%features
nodesFeatures=1:7;
t_node.Tif=8*[-1 -1 1 1; -1 1 -1 1];
NFeaturesPerNode=size(t_node.Tif,2);
t_node.Ef=[reshape(ones(NFeaturesPerNode,1)*nodesFeatures,1,[]); repmat(1:NFeaturesPerNode,1,length(nodesFeatures))]';

function R=fitFov(t_node,x)
R=bearingNetworkFitFov(x,t_node.E,t_node.y0);
