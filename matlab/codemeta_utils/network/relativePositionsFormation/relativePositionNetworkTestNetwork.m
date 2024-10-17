function t_node=relativePositionNetworkTestNetwork(NNodes,varargin)
flagEdgesProvided=false;
flagEdgesControl=false;

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'edges'
            ivarargin=ivarargin+1;
            E=varargin{ivarargin};
            flagEdgesProvided=true;
        case 'edgescontrol'
            ivarargin=ivarargin+1;
            EControl=varargin{ivarargin};
            flagEdgesControl=true;
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

if ~flagEdgesProvided
    A=adjgallery(NNodes,'banded',2);
    E=adj2edges(A);
end

t_node=testNetworkCreateStruct(E,'type','single','Edges',NNodes);

theta=0:2*pi/NNodes:2*pi*(NNodes-1)/NNodes;
TiTruth=6*[cos(theta); sin(theta)];
t_node.TiTruth=TiTruth;
t_node.TijTruth=relativePositionNetworkCompute(TiTruth,E);

if flagEdgesControl
    t_node.EControl=EControl;
    t_node.TijControlTruth=relativePositionNetworkCompute(TiTruth,EControl);
end