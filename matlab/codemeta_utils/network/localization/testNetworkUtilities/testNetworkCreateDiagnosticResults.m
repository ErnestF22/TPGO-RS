%function t_node=testNetworkCreateFakeResults(t_node)
%Create a set of estimated poses t_node().gi by using a transformed, noisy
%version of the ground truth t_node().gitruth.
%This function is for debugging and testing purposes only.

%%AUTORIGHTS%%

function t_node=testNetworkCreateDiagnosticResults(t_node,varargin)

sigmaGlobalR=1000;
sigmaGlobalT=3;
sigmaGlobalScale=0.5;
sigmaLocalR=0.2;
sigmaLocalT=0.5;

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'sigmaglobalr'
            ivarargin=ivarargin+1;
            sigmaGlobalR=varargin{ivarargin};
        case 'sigmaglobalt'
            ivarargin=ivarargin+1;
            sigmaGlobalT=varargin{ivarargin};
        case 'sigmaglobalscale'
            ivarargin=ivarargin+1;
            sigmaGlobalScale=varargin{ivarargin};
        case 'sigmalocalr'
            ivarargin=ivarargin+1;
            sigmaLocalR=varargin{ivarargin};
        case 'sigmalocalt'
            ivarargin=ivarargin+1;
            sigmaLocalT=varargin{ivarargin};
        otherwise
                error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

N=length(t_node);

gGlobal=noiserigid(eye(4),sigmaGlobalR,sigmaGlobalT);
scaleGlobal=1-sigmaGlobalScale*log(rand);

for iNode=1:N
    t_node(iNode).gi=t_node(iNode).gitruth;
    t_node(iNode).gi=noiserigid(t_node(iNode).gitruth,sigmaLocalR,sigmaLocalT);
    t_node(iNode).gi(1:3,4)=scaleGlobal*t_node(iNode).gi(1:3,4);
    t_node(iNode).gi=gGlobal*t_node(iNode).gi;
end
