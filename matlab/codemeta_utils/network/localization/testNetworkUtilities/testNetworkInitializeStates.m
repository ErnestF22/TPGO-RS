%function t_node=testNetworkInitializeStates(t_node,varargin)
%Initializes the fields t_node().gi
%Optional Arguments
%   'MethodR'       define which method to use for relative rotations
%       'Truth'         Computed from G
%       'Noisy'         Noisy version of G
%       'Eye'           Identity
%       'Rand'          Uniform random
%
%   'R',R           use content of the 3x3xN matrix R
%
%   'MethodT'       define which method to use for translations
%       'Truth'         Computed from G
%       'Noisy'         Noisy version of G
%       'Zero'          Zero
%       'Rand'          Random (zero mean, unit variance Gaussian)
%
%   'T',T       Use content of the 3xN or 3x1xN matrix T
%   'G',G       Use content of the 4x4xN matrix G
%
%   'MethodScale'   define which method to use for scales
%       'Truth'         Computed from G
%       'Noisy'         Noisy version of ground truth
%       'NoScales'      Set to 1.1
%       'Rand'          Random (using 1/exp(rand)
%       'Ignore'        Do not add the lambdaij member to the struct
%
%   'Scale',S   Use content of the NxN matrix S
%
%   'SigmaT',S      Std of noise for translation with 'Noisy' (default: 0.1)
%
%   'SigmaR',S      Std of noise for rotations with 'Noisy' (default: 0.1)
%
%   'SigmaScale',S  Std of noise for scales with 'Noisy' (default: 0.1)
%       
% (*) NOTE: for any 'MethodScale' other than 'None' or 'Unnormalized', the
%     translations of t_node().gij will be replaced by their normalized versions

%%AUTORIGHTS%%

function t_node=testNetworkInitializeStates(t_node,varargin)
initMethodR='truth';
initMethodT='truth';
initMethodScale='noscales';
sigmaNoiseR=0.1;
sigmaNoiseT=0.1;
sigmaNoiseScale=0.1;
S=[];

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'methodr'
            ivarargin=ivarargin+1;
            initMethodR=varargin{ivarargin};
        case 'methodt'
            ivarargin=ivarargin+1;
            initMethodT=varargin{ivarargin};
        case 'methodscale'
            ivarargin=ivarargin+1;
            initMethodScale=varargin{ivarargin};
        case 'r'
            ivarargin=ivarargin+1;
            R=varargin{ivarargin};
            initMethodR='given';
        case 't'
            ivarargin=ivarargin+1;
            T=varargin{ivarargin};
            if length(size(T))>2
                T=squeeze(T);
            end
            initMethodT='given';
        case 'g'
            ivarargin=ivarargin+1;
            G=varargin{ivarargin};
            R=G2R(G);
            initMethodR='given';
            T=G2T(G);
            initMethodT='given';
        case 'scale'
            ivarargin=ivarargin+1;
            S=varargin{ivarargin};
            initMethodScale='given';
        case 'sigmar'
            ivarargin=ivarargin+1;
            sigmaNoiseR=varargin{ivarargin};
        case 'sigmat'
            ivarargin=ivarargin+1;
            sigmaNoiseT=varargin{ivarargin};
        case 'sigmascale'
            ivarargin=ivarargin+1;
            sigmaNoiseScale=varargin{ivarargin};
        otherwise
                error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end


structType=testNetworkDetectStructType(t_node);
N=testNetworkGetNumberOfNodes(t_node);

%allocate pose field
switch structType
    case 'single'
        t_node.gi=repmat(eye(4),[1 1 N]);
    case 'array'
        [t_node.gi]=deal(eye(4));
        [t_node.lambdaij]=deal(zeros(1,N));
end


for iNode=1:N
    gi=eye(4);
    
    %get truth pose in case we need it and if available
    if isfield(t_node,'gitruth')
        switch structType
            case 'single'
                gitruth=t_node.gitruth(:,:,iNode);
            case 'array'
                gitruth=t_node(iNode).gitruth;
        end
    end
    
    %init rotation
    switch(lower(initMethodR))
        case 'truth'
            gi(1:3,1:3)=gitruth(1:3,1:3);
        case {'noisy','noisytruth'}
            gi(1:3,1:3)=noiserot(gitruth(1:3,1:3),sigmaNoiseR);
        case 'eye'
            gi(1:3,1:3)=eye(3);
        case 'rand'
            gi(1:3,1:3)=unif_random_rot(1);
        case 'given'
            gi(1:3,1:3)=R(:,:,iNode);
        otherwise
            error('initMethodR not valid');
    end
    
    %init translation
    switch(lower(initMethodT))
        case 'truth'
            gi(1:3,4)=gitruth(1:3,4);
        case {'noisy','noisytruth'}
            gi(1:3,4)=gitruth(1:3,4)+sigmaNoiseT*randn(3,1);
        case 'zero'
            gi(1:3,4)=zeros(3,1);
        case 'rand'
            gi(1:3,4)=8*randn(3,1);
        case 'given'
            gi(1:3,4)=T(:,iNode);
        otherwise
            error('initMethodT not valid');
    end
    
    %put generated pose
    switch structType
        case 'single'
            t_node.gi(:,:,iNode)=gi;
        case 'array'
            t_node(iNode).gi=gi;
    end
    
end
 
if ~strcmpi(initMethodScale,'ignore')
    switch structType
        case 'single'
            t_node=testNetworkExecFunOnEdges(t_node,...
                @(t,i) edgeFunInitLambdaSingle(t,i,initMethodScale,S,sigmaNoiseScale));
        case 'array'
            t_node=testNetworkExecFunOnEdges(t_node,...
                @(t,i,j) edgeFunInitLambdaArray(t,i,j,initMethodScale,S,sigmaNoiseScale));
    end
end

%callbacks for testNetworkExecFunOnEdges
function t_node=edgeFunInitLambdaSingle(t_node,iEdge,initMethodScale,S,sigmaNoiseScale)
switch(lower(initMethodScale))
    case 'truth'
        [tij,lij]=cnormalize(t_node.gijtruth(1:3,4,iEdge));
        t_node.lambdaij(iEdge)=lij;
    case 'realscales'
        t_node(iNode).lambdaij(jNode)=1;
    case {'noisy','noisytruth'}
        [tij,lij]=cnormalize(t_node.gijtruth(1:3,4,iEdge));
        t_node.lambdaij(iEdge)=lij+sigmaNoiseScale*rand;
    case 'noscales'
        t_node.lambdaij(iEdge)=1.1;
    case 'given'
        t_node.lambdaij(iEdge)=S(iEdge);
    case 'rand'
        t_node.lambdaij(iEdge)=1/exp(rand);
    otherwise
        error('initMethodScales not valid');
end

function t_node=edgeFunInitLambdaArray(t_node,iNode,jNode,initMethodScale,S,sigmaNoiseScale)
switch(lower(initMethodScale))
    case 'truth'
        [tij,lij]=cnormalize(t_node(iNode).gijtruth(1:3,4,jNode));
        t_node(iNode).lambdaij(jNode)=lij;
    case 'realscales'
        t_node(iNode).lambdaij(jNode)=1;
    case {'noisy','noisytruth'}
        [tij,lij]=cnormalize(t_node(iNode).gijtruth(1:3,4,jNode));
        t_node(iNode).lambdaij(jNode)=lij+sigmaNoiseScale*rand;
    case 'noscales'
        t_node(iNode).lambdaij(jNode)=1.1;
    case 'given'
        t_node(iNode).lambdaij(jNode)=S(iNode,jNode);
    case 'rand'
        t_node(iNode).lambdaij(jNode)=1/exp(rand);
    otherwise
        error('initMethodScales not valid');
end
