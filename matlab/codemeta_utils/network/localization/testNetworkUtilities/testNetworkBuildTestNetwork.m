function t_node=testNetworkBuildTestNetwork(varargin)
methodInit='truth';
varNoisyTruth=0.5;
methodAbsolutePoses='reference';
structType='single';
N=8;
flagAdjacencyProvided=false;

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'references'
            methodAbsolutePoses='reference';
        case 'poses'
            methodAbsolutePoses='pose';
        case 'methodabsoluteposes'
            ivarargin=ivarargin+1;
            methodAbsolutePoses=lower(varargin{ivarargin});
        case 'methodinit'
            ivarargin=ivarargin+1;
            methodInit=lower(varargin{ivarargin});
        case 'varnoisytruth'
            ivarargin=ivarargin+1;
            varNoisyTruth=lower(varargin{ivarargin});
        case 'n'
            ivarargin=ivarargin+1;
            N=lower(varargin{ivarargin});
        case 'a'
            ivarargin=ivarargin+1;
            flagAdjacencyProvided=true;
            A=lower(varargin{ivarargin});
        otherwise
                error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

if length(varNoisyTruth)==1
    varNoisyTruth=[1 1]*varNoisyTruth;
end

if flagAdjacencyProvided
    N=size(A,1);
else
    A=adjgallery(N,'banded',2);
end
[G,X]=testNetworkCreateAbsolutePoses(N);

t_node=testNetworkCreateStruct(A,'Type',structType);  
t_node=testNetworkAddGroundTruth(t_node,G,'methodAbsolutePoses',methodAbsolutePoses,'flagInvertG',true);
t_node=testNetworkProjectImages(t_node,X,'methodAbsolutePoses',methodAbsolutePoses);
t_node=testNetworkAddMeasurements(t_node,'Method','Essential');
switch methodInit
    case 'truth'
        t_node=testNetworkInitializeStates(t_node,'MethodR','truth',...
            'MethodT','truth','MethodScale','truth');
    case 'noisytruth'
        t_node=testNetworkInitializeStates(t_node,'MethodR','NoisyTruth',...
            'MethodT','NoisyTruth','MethodScale','noisyTruth',...
            'SigmaR',varNoisyTruth(1),'SigmaT',varNoisyTruth(2));
    case 'rand'
        t_node=testNetworkInitializeStates(t_node,'MethodR','rand',...
            'MethodT','rand','MethodScale','noscales');
end
t_node=testNetworkAddEssentialMatricesGroundTruth(t_node,'methodAbsolutePoses',methodAbsolutePoses);
t_node=testNetworkAddMeasurementsEssential(t_node,'method','truth');
t_node=splitgi(t_node);
t_node=splitgij(t_node);
t_node=splitgi(t_node,'gitruth','Ritruth','Titruth');
t_node=splitgij(t_node,'gijtruth','Rijtruth','Tijtruth');
