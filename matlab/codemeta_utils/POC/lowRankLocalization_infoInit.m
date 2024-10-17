function WInfo=lowRankLocalization_infoInit(nbNodes,varargin)
WInfo.nbNodes=nbNodes;
WInfo.dimAmbient=3;
WInfo.flagRotationAugmented=false;
WInfo.iNodeFix=1;

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'inodefix'
            %index of node to fix
            ivarargin=ivarargin+1;
            WInfo.iNodeFix=varargin{ivarargin};
        case 'dim'
            %ambient dimension (typically 2 or 3)
            ivarargin=ivarargin+1;
            WInfo.dimAmbient=varargin{ivarargin};
        case 'rotationaugmented'
            %include columns with explicit rotations
            WInfo.flagRotationAugmented=true;
        otherwise
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end
WInfo.size=[nbNodes*WInfo.dimAmbient nbNodes];
if WInfo.flagRotationAugmented
    numCols=WInfo.size(2);
    WInfo.jNodeR=numCols+1:numCols+WInfo.dimAmbient;
    WInfo.size(2)=numCols+WInfo.dimAmbient;
end
WInfo.numel=prod(WInfo.size);
WInfo.idxMatrix=lowRankLocalization_idxMatrix(WInfo);

