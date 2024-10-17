%Add dispersion matrices for the generation of relative measurements
%       'GivenR',dispMats   dispMats contains the pairwise dispersion matrices
%                           for rotations. If t_node is of type 'array',
%                           dispMats must  be a [NCameras x NCameras] cell
%                           array. If t_node is of type 'single', dispMats
%                           must be a [3 x 3 x NEdges] array.
%       'GivenT',dispMats   Similar to 'GivenR', but for the rotations
%       'Given',dispMats    Similar to the above, but for both rotations
%                           and translations.
function t_node=testNetworkAddDispersionMatricesRT(t_node,varargin)
methodR='anisotropic';
methodT='anisotropic';
methodCoupling='randn';

varR=[0.01; 0.5];   %0.01rad=0.5deg, 0.5rad=28.6deg
varT=[0.1; 2];
varCoupling=[100 200];

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'givenr'
            methodR='given';
            ivarargin=ivarargin+1;
            dispMatR=varargin{ivarargin};
        case 'givent'
            methodT='given';
            ivarargin=ivarargin+1;
            dispMatT=varargin{ivarargin};
        case 'givencoupling'
            methodCoupling='given';
            ivarargin=ivarargin+1;
            dispMatCoupling=varargin{ivarargin};
        case 'given'
            methodR='fromGiven';
            methodT='fromGiven';
            methodCoupling='fromGiven';
            ivarargin=ivarargin+1;
            dispMat=varargin{ivarargin};
        case 'givenisotropic'
            %substitute each dispersion matrix with identity times the
            %average of the evals of the original dispersion matrix
            methodR='fromGiven';
            methodT='fromGiven';
            methodCoupling='fromGiven';
            ivarargin=ivarargin+1;
            dispMat=varargin{ivarargin};
            dDispMat=size(dispMat,1);
            for iEdge=1:size(dispMat,3)
                dispMat(:,:,iEdge)=trace(dispMat(:,:,iEdge))/dDispMat*eye(dDispMat);
            end
        case 'identity'
            %substitute each dispersion matrix with identity times the
            %average of the evals of the original dispersion matrix
            methodR='fromGiven';
            methodT='fromGiven';
            methodCoupling='fromGiven';
            NEdges=testNetworkGetNumberOfEdges(t_node);
            dispMat=repmat(eye(6),[1 1 NEdges]);
        case 'methodt'
            ivarargin=ivarargin+1;
            methodT=lower(varargin{ivarargin});
        case 'methodr'
            ivarargin=ivarargin+1;
            methodR=lower(varargin{ivarargin});
        case 'methodcoupling'
            ivarargin=ivarargin+1;
            methodCoupling=lower(varargin{ivarargin});
        case 'vart'
            ivarargin=ivarargin+1;
            varT=lower(varargin{ivarargin});
        case 'varr'
            ivarargin=ivarargin+1;
            varR=lower(varargin{ivarargin});
        case 'varcoupling'
            ivarargin=ivarargin+1;
            varCoupling=lower(varargin{ivarargin});
        otherwise
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

E=testNetworkGetEdges(t_node);
NEdges=size(E,1);

%generate dispersion matrices if not provided
switch lower(methodR)
    case 'given'
        %nothing to do here
    case 'fromgiven'
        dispMatR=dispMat(1:3,1:3,:);
    otherwise
        dispMatR=generateDispersionMatrices(NEdges,methodR,varR);
end
switch lower(methodT)
    case 'given'
        %nothing to do here
    case 'fromgiven'
        dispMatT=dispMat(4:6,4:6,:);
    otherwise
        dispMatT=generateDispersionMatrices(NEdges,methodT,varT);
end
switch lower(methodCoupling)
    case 'given'
        %nothing to do here
    case 'zero'
        dispMatCoupling=zeros(3,3,NEdges);
    case 'randn'
        dispersion=fliplr(1./varCoupling);
        s=dispersion(1)+rand*(dispersion(2)-dispersion(1));
        dispMatCoupling=s*randn(3,3,NEdges);
    case 'fromgiven'
        dispMatCoupling=dispMat(1:3,4:6,:);
end

structType=testNetworkDetectStructType(t_node);

switch structType
    case 'single'
        %copy data into structure
        t_node.dispersionMatR=dispMatR;
        t_node.dispersionMatT=dispMatT;
        t_node.dispersionMat=[dispMatR dispMatCoupling; permute(dispMatCoupling,[2 1 3]) dispMatT];
        
    case 'array'
        N=length(t_node);
        [t_node.dispersionMatR]=deal(zeros(3,3,N));
        [t_node.dispersionMatT]=deal(zeros(3,3,N));
        
        for iEdge=1:NEdges
            iNode=E(iEdge,1);
            jNode=E(iEdge,2);
            
            t_node(iNode).dispersionMatR(:,:,jNode)=dispMatR(:,:,iEdge);
            t_node(iNode).dispersionMatT(:,:,jNode)=dispMatT(:,:,iEdge);
            t_node(iNode).dispersionMat(:,:,jNode)=[dispMatR(:,:,iEdge) dispMatCoupling(:,:,iEdge); dispMatCoupling(:,:,iEdge)' dispMatT(:,:,iEdge)];
        end
end
