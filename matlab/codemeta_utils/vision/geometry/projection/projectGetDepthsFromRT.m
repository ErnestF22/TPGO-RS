%Get the depths (lambda) of points seen by the given camera
%function [lambda,JRT]=projectGetDepthsFromRT(R,T,X,varargin)
%Note: JRT is always given in terms of the derivatives of R,T in the 'pose'
%interpretation
function [lambda,JRT]=projectGetDepthsFromRT(R,T,X,varargin)
flagComputeJacobian=nargout>1;

NPoses=size(R,3);
NX=size(X,2);

if NPoses>1
    lambda=zeros(NPoses,NX);
    if flagComputeJacobian
        JRT=zeros(NPoses,6,NX);
    end
    for iPose=1:NPoses
        if flagComputeJacobian
            [lambda(iPose,:),JRT(iPose,:,:)]=projectGetDepthsFromRT(R(:,:,iPose),T(:,iPose),X,varargin{:});
        else
            lambda(iPose,:)=projectGetDepthsFromRT(R(:,:,iPose),T(:,iPose),X,varargin{:});
        end            
    end
else
    methodAbsolutePoses='pose';
    
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
            otherwise
                error(['Argument ' varargin{ivarargin} ' not valid!'])
        end
        ivarargin=ivarargin+1;
    end

    XTransformed=rigidTransform(R,T,X,'methodAbsolutePoses',methodAbsolutePoses);
    lambda=XTransformed(3,:);

    flagInvertG=false;
    switch methodAbsolutePoses
        case 'reference'
            flagInvertG=~flagInvertG;
        case 'pose'
            %leave flagInvertG the same
    end
    if flagInvertG
        R=invR(R);
    end
    if flagComputeJacobian
        ez=[0;0;1];
        ezRXHat=multiprod(ez',multiprod(R,hat3(X)));
        JRT=[ezRXHat repmat(ez',[1 1 size(X,2)])];
    end
end
