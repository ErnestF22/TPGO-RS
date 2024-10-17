function d=sphere_dist(yi,yj,varargin)
if(size(yi,2)==1)
    yi=permute(yi,[1 3 2]);
end
if(size(yj,2)==1)
    yj=permute(yj,[1 3 2]);
end
N1=size(yi,2);
N2=size(yj,2);

outputType='matrix';

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case {'vector','matrix'}
            outputType=lower(varargin{ivarargin});
        case 'outputtype'
            ivarargin=ivarargin+1;
            outputType=lower(varargin{ivarargin});
        otherwise
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

switch lower(outputType)
    case 'matrix'
        d=zeros(N2,N1);
        for iN1=1:N1
            for iN2=1:N2
                d(iN2,iN1)=vctAngle(yi(:,iN1),yj(:,iN2));
            end
        end
        
    case 'vector'
        if(N1~=N2)
            error('Length of inputs must be the same');
        end
        d=zeros(N1,1);
        for iN1=1:N1
            d(iN1)=vctAngle(yi(:,iN1),yj(:,iN1));
        end
        
    otherwise
        error('outputType must be ''matrix'' or ''vector''');
end