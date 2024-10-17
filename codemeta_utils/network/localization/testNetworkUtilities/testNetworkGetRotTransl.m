function [Ri,Ti]=testNetworkGetRotTransl(t_node,varargin)

fieldName='gi';
flagInvertG=false;

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'fieldname'
            ivarargin=ivarargin+1;
            fieldName=varargin{ivarargin};
        case 'flaginvertg'
            ivarargin=ivarargin+1;
            flagInvertG=varargin{ivarargin};
        otherwise
            error('MATLAB:ArgumentInvalid',...
                ['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

structType=testNetworkDetectStructType(t_node);

switch structType
    case 'single'
        Gi=t_node.(fieldName);
    case 'array'
        Gi=cat(3,t_node.(fieldName));
end

if flagInvertG
    Gi=invg(Gi);
end

Ri=Gi(1:3,1:3,:);
Ti=squeeze(Gi(1:3,4,:));
