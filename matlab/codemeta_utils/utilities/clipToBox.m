%Rescale points so that they fit in a box
function x=clipToBox(x,a,varargin)
flagExtend=false;

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'extend'
            flagExtend=true;
        otherwise
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

dimX=size(x,1);
if length(a)~=dimX
    a=a*ones(1,dimX);
end
m=max(diag(a)\abs(x),[],1);
if ~flagExtend
    m=max(m,1);
end
x=x./([1;1]*m);
