%Convert points to homogeneous coordinates
%function x=homogeneous(x,dDesired)
%If size(x,2)==dDesired-1, adds a row of all ones to x.
%If size(x,2)==dDesired, divides columns by the last row
%If size(x,2)==dDesired+1, divides columns by the last row and eliminate it
%For any other value of dDesired, throws an error
%If dDesired is omitted, default to size(x,2)+1
function x=homogeneous(x,dDesired,varargin)
sz=size(x);
d=sz(1);
if ~exist('dDesired','var') || isempty(dDesired)
    dDesired=d+1;
end

dataType='points';
%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case {'points','velocities'}
            dataType=lower(varargin{ivarargin});
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

switch dataType
    case 'velocities'
        switch d
            case dDesired-1
                x=[x;zeros([1 sz(2:end)])];
            case dDesired
                if any(x(end,:)~=0)
                    warning('Velocities already in homogeneous coordinates, but last entry is not zero. Setting it to zero');
                    x(end,:)=0;
                end
            otherwise
                error('Dimension of x and dDesired are not compatible')
        end
    case 'points'
        switch d
            case dDesired-1
                x=[x;ones([1 sz(2:end)])];
            case dDesired
                x=x./repmat(reshape(x(d,:),[1 sz(2:end)]),d,1);
            case dDesired+1
                x=x./repmat(reshape(x(d,:),[1 sz(2:end)]),d,1);
                x=reshape(x(1:dDesired,:),[dDesired sz(2:end)]);
            otherwise
                error('Dimension of x and dDesired are not compatible')
        end
        
end
