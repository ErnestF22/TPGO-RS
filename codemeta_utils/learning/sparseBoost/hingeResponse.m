%Evaluate response of a linear classifier with hinge loss
%r=hingeResponse(x,w,b,varargin)
%
%Inputs
%   x   [d x Nx] matrix of points at which we want to evaluate the function
%   w   [d x Nw] matrix of normals
%   b   [1 x Nw] vector of biases
%Optional inputs
%   'appendOpposite' append also response for -w'*x-b
%Outputs
%   r   [Nx x Nw] matrix of responses

function r=hingeResponse(x,w,b,varargin)
flagAppendOpposite=false;
%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'flagappendopposite'
            ivarargin=ivarargin+1;
            flagAppendOpposite=varargin{ivarargin};
        case 'appendopposite'
            flagAppendOpposite=true;
        otherwise
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end


[d,Nx]=size(x);

s=x'*w+ones(Nx,1)*b;
r=max(0,s);
if flagAppendOpposite
    r=[r max(0,-s)];
end
