function [xt,dxt,x0,dx]=sphere_randGeodFun(x0,varargin)
s=1;
%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'speed'
            ivarargin=ivarargin+1;
            s=varargin{ivarargin};
        otherwise
                error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

if ~exist('x0','var') || isempty(x0)
    x0=sphere_randn(eye(3,1));
end

dxNorm=sphere_randTangentNormVector(x0);
dx=s*dxNorm;
xt=@(t) x0*cos(s*t)+dxNorm*sin(s*t);
dxt=@(t) s*(-x0*sin(s*t)+dxNorm*cos(s*t));
