%function [xt,dxt,x0,dx0,ddxt]=real_geodFun(x0,v,varargin)
%Creates a geodesic (straight line) starting from x0 with direction v
function [xt,dxt,x0,dx0,ddxt]=real_geodFun(x0,v,varargin)
s=1;
%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'speed'
            ivarargin=ivarargin+1;
            s=varargin{ivarargin};
            if ischar(s)
                switch s
                    case 'quadratic'
                        s=@(t) t^2/2;
                        ds=@(t) t;
                        dds=@(t) 1;
                    case 'cubic'
                        s=@(t) t^3/3;
                        ds=@(t) t^2;
                        dds=@(t) 2*t;
                    case 'rand'
                        s=rand;
                    otherwise
                        error('Speed profile not recognized')
                end
            end
        otherwise
                error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end
sz=size(x0);
d=sz(1);
NDim=length(sz);
switch NDim
    case {1,2}
        NGeodesics=sz(2);
        u=ones(1,NGeodesics);
    case 3
        NGeodesics=sz(3);
        u=ones(1,1,NGeodesics);
end

%vectorize s if multiple geodesics are required
if ~isa(s,'function_handle') && length(s)==1 && NGeodesics>1
    s=s*u;
elseif isa(s,'function_handle') && length(s(0))==1 && NGeodesics>1
    s=@(t) s(t)*u;
    ds=@(t) ds(t)*u;
    dds=@(t) dds(t)*u;
end

%ensure that s, ds, and dds are function handles
if ~isa(s,'function_handle')
    dds=@(t) 0;
    ds=@(t) s;
    s=@(t) s*t; %Note: we overwrite s, so we cannot swap with line before
end

switch NDim
    case {1,2}
        xt=@(t) x0+v.*(ones(d,1)*s(t));
        dxt=@(t) v.*(ones(d,1)*ds(t));
        dx0=dxt(0);
        ddxt=@(t) v.*(ones(d,1)*dds(t));
    case 3
        xt=@(t) x0+multiprod(v,s(t));
        dxt=@(t) multiprod(v,ds(t));
        dx0=dxt(0);
        ddxt=@(t) multiprod(v,dds(t));
end        

end %file function
