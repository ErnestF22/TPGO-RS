%function [funs]=consensus_rot3_almostGlobal_functions(varargin)
%Shape functions for consensus_rot3 for almost global convergence.
%Input arguments can be either:
%   'b', b                  inverse time constant for exponential
%                           If b is a function handle, the funs will be
%                           time-dependent
%   'thetaq', theta0,q      angle of inversion of the second derivative and
%                           minimum correlation between axes (see theory)
%   'Type', type        Type of function. The following are available.
%       'Tron' (default), 'L1-L2', 'Fair'
%
%Note: the functions are valid only for b>0 and t>0.

%%AUTORIGHTS%%

% See
% http://research.microsoft.com/en-us/um/people/zhang/inria/publis/tutorial-estim/node24.html
% for more infos and references

function [funs]=consensus_rot3_almostGlobal_functions(varargin)

type='tron';

b=1;

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'b'
            ivarargin=ivarargin+1;
            b=varargin{ivarargin};
        case 'thethaq'
            ivarargin=ivarargin+1;
            t0=varargin{ivarargin};
            ivarargin=ivarargin+1;
            q=varargin{ivarargin};
            b=(q*t0)^-1;
        case 'type'
            ivarargin=ivarargin+1;
            type=varargin{ivarargin};
        otherwise
            disp('Offending optional argument')
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end


%a scales reshaped cost function such that f(pi)=0.5*pi^2  

if ~isa(b,'function_handle')
    switch lower(type)
        case 'squared'
            funs.f=@(t) 0.5*t.^2;
            funs.df=@(t) t;
            funs.dfnorm=@(t) ones(size(t));
            funs.ddf=@(t) ones(size(t));
            funs.mumax=2;
            a=1;
        case {'tron','proposed'}
            funs.f0=@(t,b) +1/b^2-(1/b^2+t/b).*exp(-b*t);

            a=0.5*pi^2/funs.f0(pi,b);  
            a1=a/b;
            a2=a/b^2;
            % funs.f=@(t) a*f0(t);
            %funs.f=@(t) a*(-(1/b^2+t/b).*exp(-b*t)+1/b^2);
            funs.f=@(t) -(a2+a1*t).*exp(-b*t)+a2;
            funs.df=@(t) a*t.*exp(-b*t);
            funs.dfnorm=@(t) a*exp(-b*t);
            funs.ddf=@(t) a*(1-b*t).*exp(-b*t);
            funs.mumax=a;

            funs.mumaxb=@(b) 0.5*pi^2/funs.f0(pi,b);
        case {'l1-l2','l1l2'}
            funs.f0=@(t,b) 2*(sqrt(1+(b*t).^2/2)-1);
            a=0.5*pi^2/funs.f0(pi,b);

            funs.f=@(t) a*funs.f0(t,b);
            funs.df=@(t) (a*b^2*t)./sqrt((b^2*t.^2)/2 + 1);
            funs.dfnorm=@(t) (a*b^2)./sqrt((b^2*t.^2)/2 + 1);
            funs.ddf=@(t) (a*b^2)./sqrt((b^2*t.^2)/2 + 1) - (a*b^4*t.^2)./(2*sqrt((b^2*t.^2)/2 + 1).^3);

            funs.mumax=a*b^2;
            funs.mumaxb=@(b) 0.5*pi^2/funs.f0(pi,b)*b^2;
        case 'fair'
            c=1.3998;
            funs.f0=@(t,b) c*b*t-c^2*log(1+b*t/c);
            a=0.5*pi^2/funs.f0(pi,b);

            funs.f=@(t) a*funs.f0(t,b);
            funs.df=@(t) a*c*t*b^2./(c + b*t);
            funs.dfnorm=@(t) a*c*b^2./(c + b*t);
            
            funs.ddf=@(t) a*c^2*b^2./(c + b*t).^2;

            funs.mumax=a*b^2;
            funs.mumaxb=@(b) 0.5*pi^2/funs.f0(pi,b)*b^2;
        case 'huber'
            c=1;
            f01=@(t,b) c^2*t.^2/2;
            f02=@(t,b) b*c^2*(t-b/2);
            funs.f0=@(t,b) selectFunThresh(...
                @(t) f01(t,b),...
                @(t) f02(t,b),...
                t,b);

            a=0.5*pi^2/funs.f0(pi,b);
            funs.f=@(t) a*funs.f0(t,b);

            df01=@(t) c^2*t;
            df02=@(t) c^2*b;
            funs.df=@(t) a*selectFunThresh(df01,df02,t,b);

            dfnorm01=@(t) c^2;
            dfnorm02=@(t) c^2*b./t;
            funs.dfnorm=@(t) a*selectFunThresh(dfnorm01,dfnorm02,t,b);

            ddf01=@(t) c^2;
            ddf02=@(t) 0;
            funs.ddf=@(t) a*selectFunThresh(ddf01,ddf02,t,b);

            funs.mumax=a*c^2;
            %If we assume b<pi:
            %funs.mumax=0.5*pi^2/(b*(t-b/2));
            funs.mumaxb=@(b) 0.5*pi^2/funs.f0(pi,b)*c^2;
        case 'cauchy'
            c=2.3849;
            funs.f0=@(t,b) c^2/2*log(1+(b*t/c).^2);

            a=0.5*pi^2/funs.f0(pi,b);
            funs.f=@(t) a*funs.f0(t,b);
            funs.df=@(t) (a*b^2*c^2*t)./(b^2*t.^2 + c^2);
            funs.dfnorm=@(t) (a*b^2*c^2)./(b^2*t.^2 + c^2);
            funs.ddf=@(t) (a*b^2*c^2*(- b^2*t.^2 + c^2))/(b^2*t.^2 + c^2)^2;

            funs.mumax=a*b^2;
            funs.mumaxb=@(b) 0.5*pi^2/funs.f0(pi,b)*b^2;
        case {'geman-mcclure','gemanmcclure'}
            funs.f0=@(t,b) (b*t).^2/2/(1+(b*t).^2);

            a=0.5*pi^2/funs.f0(pi,b);
            funs.f=@(t) a*funs.f0(t,b);
            funs.df=@(t) a*(b^2*t)./(b^2*t.^2 + 1).^2;
            funs.dfnorm=@(t) (a*b^2)./(b^2*t.^2 + 1).^2;
            funs.ddf=@(t) -a*(b^2*(3*b^2*t.^2 - 1))./(b^2*t.^2 + 1).^3;

            funs.mumax=a*b^2;
            funs.mumaxb=@(b) 0.5*pi^2/funs.f0(pi,b)*b^2;
        case 'welsh'
            c=2.9846;
            funs.f0=@(t,b) c^2/2*(1-exp(-(b*t/c).^2));

            a=0.5*pi^2/funs.f0(pi,b);
            funs.f=@(t) a*funs.f0(t,b);

            funs.df=@(t) a*b^2*t.*exp(-(b^2*t.^2)/c^2);
            funs.dfnorm=@(t) a*b^2*exp(-(b^2*t.^2)/c^2);
            funs.ddf=@(t) (a*b^2*exp(-(b^2*t.^2)/c^2)*(- 2*b^2*t.^2 + c^2))/c^2;

            funs.mumax=a*b^2;
            funs.mumaxb=@(b) 0.5*pi^2/funs.f0(pi,b)*b^2;
        otherwise
            error('Type of function is not valid');
    end                
    funs.a=a;
    funs.b=b;

else
    f0=@(t,it) -(1/b(it)^2+t/b(it)).*exp(-b(it)*t)+1/b(it)^2;

    a=@(it) 0.5*pi^2/f0(pi,it);    %scale reshaped cost function such that f(pi)=0.5*pi^2  

%         funs.f=@(t,it) a(it)*f0(t,it);
    funs.f=@(t,it) a(it)*(-(1/b(it)^2+t/b(it)).*exp(-b(it)*t)+1/b(it)^2);
    funs.df=@(t,it) a(it)*t.*exp(-b(it)*t);
    funs.dfnorm=@(t,it) a(it)*exp(-b(it)*t);
    funs.ddf=@(t,it) a(it)*(1-b(it)*t).*exp(-b(it)*t);
    funs.mumax=@(it) 2*a(it);
    funs.a=a;
    funs.b=b;
    funs.f0=f0;
end

function fevalx=selectFunThresh(f1,f2,x,thresh)
fevalx=zeros(size(x));
flagThresh=(x<=thresh);
fevalx(flagThresh)=f1(x(flagThresh));
fevalx(~flagThresh)=f2(x(~flagThresh));
