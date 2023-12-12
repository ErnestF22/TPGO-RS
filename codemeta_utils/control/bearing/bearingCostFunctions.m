%Base functions for the bearing localization cost functions
%funs=bearingCostFunctions(name,parameter)
%Inputs
%   name    type of the cost. With c being the inner product between the
%           desired and measured bearings, available choices are 
%       'cosine'      f(c)=(c-1), cosine of the angle
%       'angle'       f(c)=acos(c), cosine of the angle
%       'cosineSq'    f(c)=(c-1)^2/2, square of the cosine of the angle
%       'angleSq'     f(c)=acos(c)^2/2, square of the angle
%       'power',name2,k     f(c)=f2(c)^k, where f2(c) is any of the other
%                           options and k is a number
%       'barrier',cs,cc     barrier function that is zero from 1 to cs and
%                           then goes to infinity at cc
%Outputs
%   funs    struct with fields that are function handles to:
%       funs.f      cost function
%       funs.df     first derivative of the cost function
%       funs.ddf    second derivative of the cost function
%
%Note:  the final cost used for the bearing localization will be 
%       sum(di.*funs.f(ci)), where di contains the distance to the i-th
%       landmark 
function funs=bearingCostFunctions(name,varargin)
NParameters=length(varargin);
switch lower(name)
    case 'power'
        name2=varargin{1};
        k=varargin{2};
        funs2=bearingCostFunctions(name2);
        funs.f=@(c) funs2.f(c)^k;
        funs.df=@(c) k*funs2.f(c)^(k-1)*funs2.df(c);
        funs.ddf=@(c) k*(k-1)*funs2.f(c)^(k-2)*funs2.df(c)^2 ...
            +k*funs2.f(c)^(k-1)*funs2.ddf(c);
    case 'cosine'
        funs.f=@(c) 1-c;
        funs.df=@(c) -ones(size(c));
        funs.ddf=@(c) zeros(size(c));
    case 'cosinehalf'
        funs.f=@(c) (1-c)/2;
        funs.df=@(c) -ones(size(c))/2;
        funs.ddf=@(c) zeros(size(c));
    case 'angle'
        funs.f=@(c) abs(acos(c));
        funs.df=@(c) -(1-c^2)^-0.5;
        funs.ddf=@(c) -c*(1-c^2)^-1.5;
    case 'cosinesq'
        funs.f=@(c) (c-1).^2/2;
        funs.df=@(c) c-1;
        funs.ddf=@(c) ones(size(c));
    case 'anglesq'
        funs.f=@(c) acos(c).^2/2;
        funs.df=@(c) dfanglesq(c);
        funs.ddf=@(c) ddfangle(c);
    case 'cosinel1l2'
        if NParameters==0
            b=0.01;
        else
            b=varargin{1};
        end
        a=2*sqrt(0.5/b^2 + 1) - 2;
        funs.f=@(c) (2*sqrt(0.5*(c/2-0.5).^2./b.^2.+1)-2)/a;
        funs.df=@(c) (c/2 - 1/2)/(2*a*b^2*((c/2 - 1/2)^2/(2*b^2) + 1)^(1/2));
        funs.ddf=@(c) 1/(4*a*b^2*((c/2 - 1/2)^2/(2*b^2) + 1)^(3/2));
    case 'anglel1l2'
        if NParameters==0
            b=0.1;
        else
            b=varargin{1};
        end
        a=2*sqrt(0.5/b^2 + 1) - 2;
        funs.f=@(c) (2*((acos(c).^2/(2*b^2) + 1).^(1/2) - 1))/a;
        funs.df=@(c) -acos(c)./(a*b^2*(1 - c.^2).^(1/2)*(acos(c).^2/(2*b^2) + 1).^(1/2));
        funs.ddf=@(c) (2*b^2*(1 - c.^2).^(3/2) + c.^3.*acos(c).^3 - c.*acos(c).^3 - 2*b^2*c.*acos(c) + 2*b^2*c.^3.*acos(c))/(2*a*b^4*(1 - c.^2).^(5/2)*(acos(c).^2/(2*b^2) + 1).^(3/2));
    case 'squared'
        funs.f=@(e) 0.5*e.^2;
        funs.df=@(e) e;
        funs.ddf=@(e) 1;
    case 'barrier'
        if NParameters==0
            cs=0.5;
        else
            cs=varargin{1};
        end
        if NParameters<2
            cc=-1;
        else
            cc=varargin{2};
        end
        if cc>cs
            error('cc must be smaller than cs')
        end
        
        funs.f=@(c) fbarrier(c,cs,cc);
        funs.df=@(c) dfbarrier(c,cs,cc);
        funs.ddf=@(c) ddfbarrier(c,cs,cc);
    otherwise
        error('Name %s not valid',name)
end

function df=dfanglesq(c)
if c==1
    df=-1;
else
    df=-acos(c)./sqrt(1-c.^2);
end

function ddf=ddfangle(c)
if c==1
    ddf=1/3;
else
    ddf=1./(1-c.^2)-(c.*acos(c))./(1-c.^2).^(3/2);
end


function f=fbarrier(c,cs,cc)
flagOut=c<=cc;
flagIn=c>cc & c<cs;
f=zeros(size(c));
f(flagOut)=Inf;
f(flagIn)=((c(flagIn)-cs).^2./(c(flagIn)-cc));

function df=dfbarrier(c,cs,cc)
flagOut=c<=cc;
flagIn=c>cc & c<cs;
df=zeros(size(c));
df(flagOut)=Inf;
df(flagIn)=-(cs-c(flagIn)).*(cs-2*cc+c(flagIn))./(cc-c(flagIn))^2;

function ddf=ddfbarrier(c,cs,cc)
flagOut=c<=cc;
flagIn=c>cc & c<cs;
ddf=zeros(size(c));
ddf(flagOut)=Inf;
ddf(flagIn)=(-2*(cs-cc)^2)./(cc-c(flagIn)).^3;
