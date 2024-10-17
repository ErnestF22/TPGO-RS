% This function computes the finite difference. It is possible to choose
% between the central method or the forward method. By default it is using
% the central method.
% You can look here for more details:
% https://en.wikipedia.org/wiki/Finite_difference
function df=funApproxDer(f,t,h,sz)
method='central';

if(nargin<2)
    t=0;
end

if(nargin<3)
    h=1e-6;
end

if(nargin<4)
    sz=size(f(t(1)));
    if(max(sz)==1)
        sz=1;
    end
end

switch method
    case 'forward'
        dfunc=@fdiff;
    case 'central'
        dfunc=@cdiff;
end

Np=length(t);
df=zeros([sz Np]);
for(p=1:Np)
    switch length(sz)
        case 1
            if(sz==1)
                df(p)=dfunc(p);
            else
                df(:,p)=dfunc(p);
            end
        case 2
            df(:,:,p)=dfunc(p);
    end
end

    function fd=fdiff(p)
        fd=(f(t(p)+h)-f(t(p)))/h;
    end

    function cd=cdiff(p)
        cd=(f(t(p)+h)-f(t(p)-h))/(2*h);
    end

end