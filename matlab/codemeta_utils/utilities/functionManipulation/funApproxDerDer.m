%function df=approx_derder(f,t,h,sz)
%Compute the approximated second derivative of f on a grid given by t. The optional arguments h is the width for the finite difference and sz is the size of the output of the function f
function df=approx_derder(f,t,h,sz)
if(nargin<3)
    h=1e-6;
end
if(nargin<4)
    sz=size(f(t(1)));
    if(max(sz)==1)
        sz=1;
    end
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

    function cd=dfunc(p)
        cd=(f(t(p)+h)-2*f(t(p))+f(t(p)-h))/h/h;
    end

end
