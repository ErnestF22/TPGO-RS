%function h=approx_RiemGrad(lf,f,x)
%Computes the Riemannian gradient of the function f at x using finite
%differences.
function h=approx_RiemGrad(lf,f,x)
    h=zeros(size(x));
    for n=1:size(x,3)
        T=lf.tangentBasis(x(:,:,n));
        for d=1:size(T,3)
            v=zeros(size(x));
            v(:,:,n)=T(:,:,d);
            fdir=@(t) evalFunGeod(lf,f,x,t*v);
            h(:,:,n)=h(:,:,n)+approx_der(fdir,0)*T(:,:,d);
        end
    end
end

%evaluate a function along a geodesic
function f=evalFunGeod(lf,fun,y0,h)
    y=multiexp(lf,y0,h);
    f=fun(y);
end

function y=multiexp(lf,y0,h)
    y=zeros(size(y0));
    for n=1:size(y0,3)
        y(:,:,n)=lf.exp(y0(:,:,n),h(:,:,n));
    end
end        

