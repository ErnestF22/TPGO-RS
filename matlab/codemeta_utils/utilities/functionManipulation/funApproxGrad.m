function df=approx_grad(f,X,h)
if(exist('h','var')==0)
    h=1e-6;
end
sz=size(X);
for(r=1:sz(1))
    for(c=1:sz(2))
        v=zeros(sz);
        v(r,c)=1;
        f1=@(t) f(X+t*v);
        df(r,c)=approx_der(f1,0,h);
    end
end
