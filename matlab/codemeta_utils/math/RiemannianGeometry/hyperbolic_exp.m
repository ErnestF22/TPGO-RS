function y=hyperbolic_exp(x,v)
flagPermute=false;
if(size(v,2)==1)
    v=permute(v,[1 3 2]);
    flagPermute=true;
end

[D,N]=size(v);

nv=sqrt(sum(v(1:end-1,:).^2)-v(end,:).^2);
p=sinh(nv)./nv;
p(nv==0)=0;
y=cosh(nv(ones(D,1),:)).*x(:,ones(1,N))+(ones(D,1)*p).*v;

if(flagPermute)
    y=permute(y,[1 3 2]);
end
