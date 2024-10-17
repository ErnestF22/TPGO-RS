function v=hyperbolic_log(x,y)
flagPermute=false;
if(size(y,2)==1)
    y=permute(y,[1 3 2]);
    flagPermute=true;
end

[D,N]=size(y);

%d=hyperbolic_lorentzMetric(x,y);
d=sum(x(1:end-1,ones(1,N)).*y(1:end-1,:))-x(end).*y(end,:);
v=acosh(-d(ones(D,1),:)).*(y+d(ones(D,1),:).*x(:,ones(1,N)))./(ones(D,1)*sqrt(d.^2-1));

if(flagPermute)
    v=permute(v,[1 3 2]);
end
