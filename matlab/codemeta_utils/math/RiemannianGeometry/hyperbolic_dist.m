function d=hyperbolic_dist(x,y)
if(size(y,2)==1)
    y=permute(y,[1 3 2]);
end

[~,N]=size(y);

v=hyperbolic_log(x,y);
d=sqrt(sum(v.^2))';

%d=acosh(-hyperbolic_lorentzMetric(x,y));
%d=acosh(-sum(x(1:end-1,ones(1,N)).*y(1:end-1,:))+x(end).*y(end,:));
