function [y1,h1]=sphere_exp_fast(y,h)
flagPermute=false;
if(size(h,2)==1)
    h=permute(h,[1 3 2]);
    flagPermute=true;
end

[D,N]=size(h);

theta=sqrt(sum(h.^2));

if(nargout>1)
    h1=-y*sin(theta)*theta+h*cos(theta);
end

y1=y;
flagtheta=(theta>1e-12);

if(sum(flagtheta)>0)
    theta=theta(flagtheta);
    y1(:,flagtheta)=(y*ones(1,sum(flagtheta))).*(ones(D,1)*cos(theta)) ...
        + h(:,flagtheta).*(ones(D,1)*(sin(theta)./theta));
end

if(flagPermute)
    y1=permute(y1,[1 3 2]);
    if(nargout>1)
        h1=permute(h1,[1 3 2]);
    end
end
