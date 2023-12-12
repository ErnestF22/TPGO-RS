function hi=sphere_log_fast(yi,yj)
flagPermute=false;

D=size(yi,1);
if(size(yj,2)==1)
    yj=permute(yj,[1 3 2]);
    flagPermute=true;
end
hi=yj;
p=min(1,max(-1,yi'*yj));
hflag=abs(acos(p))>1e-12;
hi(:,~hflag)=zeros(D,sum(~hflag));
if(sum(hflag)>0)
    hi(:,hflag)=(ones(D,1)*(acos(p(hflag))./sqrt(1-p(hflag).^2))).*sphere_tangentProj(yi,yj(:,hflag));
end
if(flagPermute)
    hi=permute(hi,[1 3 2]);
end
