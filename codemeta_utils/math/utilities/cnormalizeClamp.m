function x=cnormalizeClamp(x,threshold)
normx=sqrt(sum(x.^2));
flagClamp=normx>threshold;
x(:,flagClamp)=x(:,flagClamp)./(repmat(normx(flagClamp),size(x,1),1));
