function [yMedian,fMedian]=sphere2_median(y)
if size(y,2)==1
    y=permute(y,[1 3 2]);
end
N=size(y,2);

f=zeros(1,N);
for iN=1:N
    yi=y(:,iN);
    f(iN)=sum(abs(sphere2_distSigned(y,yi)));
end

[fMedian,idxMedian]=min(f);
yMedian=y(:,idxMedian);
