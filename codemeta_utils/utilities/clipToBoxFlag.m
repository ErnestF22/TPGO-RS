%Returns flag indicating if each point is in the box or not
function flag=clipToBoxFlag(x,a)
dimX=size(x,1);
if length(a)~=dimX
    a=a*ones(1,dimX);
end
m=max(diag(a)\abs(x),[],1);
flag=m<1;
