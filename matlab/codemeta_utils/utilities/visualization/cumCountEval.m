%Evaluate cumulative counts
%function [c,x]=cumCountEval(data,x)
%Unnormalized version of cumDistEval
function [c,x]=cumCountEval(data,x)
data=shiftdim(data);
NX=length(x);

c=zeros(1,NX);
for iX=1:NX
    c(iX)=sum(data<=x(iX));
end
