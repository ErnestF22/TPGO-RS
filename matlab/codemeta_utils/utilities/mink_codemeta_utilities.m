%Find top values of columns by repeatedly applying min
function [DBest,idx]=mink_codemeta_utilities(D,NBest)
[m,n]=size(D);
DBest=zeros(NBest,n);
idx=zeros(NBest,n);
for iBest=1:NBest
    [DBest(iBest,:),idx(iBest,:)]=min(D);
    D(sub2ind([m,n],idx(iBest,:),1:n))=Inf;
end
