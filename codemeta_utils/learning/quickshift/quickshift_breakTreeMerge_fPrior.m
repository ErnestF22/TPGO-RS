%Verify whether to combine or not two clusters using only the membership
%prior criterion
%Used as argument to quickshift_breakTreeMerge
function [flag, componentDataMerged]=quickshift_breakTreeMerge_fPrior(cmp1,cmp2,edge1,edge2)
%if a component is empty, initialize with the corresponding edge value
if isempty(cmp1)
    cmp1=edge1;
end
if isempty(cmp2)
    cmp2=edge2;
end

%flag=isempty(intersect(cmp1,cmp2));
flag=~overlapInt(cmp1,cmp2);
if flag
    componentDataMerged=[cmp1 cmp2];
else
    componentDataMerged=NaN;
end

%Return true if the intersection between two sets is non-empty
%The two sets are assumed to contain only positive integers
%Essentially, check the sum between two sparse indicator vectors to look
%for overlaps.
function flag=overlapInt(cmp1,cmp2)
l=max([cmp1 cmp2]);
l1=length(cmp1);
l2=length(cmp2);
u1=ones(1,l1);
u2=ones(1,l2);
flag=any((sparse(u1,cmp1,u1,1,l)+sparse(u2,cmp2,u2,1,l))>1);
