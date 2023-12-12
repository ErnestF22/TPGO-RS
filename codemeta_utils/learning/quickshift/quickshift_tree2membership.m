%Transform quickshift forest into labels
function membership=quickshift_tree2membership(treeVectorMembership)
%init memberships
membership=zeros(size(treeVectorMembership));

%find roots as those points for which the tree vector points to itself
flagRoots=treeVectorMembership==1:length(treeVectorMembership);

%set labels for roots
membership(flagRoots)=1:sum(flagRoots);

%propagate labels by copying label of target of tree vector pointer into
%the source
while any(membership==0)
    membership=membership(treeVectorMembership);
end
