%Return vector of indicator for leaves of the tree
%function flag=quickshift_findLeaves(vTree)
function flag=quickshift_findLeaves(vTree)
flag=true(size(vTree));
flag(vTree)=false;

