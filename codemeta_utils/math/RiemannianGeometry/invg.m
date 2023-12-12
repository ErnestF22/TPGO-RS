%function g1=invg(g)
%Compute the inverse of a rigid body transformation represented as a 4x4
%matrix

%%AUTORIGHTS%%

function g1=invg(g)
g1=g;
d=size(g,2)-1;
R=g(1:d,1:d,:);
for ig=1:size(g,3)
    g1(1:d,1:d,ig)=R(:,:,ig)';
    g1(1:d,d+1,ig)=-R(:,:,ig)'*g(1:d,d+1,ig);
end
