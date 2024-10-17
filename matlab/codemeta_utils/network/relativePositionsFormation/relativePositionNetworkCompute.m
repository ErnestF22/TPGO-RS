%function rij=relativePositionCompute(Ti,E)
%Compute the relative position measurements between nodes at locations Ti
%and along edges E
function rij=relativePositionCompute(Ti,E)
rij=Ti(:,E(:,2),:)-Ti(:,E(:,1),:);
