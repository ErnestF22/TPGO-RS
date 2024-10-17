%Invert rotation and translation as a rigid body motion
%function [R,T]=invRT(R,T)
function [R,T]=invRT(R,T)
N=size(R,3);
for iN=1:N
    T(:,iN)=-R(:,:,iN)'*T(:,iN);
    R(:,:,iN)=R(:,:,iN)';
end

