%function gi1=clampLogRot(w)
%Let w be a rotation axis times the rotation agle. If norm(w)>pi, this
%function fixes w such that norm(w)<=pi but ti represents the same
%rotations 
function w1=clampLogRot(w)
normw=sqrt(sum(w.^2));
fixIdx=find(normw>pi);

normw1=pi-mod(normw(fixIdx),pi);

f=normw(fixIdx)./normw1;

w1=w;
if(~isempty(f))
    w1(:,fixIdx)=-w(:,fixIdx)./repmat(f,3,1);
end
