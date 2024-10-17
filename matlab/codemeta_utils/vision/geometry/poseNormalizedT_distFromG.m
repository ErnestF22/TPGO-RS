function d=poseNormalizedT_distFromG(G1,G2)
[R1,T1]=G2RT(G1);
[R2,T2]=G2RT(G2);

dR=rot_dist(R1,R2,'outputType','vector');
dT=computeAngles(T1,T2);
d=(dR+dT)';


function a=computeAngles(T1,T2)
NT=size(T1,2);
c=sum(T1.*T2);
s=zeros(1,NT);
for iT=1:NT
    s(iT)=norm(cross(T1(:,iT),T2(:,iT)));
end
a=atan2(s,c)';
