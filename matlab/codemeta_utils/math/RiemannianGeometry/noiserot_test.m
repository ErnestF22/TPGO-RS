clear all
N=500;
Rtruth=rot([0;0;1]); % true value on which to reach consensus
sigmanoise=0.1;

x=[1;0;0];
for irot=1:N
    Ri{irot}=noiserot(Rtruth,sigmanoise);
    if(norm(Ri{irot}'*Ri{irot}-eye(3),'fro')>1e-6)
        warning('Check the code!')
    end
    xrot(:,irot)=Ri{irot}*x;
end

plot3(xrot(1,:),xrot(2,:),xrot(3,:),'.', 0,0,0,'o')
grid on
axis equal