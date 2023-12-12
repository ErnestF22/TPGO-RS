function POCessentialGeodesicDisplay
flagRecordVideo=false;
tMax=2*pi;
NFrames=80;

resetRands();
switch 4
    case 1
        Q1=essential_eye();
        v1=essential_tangentProj(Q1,essential_hat(Q1,[0;0;1;0;0;0]));
    case 2
        v=[pi/2;0;0];
        Q1=[rot(v);rot(v)];
        v1=essential_tangentProj(Q1,essential_hat(Q1,[0;1;0;0;0;0]));
    case 3
        v=[pi/2;0;0];
        Q1=[rot(v);rot(v)];
        v1=essential_tangentProj(Q1,essential_hat(Q1,[1;0;0;0;0;0]));
    case 4
        v=[pi/2;0;0];
        Q1=[rot(v);rot(v)];
        v1=essential_tangentProj(Q1,essential_hat(Q1,[1;0;0;-1;0;0]));
end
Q2t=essential_geodFun(Q1,v1);
%Q2t=essential_randGeodFun(Q1);
Q2rt=@(t) essential_closestRepresentative(Q1,Q2t(t));
Q2mt=@(t) inverseRotations(Q1,Q2t(t));
Q2ft=@(t) flipSecondRotation(Q2t(t));
Q2mft=@(t) flipSecondRotation(Q2mt(t));
Q2fmt=@(t) inverseRotations(Q1,Q2ft(t));

if flagRecordVideo
    vWrt=VideoWriter([mfilename '.avi']);
    vWrt.FrameRate=5;
    open(vWrt);
end

figure(1)
for t=linspace(0,tMax,NFrames)
    subplot(2,3,1)
    essential_disp(Q2t(t))
    axis equal
    axis([-1 5 -2 2 -2 2])
    title('Geodesic')
    subplot(2,3,2)
    essential_disp(Q2rt(t))
    axis equal
    axis([-1 5 -2 2 -2 2])
    title('Closest representative')
    subplot(2,3,3)
    essential_disp(Q2mt(t))
    axis equal
    axis([-1 5 -2 2 -2 2])
    title('Inverse geodesic')
    subplot(2,3,4)
    essential_disp(Q2ft(t))
    axis equal
    axis([-1 5 -2 2 -2 2])
    title('Geodesic with second rotation flipped')
    subplot(2,3,5)
    essential_disp(Q2mft(t))
    axis equal
    axis([-1 5 -2 2 -2 2])
    title('Inverse geodesic with second rotation flipped')
    subplot(2,3,6)
    essential_disp(Q2mft(t))
    axis equal
    axis([-1 5 -2 2 -2 2])
    title('Inverse of geodesic with second rotation flipped')
    pause(0.05)
    if flagRecordVideo
        rect=get(gcf,'Position');
        rect(3)=1024;
        rect(4)=768;
        set(gcf,'Position',rect);
        writeVideo(vWrt,getframe(gcf));
    end
end

function Q2m=inverseRotations(Q1,Q2)
Q2m(1:3,1:3)=rot_exp(Q1(1:3,:),-rot_log(Q1(1:3,:),Q2(1:3,:)));
Q2m(4:6,1:3)=rot_exp(Q1(4:6,:),-rot_log(Q1(4:6,:),Q2(4:6,:)));

function Q2f=flipSecondRotation(Q2)
Q2f=Q2;
Q2f(4:6,1:3)=diag([-1 -1 1])*Q2(4:6,1:3);


