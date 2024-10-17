function POCTrifocal

R1=rot([pi/2;0;0]);
R2=rot([0;pi/8;0])*rot([pi/8;0;0])*R1;
R3=rot([0;-pi/8;0])*rot([-pi/8;0;0])*R1;
T1=-R1*[0;0;0];
T2=-R2*[2;0;0];
T3=-R3*[-2;2;0];

l3d=[1 0 0 0; 0 1 0 -4]';

l2d1=project3DLineTo2DLine(R1,T1,l3d);
l2d2=project3DLineTo2DLine(R2,T2,l3d);
l2d3=project3DLineTo2DLine(R3,T3,l3d);

R0=rot_randGeodFun();
R1p=@(t) R1*R0(t);
R2p=@(t) R2*R0(t);
R3p=@(t) R3*R0(t);

figure(1)
drawScene(R1,T1,R2,T2,R3,T3,l2d1,l2d2,l2d3,l3d)

figure(2)
subplot(1,3,1)
draw2dLine(l2d1)
axis([-1 1 -1 1])
subplot(1,3,2)
draw2dLine(l2d2)
axis([-1 1 -1 1])
subplot(1,3,3)
draw2dLine(l2d3)
axis([-1 1 -1 1])

% t=linspace(-1,1,100);
% t=[t fliplr(t)];
% Nt=length(t);
% 
% figure(3)
% for it=1:Nt
%     t1=t(it);
%     drawScene(R1p(t1),T1,R2p(t1),T2,R3p(t1),T3,l2d1,l2d2,l2d3,l3d)
%     axis([-8 8 -8 8 -8 8])
%     pause(0.01)
% end

Te=computeTrifocal(RT2G(cat(3,R1,R2,R3),cat(2,T1,T2,T3)));
display(Te)

function drawScene(R1,T1,R2,T2,R3,T3,l2d1,l2d2,l2d3,l3d)
optsDrawPlanes={'optsdraw3dplane',{'side',4}};
G=RT2G(cat(3,R1,R2,R3),[T1 T2 T3]);

%compute intersection of lines
l3d1=project2DLineTo3DPlane(R1,T1,l2d1);
l3d2=project2DLineTo3DPlane(R2,T2,l2d2);
l3d3=project2DLineTo3DPlane(R3,T3,l2d3);

[U,S,V]=svd([l3d1 l3d2 l3d3]);
% disp(['Next-to-last singular value for intersection of lines: ' num2str(S(3,3))])

testNetworkDisplay(G,'poses')
hold on
draw3dLine(l3d,'side',2)
draw3dCameraLineBackprojection(R1,T1,l2d1,optsDrawPlanes{:});
draw3dCameraLineBackprojection(R2,T2,l2d2,optsDrawPlanes{:});
draw3dCameraLineBackprojection(R3,T3,l2d3,optsDrawPlanes{:});


draw3dLine(U(:,1:2),'style','g','side',8)

hold off
axis square
axis equal


function Te=computeTrifocal(G)
%assumes T1=zeros(3,1)
R=G2R(G);
T=G2T(G);
E=eye(3);
Ra=R(:,:,2)*R(:,:,1)';
Rb=R(:,:,1)*R(:,:,3)';
Te=zeros(3,3,3);
for k=1:3
    Te(:,:,k)=Ra*E(:,k)*T(:,3)'-T(:,2)*E(:,k)'*Rb;
end
