function data=trifocal_dataset(NLines)
if ~exist('NLines','var')
    NLines=6;
end

R1=rot([pi/2;0;0]);
R2=rot([0;pi/8;0])*rot([pi/8;0;0])*R1*rot_randn(eye(3),0.1);
R3=rot([0;-pi/8;0])*rot([-pi/8;0;0])*R1*rot_randn(eye(3),0.1);
T1=-R1*[0;0;0];
T2=-R2*[2;0;0];
T3=-R3*[-2;2;0];

%l3d=[1 0 0 0; 0 1 0 -4]';
c=[0;4;0];
l3d=zeros(4,2,NLines);
for iLine=1:NLines
    ci=c+0.5*randn(3,1);
    b=sphere_randn();
    bOrth=householderRotation(b)*[zeros(1,2); eye(2)];
    l3d(:,:,iLine)=[bOrth;(-bOrth'*ci)'];
end

l2d1=project3DLineTo2DLine(R1,T1,l3d);
l2d2=project3DLineTo2DLine(R2,T2,l3d);
l2d3=project3DLineTo2DLine(R3,T3,l3d);

p3d1=project2DLineTo3DPlane(R1,T1,l2d1);
p3d2=project2DLineTo3DPlane(R2,T2,l2d2);
p3d3=project2DLineTo3DPlane(R3,T3,l2d3);

data.R=cat(3,R1,R2,R3);
data.T=[T1 T2 T3];
data.G=RT2G(data.R,data.T);
data.l2d=cat(3,l2d1,l2d2,l2d3);
data.p3d=cat(3,p3d1,p3d2,p3d3);
data.l3d=l3d;
