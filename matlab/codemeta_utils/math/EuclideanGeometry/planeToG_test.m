function planeToG_test
NPoints=30;
%methodAbsolutePoses='pose';
methodAbsolutePoses='reference';

N=sphere_randn([0;0;1],1);
d=5;

NVec=planeNdToNVec(N,d);
X=planeGeneratePoints(NVec,NPoints);

GPlane=planeToG(NVec,'methodAbsolutePoses',methodAbsolutePoses);

figure(1)
draw3dcameraFromG(eye(4))
hold on
draw3dcameraFromG(GPlane,'shape','plane','methodAbsolutePoses',methodAbsolutePoses)
plotPoints(X)
hold off
axis equal
