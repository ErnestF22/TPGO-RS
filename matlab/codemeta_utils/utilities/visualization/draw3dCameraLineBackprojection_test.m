function draw3dCameraLineBackprojection_test
resetRands();
R=rot_randn();
T=randn(3,1);
l2d=[0.5;1;0.5];
subplot(1,2,1)
draw2dLine(l2d)
axis([-1 1 -1 1])
% axis equal
% axis square
% grid on

subplot(1,2,2)
draw3dcameraFromPose(R,T)
hold on
draw3dCameraLineBackprojection(R,T,l2d)
hold off
axis equal
axis square
grid on
