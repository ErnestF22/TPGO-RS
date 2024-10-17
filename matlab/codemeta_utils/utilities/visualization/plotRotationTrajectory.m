function plotRotationTrajectory(R,style)
if ~exist('style','var')
    style='-';
end
plotPoints(squeeze(R(:,1,:)),['r' style])
hold on
plotPoints(squeeze(R(:,2,:)),['g' style])
plotPoints(squeeze(R(:,3,:)),['b' style])
plot3dframe(zeros(3,1),R(:,:,1))
hold off
axis equal
axis([-1 1 -1 1 -1 1])
