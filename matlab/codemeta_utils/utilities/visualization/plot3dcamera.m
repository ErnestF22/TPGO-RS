%PLOT3DCAMERA(O,R)
%
%O  column vector for the origin of the camera frame
%R  the column of R are the axis vectors (they should be orthonormal)
function plot3dcamera(o,r)
%axes
plot3dframe(o,r)
%image plane
c1=o+r(:,3)-r(:,1)-r(:,2);
c2=o+r(:,3)+r(:,1)-r(:,2);
c3=o+r(:,3)+r(:,1)+r(:,2);
c4=o+r(:,3)-r(:,1)+r(:,2);
c=[c1 c2 c3 c4];
hold on
h=patch(c(1,:),c(2,:),c(3,:), 'b');
set(h,'FaceAlpha',0.5)
hold off