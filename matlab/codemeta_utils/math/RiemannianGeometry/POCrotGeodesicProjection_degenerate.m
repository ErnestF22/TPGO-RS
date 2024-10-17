function POCrotGeodesicProjection_degenerate
ez=[0;0;1];
Rz=@(t) rot(t*ez);
R=rot(pi*cnormalize([randn(2,1);0]));
%R=rot(pi*cnormalize(randn(3,1)));

c1=R(1,1)+R(2,2);
c2=R(2,1)-R(1,2);

disp('[c1 c2]')
disp([c1 c2])
f=@(t) rot_dist(Rz(t),R);

funPlot(f,'angle')

