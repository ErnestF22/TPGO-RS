function sphere2_median_test

N=11;
sigma=1;

y=sphere_randn(eye(2,1),sigma,N);
z=zeros(size(y));
yMedian=sphere2_median(y);

quiver(z(1,:),z(2,:),y(1,:),y(2,:),1)
hold on
quiver(0,0,yMedian(1),yMedian(2),1,'r')
hold off
axis equal


