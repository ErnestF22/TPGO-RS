load triangulate_test_dataset_datacalibrated
x1=[x(:,2:3,2);1 1];
x2=[x(:,2:3,1);1 1];
R12=R(:,:,2);
T12=T(:,2);
E=hat(T12)*R12;
%x1(:,1)'*E*x2(:,1)
e3hat=hat([0;0;1]);
Rzpihalf=rot([0;0;pi/2]);
[U,S,V]=svd(E);
R1=(U*Rzpihalf')';
R2=V';
x1prime=R1*x1;
x2prime=R2*x2;
disp('Point epipolar constraint')
disp(sum((x1prime).*(e3hat*x2prime)))

l1=null(x1');
l2=null(x2');

l1prime=R1*l1;
l2prime=R2*l2;

disp('Line check in local coordinates')
disp([l1'*x1 l2'*x2]);
disp('Line check in world coordinates')
disp([l1prime'*x1prime l2prime'*x2prime]);

plothyperplanes([l1prime l2prime],'patchOptions',{'FaceAlpha',0.5});
hold on
plot3dvect(zeros(3,1),[0;0;1],'z')
plot3dvect(zeros(3,1),l1prime,'n1')
plot3dvect(zeros(3,1),l2prime,'n2')
plot3dvect(zeros(3,1),x1prime(:,1),'x11')
plot3dvect(zeros(3,1),x1prime(:,2),'x12')
plot3dvect(zeros(3,1),x2prime(:,1),'x21')
plot3dvect(zeros(3,1),x2prime(:,2),'x22')
hold off
axis equal
