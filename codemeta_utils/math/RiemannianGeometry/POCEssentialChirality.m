function POCEssentialChirality
load triangulate_test_dataset_datacalibrated.mat
e3=[0;0;1];
Rxpi=diag([1;-1;-1]);
Rypi=diag([-1;1;-1]);
Rzpi=diag([-1;-1;1]);

Q=essential_fromG(G(:,:,1),G(:,:,2),'poses');
R1=Q(1:3,:);
R2=Q(4:6,:);

E=essential_toE(Q);

x1=homogeneous(x(:,1,1),3);
x2=homogeneous(x(:,1,2),3);

A=[R1*x1 -R2*x2];

l1=A\e3;
disp('(I,I)')
disp(l1(1)*R1*x1-l1(2)*R2*x2-e3)

l2=-l1;
disp('(Rxpi,Rxpi)')
disp(l2(1)*Rxpi*R1*x1-l2(2)*Rxpi*R2*x2-e3)

l3=[l1(1); -l1(2)]/(l1(1)*R1(3,:)*x1+l1(2)*R2(3,:)*x2);
disp('(I,Rzpi)')
disp(l3(1)*R1*x1-l3(2)*Rzpi*R2*x2-e3)

l4=-l3;
disp('(Rxpi,Rypi)')
disp(l4(1)*Rxpi*R1*x1-l4(2)*Rypi*R2*x2-e3)

