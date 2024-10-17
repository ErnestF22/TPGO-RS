function POCepipolarRank
resetRands()

ex=eye(3,1);
R1=rot_randn();
R2=rot_randn();
G1=RT2G(R1,zeros(3,1));
G2=RT2G(R2,ex);

NX=10;
X=homogeneous([randn(2,NX); 10+randn(1,NX)],4);
P=[eye(3) zeros(3,1)];

x=projectFromRT(cat(3,R1,R2),[zeros(3,1) ex],X,'references');
l=projectGetDepthsFromRT(cat(3,R1,R2),[zeros(3,1) ex],X,'references');

xhom=homogeneous(x,3);

xhom1=xhom(:,:,1);
xhom2=xhom(:,:,2);
l1=ones(3,1)*l(1,:);
l2=ones(3,1)*l(2,:);
T=ex*ones(1,NX);
%l1.*(hat(ex)*R1*xhom1)-l2.*(hat(ex)*R2*xhom2)
%[(R1(2,:)*xhom1)./(R1(3,:)*xhom1);(R2(2,:)*xhom2)./(R2(3,:)*xhom2)]
%rank([homogeneous(R1*xhom1,2);homogeneous(R2*xhom2,2)]) 
[homogeneous(R1*xhom1,2);homogeneous(R2*xhom2,2)] 



%disp(rank([R1'*xhom(:,:,1); R2'*xhom(:,:,2)]))

