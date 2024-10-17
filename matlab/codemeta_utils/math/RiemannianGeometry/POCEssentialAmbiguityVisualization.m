function POCEssentialAmbiguityVisualization
N=100;
Qeye=essential_eye();
Q0=essential_exp(Qeye,essential_hat(Qeye,reshape(pi/2*cnormalize(ones(3,2)),[],1)));
vzVec=[0;0;1]*linspace(-2*pi,2*pi,N);
vQzVec=[vzVec;vzVec];
vQz=essential_hat(Qeye,vQzVec);

I=eye(3);
Rxpi=diag([1;-1;-1]);
Rypi=diag([-1;1;-1]);
Rzpi=diag([-1;-1;1]);

Ha=repmat(cat(3,Rxpi,Rxpi),[1,1,N]);
Hb=repmat(cat(3,I,Rzpi),[1,1,N]);
Hab=repmat(cat(3,Rxpi,Rypi),[1,1,N]);

Q=essential_exp(Q0,vQz);
Qa=essential_sharp(multiprod(essential_flat(Q),Ha));
Qb=essential_sharp(multiprod(essential_flat(Q),Hb));
Qab=essential_sharp(multiprod(essential_flat(Q),Hab));

vQVec=essential_vee(Qeye,doublerot_log(Q));
vQaVec=essential_vee(Qeye,doublerot_log(Qa));
vQbVec=essential_vee(Qeye,doublerot_log(Qb));
vQabVec=essential_vee(Qeye,doublerot_log(Qab));

[x,y,z] = sphere(2^5-1);
x=pi*x;
y=pi*y;
z=pi*z;

offset=0.1*[0;0;-1;0;0;0]*ones(1,N);
vQbVec=vQbVec+offset;
vQabVec=vQabVec+offset;


subplot(1,2,1)
plotPoints(vQVec(1:3,:),'r.');
hold on
plotPoints(vQaVec(1:3,:),'g.');
plotPoints(vQbVec(1:3,:),'b.');
plotPoints(vQabVec(1:3,:),'k.');
surf(x,y,z,'FaceColor','None','EdgeColor',1-[0.1 0.1 0.1]);
hold off
axis equal
view(30,30)

subplot(1,2,2)
plotPoints(vQVec(4:6,:),'r.');
hold on
plotPoints(vQaVec(4:6,:),'g.');
plotPoints(vQbVec(4:6,:),'b.');
plotPoints(vQabVec(4:6,:),'k.');
surf(x,y,z,(1-[0.1 0.1 0.1])*ones(size(x)),'FaceColor','None');
hold off
axis equal
view(30,30)

%keyboard


function v=doublerot_log(Q)
v=essential_sharp(rot_hat(eye(3),logrot(essential_flat(Q))));
