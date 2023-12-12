function POCessential_distMinAnglePair_relations
resetRands(3)
e3=[0;0;1];
Q1=cat(3,eye(3),eye(3));
Q2=rot_randn([],[],2);
Rzt=@(t) rot(t*e3);

Rp1=zeros(3,3,4);
Rp2=zeros(3,3,4);
for k=1:4
    Rp1(:,:,k)=essential_flipAmbiguity_R1(eye(3),k);
    Rp2(:,:,k)=essential_flipAmbiguity_R2(eye(3),k);
end

[tMin,fMin,~,output]=essential_distMinAngle([eye(3);eye(3)],[Q2(:,:,1);Q2(:,:,2)]);

k=output.idxMin;
t=tMin;

Q211t=zeros(3,3,4);
Q212t=zeros(3,3,4);
for k=1:4
    t=output.tMin(k);
    Q211t(:,:,k)=Rzt(t)*Rp1(:,:,k)*Q2(:,:,1);
    Q212t(:,:,k)=Rzt(t)*Rp2(:,:,k)*Q2(:,:,2);
    dft=e3'*(logrot(Q211t(:,:,k))+logrot(Q212t(:,:,k)));
    disp(dft)
end

disp([Rp1 Rp2])

