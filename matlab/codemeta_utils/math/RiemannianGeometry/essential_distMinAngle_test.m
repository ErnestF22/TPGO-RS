function essential_distMinAngle_test
resetRands(3)
flagDegenerateCase=true;

e3=[0;0;1];
Q1=rot_randn([],[],2);
if flagDegenerateCase
    Q1b=[Q1(:,:,1);Q1(:,:,2)];
    Q2b=essential_randomVerticalMotion(Q1b);
    Q2=cat(3,Q2b(1:3,:),Q2b(4:6,:));
else
    Q2=rot_randn([],[],2);
end
Rzt=@(t) rot(t*e3);

Q21tk=@(t,k) Rzt(t)*essential_flipAmbiguity_R1(Q2(:,:,1),k);
Q22tk=@(t,k) Rzt(t)*essential_flipAmbiguity_R2(Q2(:,:,2),k);

figure(1)
[tMin,fMin,~,output]=essential_distMinAngle([Q1(:,:,1);Q1(:,:,2)],[Q2(:,:,1);Q2(:,:,2)]);
tMin=modAngle(tMin);
output.tMin=modAngle(output.tMin);
output.tBreak1=modAngle(output.tBreak1);
output.tBreak2=modAngle(output.tBreak2);
for k=1:4
    ft=@(t) (rot_dist(Q1(:,:,1),Q21tk(t,k))^2+rot_dist(Q1(:,:,2),Q22tk(t,k))^2);
    dft=@(t) 2*e3'*(Q1(:,:,1)*logrot(Q1(:,:,1)'*Q21tk(t,k))+Q1(:,:,2)*logrot(Q1(:,:,2)'*Q22tk(t,k)));
    check_der(ft,dft,'angle')
    hold on
    plot(output.tBreak1(k),ft(output.tBreak1(k)),'r+')
    plot(output.tBreak2(k),ft(output.tBreak2(k)),'g+')
end

plot(output.tMin,output.fMin,'ko')
plot(tMin,fMin,'kx','MarkerSize',20)

hold off

k=output.idxMin;
t=tMin;
Q11=Q1(:,:,1);
Q12=Q1(:,:,2);
Q21t=Q21tk(t,k);
Q22t=Q22tk(t,k);
Q211t=Q21t*Q11';
Q212t=Q22t*Q12';
ft=(rot_dist(eye(3),Q211t)^2+rot_dist(eye(3),Q212t)^2);
dft=2*e3'*(logrot(Q211t)+logrot(Q212t));
[~,ft2]=essential_distMinAnglePair_base(Q211t,Q212t);

disp([fMin ft ft2 fMin-ft fMin-ft2])
disp(dft)
