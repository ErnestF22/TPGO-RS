function essential_disp_test

R1=rot_randn();
T1=cnormalize(-rand(3,1));
R2=rot_randn();
T2=cnormalize(rand(3,1));

subplot(2,1,1)
draw3dcameraFromAxesAndCenter(R1,T1)
hold on
draw3dcameraFromAxesAndCenter(R2,T2)
quiver3(T1(1),T1(2),T1(3),T2(1)-T1(1),T2(2)-T1(2),T2(3)-T1(3),0)
hold off
axis equal

Q=essential_build(R1,T1,R2,T2);
subplot(2,1,2)
essential_disp(Q)
axis equal