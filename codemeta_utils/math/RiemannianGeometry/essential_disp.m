function essential_disp(Q)
if ~exist('Q','var')
    Q=zeros(6,3);
    Q(1:3,:)=rot(3*pi/8*[0;1;0]);
    Q(4:6,:)=rot(5*pi/8*[0;1;0]);
end

bLength=3;

[R1,T1,R2,T2]=essential_getRTPair(Q,'references');
T2=bLength*T2;

flagHold=ishold;
draw3dcameraFromRT(R1,T1,'references')
hold on
draw3dcameraFromRT(R2,T2,'references')
quiver3(T1(1),T1(2),T1(3),T2(1),T2(2),T2(3),0)
if ~flagHold
    hold off
end
xlabel('x')
ylabel('y')
zlabel('z')
axis square
axis equal