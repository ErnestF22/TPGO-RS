function [Q,tRand]=essential_randomVerticalMotion(Q)
tRand=2*pi*rand-pi;
R=Rz(tRand);
Q(1:3,:)=R*Q(1:3,:);
Q(4:6,:)=R*Q(4:6,:);

function R=Rz(t)
R=[cos(t) -sin(t) 0; sin(t) cos(t) 0; 0 0 1];
