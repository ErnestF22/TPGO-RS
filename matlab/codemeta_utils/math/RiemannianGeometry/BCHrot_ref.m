%function [w,alpha,beta,gamma]=BCHrot_ref(u1,u2)
%Implement the closed form solution of the BCH formula for rotations
function [w,alpha,beta,gamma]=BCHrot_ref(u1,u2)
N=size(u1,2);
w=zeros(3,N);
for(n=1:N)
    u11=u1(:,n);
    u21=u2(:,n);
    w(:,n)=logrot(rot(u11)*rot(u21));
    abc(:,n)=[u11 u21 cross(u11,u21)]\w(:,n);
end

alpha=abc(1,:);
beta=abc(2,:);
gamma=abc(3,:);
