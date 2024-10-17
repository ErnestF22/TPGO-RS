function [A,U,V]=bidiagonalize3x3(A)
%zeros on first column
[u3,bu3]=house(A(:,1));
A=house_leftMult(A,u3,bu3);
%zero on first row
[v2,bv2]=house(A(1,2:3)');
A(:,2:3)=house_rightMult(A(:,2:3),v2,bv2);
%zero on second column
[u2,bu2]=house(A(2:3,2));
A(2:3,:)=house_leftMult(A(2:3,:),u2,bu2);

%Accumulate U and V
U=blkdiag(1,eye(2)-bu2*(u2*u2'));
U=house_rightMult(U,u3,bu3);
V=blkdiag(1,eye(2)-bv2*(v2*v2'));
