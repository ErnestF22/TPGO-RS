function [R,Q1]=POCRQ2(A)
P=blkdiag(1,fliplr(eye(2)));
[R,Q]=POCRQ(A*P);
Q1=Q*P';
