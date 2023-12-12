clear all
close all
A =  [2  2 -1 1 0 0 0 17;
      0  2 -1 0 1 0 0 13;
      2  0 -1 0 0 1 0 13;
      0  0 -1 0 0 0 1 -1;
      0  0  0 0 0 0 0  0];
basic = [4 5 6 7];
[basic,result,A]=dual_simplex(A,basic);
V = find_vertices(A,basic);
%AX+B X=[x1;x2;x3]
% A[x1;x2]+b=x3
%yellow plane: A=[2;2] b=-17 
% [2;2]*[x1;x2]-17=x3

% A'[x1;x2;x3]<=-b adding slack variables 
%A'[x1;x2;x3]+S1== -b
