clear all
close all
A=[0 -1 1 0 0 -3;1 0 0 1 0 15;-1 1 0 0 1 -5;0 0 0 0 0 0];
basic = [3 4 5];
[basic,result,A]=dual_simplex(A,basic);
V = find_vertices(A,basic);