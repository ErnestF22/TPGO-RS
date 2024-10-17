function [R,Q]=POCRQ(A)
[Qt,Rt]=qr(A');
R=Rt';
Q=Qt';
