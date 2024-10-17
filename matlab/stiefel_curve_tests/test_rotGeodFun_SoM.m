clear;
close all;
clc;

R_start = RandOrthMat(5);

[Rt,dRt,R0,dR0,vVec,ddRt,dvVec]=rot_geodFun(R_start,[],'speed','cubic');

figure(1)
funCheckDer(Rt,dRt,'angle')
% figure(2)
% funCompare(vVec, @(t) rot_vee(Rt(t),dRt(t)))
figure(3)
funCheckDer(dRt,ddRt,'angle')
% figure(4)
% funCheckDer(vVec,dvVec,'angle')

