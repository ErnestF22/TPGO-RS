% Test Pos Def bounding from "When is Rotations Averaging Hard?" K. Wilson
% et al.

close all; clear all; clc

R = rot_randn;

axang = rotm2axang(R);
w = axang(1:3)';
theta = axang(4);
u = theta*cot(theta/2);

H = (2-u)*[w*w' -w*w';-w*w' w*w']...
    + u*[eye(3) -eye(3);-eye(3) eye(3)]...
    + theta*[zeros(3) -hat3(w);hat3(w) eye(3)];

A = u*[eye(3) -eye(3);-eye(3) eye(3)]...
    + theta*[zeros(3) -hat3(w);hat3(w) zeros(3)];

B = u*[eye(3) -eye(3);-eye(3) eye(3)]...
    - theta*[eye(3) zeros(3);zeros(3) eye(3)];

fprintf("theta: %0.3f\n", theta);
fprintf("eigs of H: %0.3f, %0.3f, %0.3f, %0.3f, %0.3f, %0.3f\n", eig(H))
fprintf("eigs of H-A: %0.3f, %0.3f, %0.3f, %0.3f, %0.3f, %0.3f\n", eig(H-A))
fprintf("eigs of H-B: %0.3f, %0.3f, %0.3f, %0.3f, %0.3f, %0.3f\n", eig(H-B))