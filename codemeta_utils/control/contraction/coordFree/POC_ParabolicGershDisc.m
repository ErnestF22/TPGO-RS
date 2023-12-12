% Check the derivatives of the diagonalized gershgorin disc of the omega
% matrix, in particular the disc with parabolic function for the
% max encompassing value

close all; clear all; clc;

% Define the variables (kv >= 0, m1*m3-m2^2 > 0, normW >= 0)
syms normW m1 m2 m3 kv real

% Define the parabolic disc
alpha = normW*(m3*kv-m2)/4;
Dw_D2 = -m2/4*normW^2 + sqrt( (-m3/8*normW^2)^2 + alpha^2);

% Find the first derivative
der_Dw_D2 = diff(Dw_D2,normW);
% evaluated at 0 from the "right"
% Can ignore from the "left" since normW is nonnegative!
der_Dw_D2_0right = limit(der_Dw_D2,normW,0,'right')