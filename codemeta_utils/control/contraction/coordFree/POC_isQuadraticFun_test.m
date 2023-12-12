% POC_isQuadraticFun_test.m
% Code to test if some arbitrary function, f, is quadratic.
% Steps:
%   1) define "x" as a linear function of time so that
%       f'(x)=(df/dx)*(dx/dt) where dx/dt = 1 since linear
%   2) Check if f'(x) is affine, if so f must be quadratic
close all; clear all; clc;

% Define an "unknown" quadratic function
a = rand*10; b = rand*10; c = rand*10;
fQuad = @(x) a*x^2+b*x+c;
% Define some "unknown" non quadratic function
fNotQuad = @(x) x^7+x^3+x^2+31;

% Define the test points
x1 = randn; x2 = randn;
x1_t = @(t) x1+t; x2_t = @(t) x2+t;
lambda = rand;

% Check if fQuad is quadratic with affine test on the derivative
A_quad=computeDer(fQuad,@(t) lambda*x1_t(t) + (1-lambda)*x2_t(t));
B_quad=lambda*computeDer(fQuad,x1_t) + (1-lambda)*computeDer(fQuad,x2_t);
if abs(B_quad-A_quad)<1e-6
    fprintf('fQuad is quadratic\n')
else
    fprintf(2,'fQuad is NOT quadratic\n')
end

% Check if fNotQuad quadratic with affine test on the derivative
A_notQuad = computeDer(fNotQuad,@(t) lambda*x1_t(t) + (1-lambda)*x2_t(t));
B_notQuad = lambda*computeDer(fNotQuad,x1_t)...
    + (1-lambda)*computeDer(fNotQuad,x2_t);
if abs(B_notQuad-A_notQuad)<1e-6
    fprintf('fNotQuad is quadratic\n')
else
    fprintf(2,'fNotQuad is NOT quadratic\n')
end


function der = computeDer(f,x)
% The new "unknown" function that we want to do an affine check on
der = funApproxDer(@(t) f(x(t)),0);
end
