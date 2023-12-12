function[analyticalBounds, critResult] = POC_FeasibleRegionBeta0(kd,maxD,magW)
% Test if the 3 limiting inequalities determines the min Kv for a given kd
% (assuming beta = 0)
% INPUTS:
%   kd := scalar gain to find min kv
%   maxD := max scalar distance (0,pi)
%   magW := max speed (0,4)


% Define the system parameters
syms m2Var m3Var kdVar kvVar thetaVar magWVar real
assumeAlso(kdVar >0)
assumeAlso(thetaVar>0 & thetaVar < 1)
assumeAlso(magWVar > 0 & magWVar < 4)
Lambda = [1;maxD/2*cot(maxD/2)];

fun1 = m2Var*(-kdVar*thetaVar-kvVar/2-magWVar/4)+m3Var*(-kdVar/2*thetaVar+magWVar^2/8+magWVar/4*kvVar)+1/2;
fun2 = m2Var*(1-kvVar/2-magWVar/4)+m3Var*(-kvVar-kdVar/2*thetaVar+magWVar^2/8+magWVar/4*kvVar)+1/2;
fun3 = m2Var*(1+kvVar/2-magWVar/4)+m3Var*(-kvVar+kdVar/2+magWVar^2/8+magWVar/4*kvVar)-1/2;

% Solve for m2,m3,kv
critResult = solve([fun1;fun2;fun3]==0,[m2Var m3Var kvVar]);
m2 = vpa(subs(critResult.m2Var,[kdVar thetaVar magWVar],[kd maxD/2*cot(maxD/2) magW]));
m3 = vpa(subs(critResult.m3Var,[kdVar thetaVar magWVar],[kd maxD/2*cot(maxD/2) magW]));
kv = vpa(subs(critResult.kvVar,[kdVar thetaVar magWVar],[kd maxD/2*cot(maxD/2) magW]));

% Define the 8 inequalities that must be satisfied
coeff = @(kd,kv,Lambda,magW) [ ...
    [-kd*Lambda-kv/2-magW/4,...
    -kd*Lambda/2+magW^2/8+magW/4*kv,...
    +1/2*[1;1] ];...
    % Gerg Disc 3 and 4
    [-kd*Lambda+kv/2-magW/4,...
    kd*Lambda/2+magW^2/8+magW/4*kv,...
    -1/2*[1;1] ];...
    % Gerg Disc 5 and 6
    [1-kv/2-magW/4*[1;1],...
    -kv-kd*Lambda/2+magW^2/8+magW/4*kv,...
    +1/2*[1;1] ];...
    % Gerg Disc 7 and 8
    [1+kv/2-magW/4*[1;1],...
    -kv+kd*Lambda/2+magW^2/8+magW/4*kv,...
    -1/2*[1;1] ]];

mTermsMat = [m2;m3;1];

% Display Analytical Results

analyticalBounds=simplify(coeff(kdVar,critResult.kvVar,[1;thetaVar],magWVar)*[critResult.m2Var;critResult.m3Var;1]);
pretty(analyticalBounds)

% Display Numerical results
fprintf('kv = %0.5f, m2 = %0.5f, m3 = %0.5f\nResults below should all be less than 0:\n',kv,m2,m3);

coeff(kd,kv,Lambda,magW)*mTermsMat
end