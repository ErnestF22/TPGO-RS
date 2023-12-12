function[analyticalBounds, critResult] = POC_FeasibleRegionBetaNonZero(kd,maxD,magW,beta)
% Test if the 3 limiting inequalities determines the min Kv for a given kd
% (assuming beta = 0)
% INPUTS:
%   kd := scalar gain to find min kv
%   maxD := max scalar distance (0,pi)
%   magW := max speed (0,4)
% OUTPUTS:
%   analyticalBounds := bounds on the 8 inequalities
%   critResult := m2,m3,kv value at critical point

% Define the system parameters
syms m2Var m3Var kdVar kvVar thetaVar magWVar betaVar real
assume(kdVar >0)
assume(thetaVar>0 & thetaVar < 1)
assume(magWVar > 0 & magWVar < 4)
assume(betaVar > 0 )
Lambda = [1;maxD/2*cot(maxD/2)];

fun1 = m2Var*(-kdVar*thetaVar+betaVar-kvVar/2-magWVar/4)+m3Var*(-kdVar/2*thetaVar+magWVar^2/8+magWVar/4*kvVar)+betaVar+1/2;
fun2 = m2Var*(1+betaVar-kvVar/2-magWVar/4)+m3Var*(betaVar-kvVar-kdVar/2*thetaVar+magWVar^2/8+magWVar/4*kvVar)+1/2;
fun3 = m2Var*(1-betaVar+kvVar/2-magWVar/4)+m3Var*(betaVar-kvVar+kdVar/2+magWVar^2/8+magWVar/4*kvVar)-1/2;

% Solve for m2,m3,kv
critResult = solve([fun1;fun2;fun3]==0,[m2Var m3Var kvVar]);
m2 = vpa(subs(critResult.m2Var,[kdVar thetaVar magWVar betaVar],[kd maxD/2*cot(maxD/2) magW beta]));
m3 = vpa(subs(critResult.m3Var,[kdVar thetaVar magWVar betaVar],[kd maxD/2*cot(maxD/2) magW beta]));
kv = vpa(subs(critResult.kvVar,[kdVar thetaVar magWVar betaVar],[kd maxD/2*cot(maxD/2) magW beta]));

% Define the 8 inequalities that must be satisfied
coeff = @(kd,kv,Lambda,magW,beta) [ ...
    [-kd*Lambda+beta-kv/2-magW/4,...
    -kd*Lambda/2+magW^2/8+magW/4*kv,...
    beta+1/2*[1;1] ];...
    % Gerg Disc 3 and 4
    [-kd*Lambda-beta+kv/2-magW/4,...
    kd*Lambda/2+magW^2/8+magW/4*kv,...
    beta-1/2*[1;1] ];...
    % Gerg Disc 5 and 6
    [1+beta-kv/2-magW/4*[1;1],...
    beta-kv-kd*Lambda/2+magW^2/8+magW/4*kv,...
    +1/2*[1;1] ];...
    % Gerg Disc 7 and 8
    [1-beta+kv/2-magW/4*[1;1],...
    beta-kv+kd*Lambda/2+magW^2/8+magW/4*kv,...
    -1/2*[1;1] ]];

mTermsMat = [m2';m3';[1,1]];

% Display Analytical Results
analyticalBounds=simplify(coeff(kdVar,critResult.kvVar,[1;thetaVar],magWVar,betaVar)*[critResult.m2Var';critResult.m3Var';[1,1]]);
% pretty(analyticalBounds)

% Display Numerical results
fprintf('kv = %0.5f, m2 = %0.5f, m3 = %0.5f\n\tor\n',kv(1),m2(1),m3(1));
fprintf('kv = %0.5f, m2 = %0.5f, m3 = %0.5f\nResults below should all be less than 0:\n',kv(2),m2(2),m3(2));

coeff(kd,kv,Lambda,magW,beta)*mTermsMat
end
