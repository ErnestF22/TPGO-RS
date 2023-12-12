function [result] = plot_ContractionInequalities(kd,kv,magW,maxD,m2,m3)
% Plot the 8 inequality contraction constraints for beta=0, as a function
% of m2 (set m2^2=m3).

% INPUTS:
%   kd := scalar gain
%   kv := scalar gain
%   magW := max speed
%   maxD := max distance
%   m2 := list of m2's to test
% OUTPUTS:
%   result := the 8 entries evaluated over m2

if ~exist('m3','var')
    m3 = 1;
end

if ~exist('m2','var')
    m2 = linspace(1e-3,m3);
end

Lambda = [1;maxD/2*cot(maxD/2)];

A = m2.*(-kd*Lambda-kv/2-magW/4) + m3*(-kd*Lambda/2+magW^2/8+magW*kv/4) + 1/2;
B = m2.*(-kd*Lambda+kv/2-magW/4) + m3*(kd*Lambda/2+magW^2/8+magW*kv/4) - 1/2;
C = m2.*(1-kv/2-magW/4) + m3*(-kv-kd*Lambda/2+magW^2/8+magW*kv/4) + 1/2;
D = m2.*(1+kv/2-magW/4) + m3*(-kv+kd*Lambda/2+magW^2/8+magW*kv/4) - 1/2;

figure
subplot(2,4,1);
plot(m2,A(1,:),'bo');
subplot(2,4,2);
plot(m2,A(2,:),'bo');
subplot(2,4,3);
plot(m2,B(1,:),'bo');
subplot(2,4,4);
plot(m2,B(2,:),'bo');
subplot(2,4,5);
plot(m2,C(1,:),'bo');
subplot(2,4,6);
plot(m2,C(2,:),'bo');
subplot(2,4,7);
plot(m2,D(1,:),'bo');
subplot(2,4,8);
plot(m2,D(1,:),'bo');

result = [A;B;C;D];
end

