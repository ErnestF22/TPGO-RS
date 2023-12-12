function [] = plotCLF(x,t,xd)
% Plot the value of the CLF function over a trajectory
% INPUTS:
%   x := state vector [3 x numTimeSteps]
%   t := time vector [numTimeSteps x 1]
%   xd := the desired point [3x1] vector, assumed to be fixed
% OUTPUTS:

figure;
hold on
% Compute CLF
CLF = zeros(size(t));
for i = 1:length(t)    
    CLF(i) = compute_CLF(x(i,:),xd);
end
plot(t,CLF);
title('CLF');
end

