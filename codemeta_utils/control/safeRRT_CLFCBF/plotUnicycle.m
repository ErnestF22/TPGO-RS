function [] = plotUnicycle(t,x,varargin)
% Plot the center of mass of the unicycle as a point mass in 2D with a
% vector indicating the orientation
% INPUTS:
%   t := simulation time vector (a column vector)
%   x := System state vectors [x,y,theta] as a column vector
% OUTPUTS:
%   

% Define Parameters
stepSize = 100;

ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'stepsize'
            ivarargin = ivarargin+1;
            stepSize = varargin{ivarargin};
        otherwise
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

% Plot a state plot with stepSize
figure
hold on
for i = 1:stepSize:length(t)
    xCurrent = x(i,:);
    x1 = xCurrent(1); y1 = xCurrent(2); theta = xCurrent(3);
    xVector = [cos(theta);sin(theta)];
    plot(x1,y1,'ro','MarkerSize',5)
end
quiver(x1,y1,xVector(1),xVector(2),0,'r'); % plot last orientation
title('Phase Portrait');
end

