function [] = animateUnicycle(t,x,varargin)
% Plot the center of mass of the unicycle as a point mass in 2D with a
% vector indicating the orientation
% INPUTS:
%   t := simulation time vector (a column vector)
%   x := System state vectors [x,y,theta] as a column vector
% OUTPUTS:
%   

% Define Parameters
stepSize = 100;
windowBuffer = 1; % Unit buffer on axis
pauseTime = 0.001;
xd = zeros(3,1);

ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'stepsize'
            ivarargin = ivarargin+1;
            stepSize = varargin{ivarargin};
        case 'windowbuffer'
            ivarargin = ivarargin+1;
            windowBuffer = varargin{ivarargin};
        case 'pausetime'
            ivarargin = ivarargin+1;
            pauseTime = varargin{ivarargin};
        case 'target'
            ivarargin = ivarargin+1;
            xd = varargin{ivarargin};
        otherwise
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

% Plot a state plot with stepSize
h=figure;
% Determine axis window
xMin = min(x(:,1))-windowBuffer;
xMax = max(x(:,1))+windowBuffer;
yMin = min(x(:,2))-windowBuffer;
yMax = max(x(:,2))+windowBuffer;
for i = 1:stepSize:length(t)
    figure(h);
    xCurrent = x(i,:);
    x1 = xCurrent(1); y1 = xCurrent(2); theta = xCurrent(3);
    xVector = 0.5*[cos(theta);sin(theta)];
    % Plot initial state
    plot(x(1,1),x(1,2),'bx','MarkerSize',5);
    hold on;
    % Plot current state
    plot(x1,y1,'ro','MarkerSize',5); 
    quiver(x1,y1,xVector(1),xVector(2),0,'r');
    % Plot desired state
    plot(xd(1),xd(2),'go','MarkerSize',5);
    xDVec = 0.5*[cos(xd(3));sin(xd(3))];
    quiver(xd(1),xd(2),xDVec(1),xDVec(2),0,'g');
    % Plot trajectory up to time t(i)
    plot(x(1:i,1),x(1:i,2),'r--');
    axis([xMin xMax yMin yMax]);
    hold off;    
    pause(pauseTime);
end
end