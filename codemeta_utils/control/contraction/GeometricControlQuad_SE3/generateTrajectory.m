function [trajectory]=generateTrajectory(t)
    trajectory.position=[0.4*t;0.4*sin(pi*t);0.6*cos(pi*t)];
    trajectory.velocity=[0.4;0.4*pi*cos(pi*t);-0.6*pi*sin(pi*t)];
    trajectory.acceleration=[0;-0.4*pi*pi*sin(pi*t);-0.6*pi*pi*cos(pi*t)];
    trajectory.bodyDir1=[cos(pi*t);sin(pi*t);0];
end