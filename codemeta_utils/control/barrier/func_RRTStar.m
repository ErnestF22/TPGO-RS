% RRT* algorithm in 2D with collision avoidance.
% 
% Author: Sai Vemprala
% 
% nodes:    Contains list of all explored nodes. Each node contains its
%           coordinates, cost to reach and its parent.
% 
% Brief description of algorithm: 
% 1. Pick a random node q_rand.
% 2. Find the closest node q_near from explored nodes to branch out from, towards
%    q_rand.
% 3. Steer from q_near towards q_rand: interpolate if node is too far away, reach
%    q_new. Check that obstacle is not hit.
% 4. Update cost of reaching q_new from q_near, treat it as Cmin. For now,
%    q_near acts as the parent node of q_new.
% 5. From the list of 'visited' nodes, check for nearest neighbors with a 
%    given radius, insert in a list q_nearest.
% 6. In all members of q_nearest, check if q_new can be reached from a
%    different parent node with cost lower than Cmin, and without colliding
%    with the obstacle. Select the node that results in the least cost and 
%    update the parent of q_new.
% 7. Add q_new to node list.
% 8. Continue until maximum number of nodes is reached or goal is hit.

function [ path ] = func_RRTStar(start, goal, boundary, obstacles, EPS, numNodes) 
% INPUTS:
%   start   := a start point (x,y) within the boundary
%   end     := an end point (x,y) within the boundary
%   boundary  := sets the limit on the workspace [xmin,xmax,ymin,ymax]
%   obstacles := a list of circular obstacles (bounded by rectangles). The
%       list is in the form [xLeft,yBottom,width,height], each row is an
%       obstacle
%   EPS       := max step size
%   numNodes  := the max number of nodes to generate
% OUTPUTS:
%   path      := a list of waypoints from q_start to q_goal in boundary

h = figure('units','normalized','outerposition',[0 0 1 1]);
filename = 'RRT_Star.gif';
xlabel('x (m)','FontSize',20);
ylabel('y (m)','FontSize',20);

x_min = boundary(1);
x_max = boundary(2);
y_min = boundary(3);
y_max = boundary(4);
% obstacles = [5,10,10,10;...
%     30,30,10,10]; % rectangle uses bottom left coord to locate, last 2 parameters are the diameters
% EPS = 2;
% numNodes = 1000;

q_start.coord = start;
q_start.cost = 0;
q_start.parent = 0;
q_goal.coord = goal;
q_goal.cost = 0;

nodes(1) = q_start;
axis([x_min x_max y_min y_max])
for i = 1:size(obstacles,1)
    obs = obstacles(i,:);
    rectangle('Position',obs,'FaceColor',[0 .5 .5],'Curvature',[1 1])
    hold on
end
% Plot the start and end locations
plot(start(1), start(2), 'rx', 'MarkerSize',10);
hold on
plot(goal(1), goal(2), 'gx', 'MarkerSize', 10);

% Start recording gif
frame = getframe(h);
im = frame2im(frame);
[imind,cm] = rgb2ind(im,256);
imwrite(imind,cm,filename,'gif','Loopcount',inf);

for i = 1:1:numNodes
    q_rand = [floor((x_max-x_min)*rand(1)+x_min), ...
        floor((y_max-y_min)*rand(1)+y_min)];
    plot(q_rand(1), q_rand(2), 'x', 'Color',  [0 0.4470 0.7410])
    
    % Break if goal node is already reached
    for j = 1:1:length(nodes)
        if nodes(j).coord == q_goal.coord
            break
        end
    end
    
    % Pick the closest node from existing list to branch out from
    ndist = [];
    for j = 1:1:length(nodes)
        n = nodes(j);
        tmp = dist(n.coord, q_rand);
        ndist = [ndist tmp];
    end
    [val, idx] = min(ndist);
    q_near = nodes(idx);
    
    q_new.coord = steer(q_rand, q_near.coord, val, EPS);
    if noCollision(q_rand, q_near.coord, obstacles)
        line([q_near.coord(1), q_new.coord(1)], [q_near.coord(2), q_new.coord(2)], 'Color', 'k', 'LineWidth', 2);
        drawnow
        hold on
        q_new.cost = dist(q_new.coord, q_near.coord) + q_near.cost;
        
        % Within a radius of r, find all existing nodes
        q_nearest = [];
        r = 60;
        neighbor_count = 1;
        for j = 1:1:length(nodes)
            if noCollision(nodes(j).coord, q_new.coord, obstacles) && dist(nodes(j).coord, q_new.coord) <= r
                q_nearest(neighbor_count).coord = nodes(j).coord;
                q_nearest(neighbor_count).cost = nodes(j).cost;
                neighbor_count = neighbor_count+1;
            end
        end
        
        % Initialize cost to currently known value
        q_min = q_near;
        C_min = q_new.cost;
        
        % Iterate through all nearest neighbors to find alternate lower
        % cost paths
        
        for k = 1:1:length(q_nearest)
            if noCollision(q_nearest(k).coord, q_new.coord, obstacles) && q_nearest(k).cost + dist(q_nearest(k).coord, q_new.coord) < C_min
                q_min = q_nearest(k);
                C_min = q_nearest(k).cost + dist(q_nearest(k).coord, q_new.coord);
%                 line([q_min.coord(1), q_new.coord(1)], [q_min.coord(2), q_new.coord(2)], 'Color', 'g');                
%                 hold on
            end
        end
        
        % Update parent to least cost-from node
        for j = 1:1:length(nodes)
            if nodes(j).coord == q_min.coord
                q_new.parent = j;
            end
        end
        
        % Append to nodes
        nodes = [nodes q_new];
    end
    % Capture a frame
    if (mod(i, floor(numNodes/10)) == 0)
        frame = getframe(h);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        imwrite(imind,cm,filename,'gif','WriteMode','append');
    end
end

D = [];
for j = 1:1:length(nodes)
    tmpdist = dist(nodes(j).coord, q_goal.coord);
    D = [D tmpdist];
end

% Search backwards from goal to start to find the optimal least cost path
[~, idx] = min(D);
% q_final = nodes(idx);
q_goal.parent = idx;
q_end = q_goal;
nodes = [nodes q_goal];
path = q_end.coord;
while q_end.parent ~= 0
    start = q_end.parent;
    line([q_end.coord(1), nodes(start).coord(1)], [q_end.coord(2), nodes(start).coord(2)], 'Color', 'r', 'LineWidth', 2);
    hold on
    q_end = nodes(start);
    path = [path; q_end.coord];
    % Capture a frame
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    imwrite(imind,cm,filename,'gif','WriteMode','append');
end

path = flipud(path);
end


