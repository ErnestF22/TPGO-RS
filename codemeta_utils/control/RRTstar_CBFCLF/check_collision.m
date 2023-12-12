% Given the obstacles (center of circle and its radisu), if there is no 
% collision between the node and its parent, this function returns True. 
% This function divides the line between the node and the parent to a set
% of points with resolution distance from each other and check if all the 
% points are outside of the obstacles.   
function collisionFlag = check_collision(parent,node,obstacle_list,resolution)
path = Line(parent,node,resolution);
collisionFlag = 1;
for o=1:size(obstacle_list,2)
    dis = [];
    for i=1:size(path,2)
        d = norm(obstacle_list(1:2,o)-path(:,i));
        dis = [dis d];
    end
    min_dis = min(dis);
    if min_dis < obstacle_list(3,o)
        collisionFlag = 0;
    end
end
end
