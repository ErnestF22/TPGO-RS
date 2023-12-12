%given the goal point and X=[x;y], this function computes the distance
%between goal point and X. 
function d = calc_dist_to_goal(goal,x,y)
d = norm([x;y]-goal,2);
end