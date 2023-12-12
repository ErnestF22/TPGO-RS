function goal_idx= search_best_goal_node(tree,expand_dis...
                                     ,goal,obstacle_list ...
                                     ,path_resolution)
goal_indx = [];
dist_to_goal_list = [];
for i = 1:size(tree,2)
    n = tree(i);
    d = calc_dist_to_goal(goal,n.position(1),n.position(2));
    dist_to_goal_list = [dist_to_goal_list d];
    if d <=expand_dis
        goal_indx = [goal_indx i];
    end
end
safe_goal_inds = [];
for i =1:length(goal_indx)
    t_node = steer(tree(goal_indx(i)).position,goal,expand_dis);
    if check_collision(tree(goal_indx(i)).position,t_node,obstacle_list,path_resolution)
        safe_goal_inds = [safe_goal_inds goal_indx(i)];
    end
end
cost = [];
node = [];
if ~isempty(safe_goal_inds)
    for i=1:length(safe_goal_inds)
        c = tree(safe_goal_inds(i)).cost;
        cost = [cost c];
    end
end
[~,idx] = min(cost);
goal_idx = safe_goal_inds(idx);
% node = steer(tree(safe_goal_inds(idx)).position,goal,expand_dis);
end
