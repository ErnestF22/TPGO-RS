
function [path,fatal_samples,tree,rands] = expanding_tree(tree,goal,max_itr,lim_x,lim_y,expand_dis,...
    path_resolution,obstacle_list)

fatal_samples = [];%borders_of_environment(lim_x,lim_y,obstacle_list);
rands = [];

for i=1:max_itr
%     fprintf('Iter: %d, number of nodes:%d\n',i,size(tree,2));
    x_rnd = get_random_node(lim_x(1),lim_x(2),lim_y(1),lim_y(2));
    rands = [rands x_rnd];
    flag_inside_obstacle = inside_obstacle(obstacle_list,x_rnd);
    
    if flag_inside_obstacle
        fatal_samples = [fatal_samples x_rnd];
    end
    plot_random(x_rnd,fatal_samples);
    
    nearest_idx = get_nearest_node_index(tree,x_rnd);
    x_nearest = tree(nearest_idx).position;
    x_new = steer(x_nearest,x_rnd,expand_dis);
    
    if check_collision(x_nearest,x_new,obstacle_list,path_resolution)
        r = optimal_disc(tree,obstacle_list,lim_x,lim_y,expand_dis);
        near_indices = find_near_nodes(tree,x_new,r);
        tree(end+1).position = x_new;
        idx_min = nearest_idx;
        c_min = tree(idx_min).cost + norm(x_new-x_nearest);
        [idx_min,c_min] = choose_parent(x_new,idx_min,c_min,near_indices,obstacle_list,tree,path_resolution);
        tree(end).parent = idx_min;
        tree(end).cost = c_min;
        
        tree = rewire(tree,x_new,near_indices,obstacle_list,path_resolution);
        
        
    end
end
% fprintf('reach max iteration\n')

last_index = search_best_goal_node(tree,expand_dis...
    ,goal,obstacle_list...
    ,path_resolution);
path = [];
% tree = modify_tree(tree,obstacle_list,path_resolution);
if ~isempty(last_index)
    path = generate_path(tree,last_index,goal);

end
end