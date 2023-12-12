function flag_inside_obstacle = inside_obstacle(obstacle_list,rnd)
flag_inside_obstacle = 0;
for i=1:size(obstacle_list,2)
    d = norm(rnd-[obstacle_list(1,i);obstacle_list(2,i)],2);
    if d<=obstacle_list(3,i)
        flag_inside_obstacle = 1;
    end
end
end