function rnd = get_random_node(min_env_x,max_enx_x,min_env_y,max_enx_y)
x = min_env_x + (max_enx_x-min_env_x) .* rand(1,1);
y = min_env_y + (max_enx_y-min_env_y) .* rand(1,1);
rnd = [x;y];
end