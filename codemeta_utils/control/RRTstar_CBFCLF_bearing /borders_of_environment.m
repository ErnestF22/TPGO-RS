%Given the boarders of the environment: lim_x = [min_x,max_x] and lim_y =
%[min_y,max_y], this function generate samples on boarders. 
function fatal_samples = borders_of_environment(lim_x,lim_y,b)
x = lim_x(1);
y = lim_y(1):2:lim_y(2);
x_d = [x*ones(1,size(y,2));y];

x = lim_x(2);
y = lim_y(1):2:lim_y(2);
x_u = [x*ones(1,size(y,2));y];

y = lim_y(1);
x = lim_x(1):2:lim_x(2);
y_l = [x;y*ones(1,size(x,2))];

y = lim_y(2);
x = lim_x(1):2:lim_x(2);
y_r = [x;y*ones(1,size(x,2))];

fatal_samples = [x_u x_d y_l y_r];
end