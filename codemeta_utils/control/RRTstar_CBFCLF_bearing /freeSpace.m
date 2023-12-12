function [free_area] = freeSpace(lim_x,lim_y,obstacles)

free_area=(lim_x(2)-lim_x(1))*(lim_y(2)-lim_y(1));
for i=1:size(obstacles,2)
    free_area=free_area-pi*(obstacles(3,i))^2;
end