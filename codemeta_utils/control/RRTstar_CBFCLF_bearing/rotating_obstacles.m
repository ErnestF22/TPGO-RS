function obs2 = rotating_obstacles(obs,c,t)
obs2 = obs;
for i=1:size(obs,2)
    x = (obs(1,i)-c(1))*cosd(t)-(obs(2,i)-c(2))*sind(t)+c(1);
    y = (obs(1,i)-c(1))*sind(t)+(obs(2,i)-c(2))*cosd(t)+c(2);
    obs2(1,i) = x;
    obs2(2,i) = y;
end
end