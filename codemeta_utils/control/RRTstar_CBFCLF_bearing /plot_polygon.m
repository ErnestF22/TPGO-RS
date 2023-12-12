function plot_polygon(obs1,obs_offset)
xobs = obs_offset(1,:);
yobs = obs_offset(2,:);
pgon = polyshape(xobs,yobs);
plot(pgon)
hold on 
xobs = obs1(1,:);
yobs = obs1(2,:);
pgon = polyshape(xobs,yobs);
plot(pgon)
hold on 
end