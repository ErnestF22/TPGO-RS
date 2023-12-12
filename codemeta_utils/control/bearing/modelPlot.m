function modelPlot()
%After running the simulink model, the trajectories are saved to the
%workspace, and this function plots that data

x = evalin('base', 'simout');
data = x.data;
l = evalin('base', 'simout1');
land = l.data(:,:,1);
plot(data(:,1),data(:,2),data(:,3),data(:,4))
hold on
plot(land(1,:),land(2,:),'xk')
hold off

end

