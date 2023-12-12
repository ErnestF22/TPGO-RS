function plot_exp(T,obs1,obs2,l1,l2,l3,l4)

% plot_tree_(T)
% hold on
plot_polygon (obs1,obs2)
hold on
plot_landmarks(l1,'blue','s')
hold on
plot_landmarks(l2,'blue','s')
hold on
plot_landmarks(l3,'blue','s')
hold on
plot_landmarks(l4,'blue','s')
hold on
end
% % start = [220;120];
% start = [236;130];
% k_path = navigate_exp(start,T,1,'r*',l1,l2,l3,l4);
