function plot_env(y,xe,L)

plot( [y(1,end) y(1,1)],[y(2,end) y(2,1)], '.-.r');
hold on
plot( [y(1,2) y(1,1)],[y(2,2) y(2,1)], '.-r','LineWidth', 2);
hold on
plot( [y(1,2) y(1,3)],[y(2,2) y(2,3)], '.-r','LineWidth', 2);
hold on
plot( [y(1,3) y(1,4)],[y(2,3) y(2,4)], '.-r','LineWidth', 2);
hold on


% plot(xe(1),xe(2),'c*')


for i=1:size(y,2)
    plot(y(1,i),y(2,i),'-bo','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize', 12)
    hold on
end

% for i=1:size(L,2)
%     plot(L(1,i),L(2,i),'-go','MarkerSize', 12)
%     hold on
% end

% grid on
axis equal

end