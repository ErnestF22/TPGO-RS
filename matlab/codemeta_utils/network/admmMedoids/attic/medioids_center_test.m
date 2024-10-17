function medioids_center_test
resetRands()
x=rand(2,40);
mu=medioids_center(x);
plotPoints(x)
hold on
plotPoints(mu,'o')
hold off
