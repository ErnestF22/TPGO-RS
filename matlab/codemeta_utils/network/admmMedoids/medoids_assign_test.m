function medoids_assign_test
X=rand(2,200);
mu=rand(2,4);
%mu=rand(2,1)*[1 1];

idx=medoids_assign(X,mu);
plotLines(X,mu(:,idx),'k')
hold on
plotPoints(X,'x')
plotPoints(mu,'o')
hold off
axis equal
display(idx)
