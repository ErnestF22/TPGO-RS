function quickshift_plotDensity(X,P)
[Xs,Ys,Zs]=interpGridScattered(X(1,:),X(2,:),P,100);

mesh(Xs,Ys,Zs)
axis tight; hold on
plot3(X(1,:),X(2,:),P,'.','MarkerSize',15) 
hold off