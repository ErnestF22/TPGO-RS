function POCReviewAlternativeFormInnerProductSphere

ui=sphere_randGeodFun();
uj=sphere_randGeodFun();
thetaij=@(t) sphere_dist(ui(t),uj(t));

funCompare(@(t) ui(t)'*uj(t), @(t) 2*cos(thetaij(t)/2)^2-1)




