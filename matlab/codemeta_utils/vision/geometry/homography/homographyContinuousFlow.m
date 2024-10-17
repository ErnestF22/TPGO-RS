function dx=homographyContinuousFlow(H,x,lambda,dlambda)
xHom=homogeneous(x,3);
dx=H*xHom-([1;1;1]*(dlambda./lambda)).*xHom;
dx=dx(1:2,:);


