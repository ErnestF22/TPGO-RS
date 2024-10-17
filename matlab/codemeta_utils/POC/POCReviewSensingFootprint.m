function POCReviewSensingFootprint
R=1;
alpha=pi/4;

phi=@(x,y) atan2(y,x);
c1=@(x,y) R^2-x.^2-y.^2;
c2=@(x,y) alpha-phi(x,y);
c3=@(x,y) alpha+phi(x,y);

C1=@(x,y) max(0,c1(x,y));
C2=@(x,y) max(0,c2(x,y));
C3=@(x,y) max(0,c3(x,y));

B=@(x,y) 1./C1(x,y)+1./C2(x,y)+1./C3(x,y);

S=@(x,y) 1./B(x,y);

x=linspace(-1.5,1.5);
[X,Y]=meshgrid(x,x);
Z=S(X,Y);
surf(X,Y,Z)


