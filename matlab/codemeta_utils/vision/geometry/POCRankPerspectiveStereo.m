function POCRankPerspectiveStereo
resetRands()

R=rot_randn(eye(3),0.2);
T=randn(3,1);
G=RT2G(R,T);

B=[0;0;1];

Gl1=RT2G(eye(3),zeros(3,1));
Gr1=RT2G(eye(3),B);
Gl2=Gl1*G;
Gr2=Gr1*G;

NPoints=10;
X=randn(3,2)*randn(2,NPoints)+[0;0;-10]*ones(1,NPoints);
xl1=projectFromG(Gl1,X);
xr1=projectFromG(Gr1,X);
xl2=projectFromG(Gl2,X);
xr2=projectFromG(Gr2,X);

W=[xl1;xr1;xl2;xr2];

disp(rank(W))
