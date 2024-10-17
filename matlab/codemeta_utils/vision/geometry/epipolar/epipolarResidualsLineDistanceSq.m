function e=epipolarResidualsLineDistanceSq(E,x1,x2)
Nx=size(x1,2);
x1=[x1;ones(1,Nx)];
x2=[x2;ones(1,Nx)];
e3hat=[0 -1 0; 1 0 0; 0 0 0];

numerator=(sum(x1.*(E*x2))).^2;
denominator1=sum((e3hat*E*x2).^2);
denominator2=sum((e3hat*E'*x1).^2);
e=numerator./denominator1+numerator./denominator2;
