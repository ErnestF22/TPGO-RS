function [e,grade,hessOpe]=epipolarResidualsSampsonSq(E,x1,x2)
x1=homogeneous(x1,3);
x2=homogeneous(x2,3);

numerator=(sum(x1.*(E*x2))).^2;
denominator1=sum((E(1:2,:)*x2).^2);
denominator2=sum((E(:,1:2)'*x1).^2);
e=numerator./(denominator1+denominator2);
