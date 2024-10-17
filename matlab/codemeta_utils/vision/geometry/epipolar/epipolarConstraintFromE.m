%Evaluate the epipolar constraint x1'*E*x2 given E
%function e=epipolarConstraintFromE(E,x1,x2)
function e=epipolarConstraintFromE(E,x1,x2)
x1=homogeneous(x1,3);
x2=homogeneous(x2,3);
e=sum(x1.*(E*x2));
