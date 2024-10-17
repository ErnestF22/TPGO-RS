function POCBilinearFactorization
x=sphere_randn();
y=randn(3,1);
A=x*y'+y*x';
Z=x*y';


F=@(A) trace(A'*A);
T=@(A,B) trace(A'*B);
disp(F(A-(Z+Z')))
disp(F(A)+F(Z+Z')-2*T(A,Z+Z'))
disp(F(A)+2*F(Z)+2*T(Z,Z')-2*T(A,Z)-2*T(A,Z'))
disp(F(A)+2*F(Z)+2*T(Z,Z')-2*T(A,Z)-2*T(A',Z))
disp(F(A)+2*F(Z)+2*T(Z,Z')-2*T(A+A',Z))
disp(F(A)+2*F(Z)+2*trace(Z)^2-2*T(A+A',Z))


%[U,S,V]=svd(A+A');
%disp([x U(:,1)])
%keyboard
