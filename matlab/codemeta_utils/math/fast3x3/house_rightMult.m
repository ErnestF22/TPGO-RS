function A=house_rightMult(A,v,beta)
%Compute A*P, where P is the Householder transformation given by the pair
%v, beta
A=A-(beta*A*v)*v';
