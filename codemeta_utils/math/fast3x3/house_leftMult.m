function A=house_leftMult(A,v,beta)
%Compute P*A, where P is the Householder transformation given by the pair
%v, beta
A=A-(beta*v)*(v'*A);
