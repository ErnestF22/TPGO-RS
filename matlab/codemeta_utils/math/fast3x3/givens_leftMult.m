function A=givens_leftMult(A,c,s)
%compute R*A, where R=[c s; -s c]
A=[c s; -s c]*A;
