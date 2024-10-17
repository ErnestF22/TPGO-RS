function A=givens_rightMult(A,c,s)
%compute A*R, where R=[c s; -s c]
A=A*[c s; -s c];
