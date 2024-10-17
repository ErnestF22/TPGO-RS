%Gives a measure of how far the columns of Y are far from being orthonormal
function e=orthErr(Y)
e=norm(Y'*Y-eye(size(Y,2)),'fro');
