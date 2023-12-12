function A = transferPDM(A,d)
Q = orth(randn(d));
D = diag(randi([1,5],[1,d]));
P = Q*D*Q';
t = randi([10,50], [d,1]);
% A(Px+t)<=b
A(:,1:d) = A(:,1:d)*P;
A(:,d+1) = A(:,d+1)+A(:,1:d)*t;
end
