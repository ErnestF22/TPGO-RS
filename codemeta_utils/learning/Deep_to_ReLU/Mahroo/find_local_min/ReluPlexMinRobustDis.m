% given matrix of the eqaution that make the region, returns the min 
% distance to move p and remain in the same region:
% A(:,1)*x+A(:,2)*y+A(:,3) = z.
function dRobust = ReluPlexMinRobustDis(A,p)
d = zeros(size(A,1),1);
for i=1:size(A,1)
    a = A(i,1);
    b = A(i,2);
    c = A(i,3);
    d(i)=PointToLineDistance(a,b,c,p); 
end
dRobust = min(d);
end