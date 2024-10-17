%function E = epipolarEstimateE8pt(x1,x2)
%Estimate E matrix given a set of pairs of matching *calibrated* points
function E = epipolarEstimateE8pt(x1,x2)

NPoints=size(x1,2);
x1 = homogeneous(x1,3);
x2 = homogeneous(x2,3);

A = [] ;
for iPoint = 1:NPoints
    A = [A; kron(x2(:,iPoint), x1(:,iPoint))'] ;
end

[U,S,V] = svd(A) ;
E = reshape(V(:,9),3,3) ;

% Project E on the space of essential matrices
[U,S,V] = svd(E);

%U = det(U)*U;
%V = det(V)*V;

E = U * diag([1 1 0]) * V';
