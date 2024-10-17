function lambdas=epipolarTriangulateDepths(R,T,x1,x2)

x1=homogeneous(x1,3);
x2=homogeneous(x2,3);

nPoints=size(x1,2);
lambdas = zeros(2,nPoints);
for pt=1:nPoints
    lambdas(:,pt) = [x1(:,pt) -R*x2(:,pt)]\T;
end
