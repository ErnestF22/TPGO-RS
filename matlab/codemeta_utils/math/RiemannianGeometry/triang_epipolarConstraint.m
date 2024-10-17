function v=triang_epipolarConstraint(x)
v=[x(1);x(2)].'*[0 -1; 1 0]*[x(3);x(4)];

