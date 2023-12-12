function xProj=triang_project(x)
[U,S,V]=svd(x);
xProj=U(:,1)*S(1,1)*V(:,1)';
