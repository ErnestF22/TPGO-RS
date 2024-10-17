function l3d=project2DLineTo3DPlane(R,T,l2d)
%pose interpretation
l3d=[R'*l2d;
     T'*l2d];
