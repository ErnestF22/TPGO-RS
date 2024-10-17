function X=triang_getX(x)
l=triang_getDepths(x);
X=([x;1 1]*l-[0;0;1])/2;
