function l=triang_getDepths(x)
l=[x(:,1) -x(:,2); 1 -1]\[0;0;-1];