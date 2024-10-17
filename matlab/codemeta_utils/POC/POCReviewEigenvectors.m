function POCReviewEigenvectors
b=randn(3,1);
T=RT2G(eye(3),b);

disp(T*T'-(eye(4)+[b*b' b; b' 0]))

[U,S,V]=svd(T*T')
keyboard
