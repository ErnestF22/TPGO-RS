function POCReviewR_check27
w=sym('w',[2 1]);
omega=sym('omega',[2 1]);
R=POCReviewRwz(w(1),w(2),0);
tR=R(1:2,1:2);
a=(omega.'*tR*omega)-omega.'*omega;
disp(simplify(a))


