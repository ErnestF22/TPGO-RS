function stiefel_log_test
randn('state',0)

disp '3,3'
test(3,3)
disp '3,2'
test(3,2)
disp '3,1'
test(3,1)

disp '8,8'
test(8,8)
disp '8,7'
test(8,7)
disp '8,4'
test(8,4)
disp '8,3'
test(8,3)
disp '8,1'
test(8,1)

disp '9,9'
test(9,9)
disp '9,7'
test(9,7)
disp '9,4'
test(9,4)
disp '9,3'
test(9,3)
disp '9,1'
test(9,1)

%Y1=orth(randn(n,p));

function test(n,p)
t=0.5*pi;
Y1=[eye(p); zeros(n-p,p)];
H12=stiefel_tangentProj(Y1,randn(n,p));
H12=H12/sqrt(stiefel_metric(Y1,H12,H12))*t;
[Y2]=stiefel_exp(Y1,H12);
disp([orthErr(Y2) sqrt(sum(sum((H12-stiefel_log(Y1,Y2)).^2))) ...
    sqrt(sum(sum((H12-stiefel_log(Y1,Y2,'method','SO(n)')).^2)))])
