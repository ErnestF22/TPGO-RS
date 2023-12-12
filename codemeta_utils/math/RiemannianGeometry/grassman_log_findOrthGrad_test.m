clear variables
randn('state',0)
n=7;
p=3;

e=stiefel_eye(randn(n,p));
h1=grassman_tangentProj(e,randn(n,p));
y1=grassman_exp(e,h1);
h12=grassman_tangentProj(y1,1*randn(n,p));
y2=grassman_exp(y1,h12);

Y1=y1;
Y2=y2;

% hest=grassman_log(Y1,Y2);
% Hest=[hest [hest(p+1:end,:)'; zeros(n-p)]];

YY1=orthCompleteBasis(Y1);
YY2=orthCompleteBasis(Y2);

[YY2,errors]=grassman_log_findOrthGrad(YY1,YY2,p);
[YY2,errorsLin]=grassman_log_findOrthGrad(YY1,YY2,p,'linearInit');

E=[eye(p); zeros(n-p,p)];
v=YY1'*rot_log(YY1,YY2);
disp(v);

