function spconvmtx_test
c=(1:5)';
x=(1:10)';

y = conv(x,c);
n = length(x);
Cc=convmtx(c,n);
Ccsp=spconvmtx(c,n);
disp(norm(Cc-Ccsp,'fro'))

yB = Ccsp*x;
disp(norm(y-yB))

ysame=conv(x,c,'same');
Ccsp=spconvmtx(c,n,'same');
ysameB = Ccsp*x;
disp(norm(ysame-ysameB))

yvalid=conv(x,c,'valid');
Ccsp=spconvmtx(c,n,'valid');
yvalidB = Ccsp*x;
disp(norm(yvalid-yvalidB))
