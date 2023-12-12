n=2;
K=5;
x=10*randn(K,6)

y=veronese(x,n);
y=y+0.1*randn(size(y));

x1=inverseveronese(y,n,K)