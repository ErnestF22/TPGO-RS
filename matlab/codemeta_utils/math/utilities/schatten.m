function n=schatten(A,p)
s=svd(A);
n=norm(s,p);
