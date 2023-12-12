syms l1 l2 l3 a b c d e f
A = [a+b+c d+e+f];
B = [l1*a+l2*b+l3*c l1*d+l2*e+l3*f];
aa = A/norm(A)
bb = B/norm(B)