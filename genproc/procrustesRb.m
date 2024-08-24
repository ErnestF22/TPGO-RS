function RbEst=procrustesRb(c,q)
[U,~,V]=svd(c*q');
RbEst=U*diag([1 det(U*V')])*V';
end
