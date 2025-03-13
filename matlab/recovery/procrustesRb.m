function RbEst=procrustesRb(c,q)
[U,~,V]=svd(c*q');
diagmat = eye(size(U,1));
diagmat(end, end) = det(U*V');
RbEst=U*diagmat*V';
end
