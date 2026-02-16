function R=procrustes_R(X,Y)
nb_dim=size(X,1);
[U,~,V]=svd(Y*X');
R=U*diag([ones(nb_dim-1,1);det(U*V')])*V';
end %function