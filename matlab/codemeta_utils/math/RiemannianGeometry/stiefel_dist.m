function d=stiefel_dist(Y1,Y2)
H12=stiefel_log(Y1,Y2);
N=size(Y2,3);
d=zeros(N,1);
for iN=1:N
    d(iN)=sqrt(stiefel_metric(Y1,H12(:,:,iN),H12(:,:,iN)));
end
